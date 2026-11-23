/**
 * BigMath: Internal thread pool for parallel NTT.
 *
 * Compiled only when BIGMATH_USE_THREADS=1. Lazy-initialized on first call.
 * Workers poll briefly for the next job (BIGMATH_POOL_SPIN_US) before parking
 * on a condition variable; one caller owns the pool at a time and any other
 * concurrent caller runs its work inline.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/common/Parallel.h"

#if BIGMATH_USE_THREADS

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <thread>
#include <vector>

// How long an idle worker (and a caller waiting for its last chunks) keeps
// polling before it parks on a condition variable. A parked thread costs
// 20-30 us to wake on macOS — more than a whole small NTT — while one
// multiply issues its dispatches only tens of microseconds apart, so a short
// spin keeps the pool hot across the phases of one operation and across
// back-to-back operations. 0 disables spinning (park immediately).
#ifndef BIGMATH_POOL_SPIN_US
#define BIGMATH_POOL_SPIN_US 100
#endif

namespace BigMath
{
#ifdef BIGMATH_TEST_POOL_CLOCK
  // The isolated temperature test advances worker deadlines explicitly.
  std::chrono::steady_clock::time_point PoolSpinNow() noexcept;
#else
  static inline auto PoolSpinNow() noexcept { return std::chrono::steady_clock::now(); }
#endif

  namespace
  {
    // Nesting guard: a ParallelDo issued from inside a chunk body — worker
    // thread or the caller thread running a chunk — runs inline serially
    // instead of re-entering the pool. Top-level orchestration (e.g. the
    // ToString subtree fan-out) relies on this: subtree workers call
    // Newton/NTT code that itself uses ParallelDo.
    thread_local int tl_chunkDepth = 0;
    std::atomic<bool> poolStarted{false};

    inline void CpuRelax()
    {
#if defined(__aarch64__) || defined(_M_ARM64)
      __asm__ __volatile__("yield");
#elif defined(__x86_64__) || defined(__i386__)
      __builtin_ia32_pause();
#endif
    }

    // Poll `done()` for up to BIGMATH_POOL_SPIN_US. Returns true if it became
    // true while spinning. The clock is read every 32 polls, so one poll is
    // an atomic load plus a pause.
    template <typename Pred>
    inline bool SpinUntil(Pred &&done)
    {
#if BIGMATH_POOL_SPIN_US > 0
      const auto deadline = PoolSpinNow() + std::chrono::microseconds(BIGMATH_POOL_SPIN_US);
      for (;;)
      {
        for (int i = 0; i < 32; ++i)
        {
          if (done()) return true;
          CpuRelax();
        }
        if (PoolSpinNow() >= deadline) return done();
      }
#else
      return done();
#endif
    }

    class ThreadPool
    {
    public:
      ThreadPool()
      {
        unsigned hw = std::thread::hardware_concurrency();
        if (hw == 0) hw = 4;
        if (hw > BIGMATH_MAX_THREADS) hw = BIGMATH_MAX_THREADS;
        if (hw < 1) hw = 1;
        numThreads = hw;

        if (numThreads > 1)
        {
          workers.reserve(numThreads - 1);
          for (SizeT i = 1; i < numThreads; ++i)
            workers.emplace_back([this] { WorkerLoop(); });
        }
        poolStarted.store(true, std::memory_order_release);
      }

      ~ThreadPool()
      {
        {
          std::lock_guard<std::mutex> lk(m);
          shutdown.store(true, std::memory_order_release);
        }
        cv.notify_all();
        for (auto &t : workers) t.join();
      }

      SizeT NumThreads() const { return numThreads; }

      // Workers that are not parked on the condition variable (spinning or
      // running a slice). A dispatch issued while enough of them are awake is
      // picked up within nanoseconds instead of after a 20-30 us wake-up.
      bool Warm() const
      {
#if BIGMATH_POOL_SPIN_US > 0
        const int needed = (int)std::min<SizeT>(numThreads - 1, 5);
        return needed > 0 && awake.load(std::memory_order_relaxed) >= needed;
#else
        return false; // workers park immediately; every dispatch is a cold one
#endif
      }

      // Wake parked workers without giving them work. They spin for
      // BIGMATH_POOL_SPIN_US and park again if nothing arrives. Non-blocking,
      // and cheap for the caller: it wakes ONE worker (a broadcast to seven
      // sleepers costs the caller ~12 us on macOS), and that worker relays
      // the wake-up to the rest from its own thread.
      void Prewarm()
      {
#if BIGMATH_POOL_SPIN_US > 0
        if (workers.empty() || awake.load(std::memory_order_relaxed) >= (int)workers.size())
          return;
        {
          std::lock_guard<std::mutex> lk(m);
          poke.fetch_add(1, std::memory_order_release);
        }
        cv.notify_one();
#endif
      }

      // Run body over `numChunks` equal slices of [0, total). Blocks until
      // every slice has completed. The caller works too: it and the workers
      // claim slices from a shared counter, so a dispatch never waits for a
      // parked worker to wake up when the caller could do the slice itself.
      void RunChunks(Int numChunks, Int total, void (*body)(Int, Int, void *), void *ctx)
      {
        if (numChunks <= 1 || tl_chunkDepth > 0)
        {
          body(0, total, ctx);
          return;
        }
        if ((SizeT)numChunks > numThreads) numChunks = (Int)numThreads;

        // The pool has one work slot. A second user thread dispatching while
        // the slot is busy runs its work inline instead of queueing behind
        // (or, as before this lock existed, overwriting) the job in flight.
        std::unique_lock<std::mutex> owner(dispatchMutex, std::try_to_lock);
        if (!owner.owns_lock())
        {
          body(0, total, ctx);
          return;
        }

        // Publish the job without a lock, so that awake workers pick it up
        // without one either (seven workers reading the description under a
        // shared mutex serialize at every dispatch, and the losers sleep in
        // the kernel). Only the slot owner writes, in this order:
        //   1. close `claim` under the new generation tag — from here on no
        //      claim against the previous job can succeed;
        //   2. write the description;
        //   3. open `claim` at slice 0;
        //   4. publish `generation`.
        // A worker reads `generation`, then the description, then claims
        // with that generation. A description read before step 4 belongs to
        // an older generation, and step 1 has already made every claim with
        // that generation fail — so a torn or stale description is never
        // acted on.
        const uint32_t gen = generation.load(std::memory_order_relaxed) + 1;
        claim.store(((uint64_t)gen << 32) | 0xFFFFFFFFu, std::memory_order_seq_cst);
        curBody.store(body, std::memory_order_relaxed);
        curCtx.store(ctx, std::memory_order_relaxed);
        curTotal.store(total, std::memory_order_relaxed);
        curChunks.store(numChunks, std::memory_order_relaxed);
        pending.store(numChunks, std::memory_order_relaxed);
        claim.store((uint64_t)gen << 32, std::memory_order_release);
        generation.store(gen, std::memory_order_seq_cst);

        // Parked workers need a wake-up. A worker announces that it is about
        // to park by decrementing `awake` (under `m`) before its last look at
        // `generation`; with both sides sequentially consistent, either it
        // sees this job or we see it missing here. Taking `m` orders the
        // notify after its wait has begun.
        if (awake.load(std::memory_order_seq_cst) < (int)workers.size())
        {
          {
            std::lock_guard<std::mutex> lk(m);
          }
          cv.notify_all();
        }

        ClaimAndRun(gen, body, ctx, total, numChunks);

        if (!SpinUntil([this] { return pending.load(std::memory_order_acquire) == 0; }))
        {
          std::unique_lock<std::mutex> lk(doneMutex);
          doneCv.wait(lk, [this] { return pending.load(std::memory_order_acquire) == 0; });
        }
      }

    private:
      // Claim slices of job `gen` until none are left. `claim` packs
      // (generation, next slice index) into one word, so a claim can only
      // succeed against the job whose body/ctx this thread actually read: a
      // thread that arrives after the job finished sees a different
      // generation tag (or an exhausted index) and touches nothing.
      void ClaimAndRun(uint32_t gen, void (*body)(Int, Int, void *), void *ctx, Int total, Int chunks)
      {
        const Int chunkSize = (total + chunks - 1) / chunks;
        uint64_t word = claim.load(std::memory_order_acquire);
        for (;;)
        {
          if ((uint32_t)(word >> 32) != gen) return;
          const uint32_t index = (uint32_t)word; // 0xFFFFFFFF while the slot is being written
          if (index >= (uint32_t)chunks) return;
          const Int c = (Int)index;
          if (!claim.compare_exchange_weak(word, word + 1, std::memory_order_acq_rel, std::memory_order_acquire))
            continue;

          Int s = c * chunkSize;
          Int e = std::min(total, s + chunkSize);
          if (s < e)
          {
            ++tl_chunkDepth;
            body(s, e, ctx);
            --tl_chunkDepth;
          }
          if (pending.fetch_sub(1, std::memory_order_acq_rel) == 1)
          {
            // Last slice done: wake the caller if it has parked. Taking the
            // lock orders this notify after the caller's predicate check.
            std::lock_guard<std::mutex> lk(doneMutex);
            doneCv.notify_one();
          }
          word = claim.load(std::memory_order_acquire);
        }
      }

      void WorkerLoop()
      {
        uint32_t lastSeen = 0;
        uint32_t lastPoke = 0;
        awake.fetch_add(1, std::memory_order_seq_cst);
        for (;;)
        {
          // Stay hot for a short while after the previous job (or a prewarm
          // poke). A job that arrives during the spin is picked up without
          // touching `m`.
          const bool hot = SpinUntil([this, lastSeen] {
            return generation.load(std::memory_order_acquire) != lastSeen ||
                   shutdown.load(std::memory_order_acquire);
          });

          bool poked = false;
          if (!hot)
          {
            // Park. `awake` drops before the last look at `generation` (see
            // RunChunks), and stays dropped only while actually waiting.
            std::unique_lock<std::mutex> lk(m);
            awake.fetch_sub(1, std::memory_order_seq_cst);
            for (;;)
            {
              if (generation.load(std::memory_order_seq_cst) != lastSeen ||
                  shutdown.load(std::memory_order_acquire))
                break;
              if (poke.load(std::memory_order_acquire) != lastPoke)
              {
                poked = true;
                break;
              }
              cv.wait(lk);
            }
            awake.fetch_add(1, std::memory_order_seq_cst);
            lastPoke = poke.load(std::memory_order_acquire);
          }

          if (shutdown.load(std::memory_order_acquire)) return;

          const uint32_t gen = generation.load(std::memory_order_acquire);
          if (gen == lastSeen)
          {
            // Poked awake with no job: relay the wake-up to any workers
            // still parked, then go back to the spin phase.
            if (poked && awake.load(std::memory_order_relaxed) < (int)workers.size())
              cv.notify_all();
            continue;
          }
          lastSeen = gen;
          void (*body)(Int, Int, void *) = curBody.load(std::memory_order_relaxed);
          void *ctx = curCtx.load(std::memory_order_relaxed);
          const Int total = curTotal.load(std::memory_order_relaxed);
          const Int chunks = curChunks.load(std::memory_order_relaxed);
          ClaimAndRun(gen, body, ctx, total, chunks);
        }
      }

      SizeT numThreads = 1;
      std::vector<std::thread> workers;

      std::mutex dispatchMutex; // held by the one caller that owns the slot
      std::mutex m;             // parking only: guards waits on `cv`
      std::condition_variable cv;
      std::mutex doneMutex;     // the slot owner's wait for the last slice
      std::condition_variable doneCv;
      std::atomic<bool> shutdown{false};

      // Job description: written by the slot owner while `claim` is closed,
      // read lock-free by workers (protocol in RunChunks). Atomics so that a
      // racing read is merely stale, never undefined.
      std::atomic<void (*)(Int, Int, void *)> curBody{nullptr};
      std::atomic<void *> curCtx{nullptr};
      std::atomic<Int> curTotal{0};
      std::atomic<Int> curChunks{0};

      std::atomic<int> awake{0};           // workers not parked on `cv`
      std::atomic<uint32_t> poke{0};       // bumped by Prewarm(); wakes workers into their spin phase
      std::atomic<uint32_t> generation{0}; // bumped once per job
      std::atomic<uint64_t> claim{0};      // (generation << 32) | next slice
      std::atomic<Int> pending{0};         // slices of the current job not yet finished
    };

    ThreadPool &Pool()
    {
      static ThreadPool p;
      return p;
    }
  }

  SizeT ParallelNumThreads()
  {
    return Pool().NumThreads();
  }

  SizeT ParallelAvailableThreads(bool allowParked)
  {
    if (tl_chunkDepth > 0 || (!allowParked && !poolStarted.load(std::memory_order_acquire)))
      return 1;
    auto &pool = Pool();
    return allowParked || pool.Warm() ? pool.NumThreads() : 1;
  }

  bool ParallelWarm()
  {
    return Pool().Warm();
  }

  void ParallelPrewarm()
  {
    Pool().Prewarm();
  }

  namespace
  {
    // End time of this thread's last operation that took the serial
    // algorithm because the pool was cold (steady-clock ns; 0 = never).
    thread_local int64_t tl_lastColdDoneNs = 0;

    inline int64_t NowNs()
    {
      return std::chrono::duration_cast<std::chrono::nanoseconds>(
                 std::chrono::steady_clock::now().time_since_epoch())
          .count();
    }
  }

  bool ParallelWarmOrPrewarm()
  {
    ThreadPool &pool = Pool();
    if (pool.Warm())
      return true;
#if BIGMATH_POOL_SPIN_US > 0
    // Cold. If this thread finished another such operation moments ago, more
    // are probably coming: wake the pool now so the next one finds it warm.
    // An isolated operation wakes nothing.
    const int64_t last = tl_lastColdDoneNs;
    if (last != 0 && NowNs() - last < (int64_t)BIGMATH_POOL_SPIN_US * 1000)
      pool.Prewarm();
#endif
    return false;
  }

  void ParallelNoteColdDone()
  {
#if BIGMATH_POOL_SPIN_US > 0
    tl_lastColdDoneNs = NowNs();
#endif
  }

  SizeT ParallelMinSize()
  {
    // Tuned from sample(1) profile of div 200k/50k under threads (was 4096):
    // psynch_cvwait dominated at 42K samples vs ~5K for actual NTT work.
    // Root cause: NTT has log2(n) layers per transform; each layer does one
    // ParallelFor dispatch costing ~30µs. CRT does 6 transforms × ~13 layers
    // ≈ 78 dispatches per Multiply, totalling ~2.3ms overhead — exceeds the
    // 1.5ms actual work at small NTT sizes.
    //
    // Raising to 65536 limits layer-parallelism to very large NTTs where the
    // per-layer work amortizes the dispatch overhead. Preserves large-mul
    // wins (1M×1M had layers ≥ 262144 ops) while eliminating regression on
    // small-to-mid NTTs.
    return 65536;
  }

  void ParallelFor(Int total, void (*body)(Int start, Int end, void *ctx), void *ctx)
  {
    if (total <= 0) return;
    SizeT n = Pool().NumThreads();
    if (n <= 1 || (SizeT)total < ParallelMinSize())
    {
      body(0, total, ctx);
      return;
    }
    Pool().RunChunks((Int)n, total, body, ctx);
  }

  void ParallelDo(Int numTasks, void (*body)(Int start, Int end, void *ctx), void *ctx)
  {
    if (numTasks <= 0) return;
    SizeT n = Pool().NumThreads();
    if (n <= 1)
    {
      body(0, numTasks, ctx);
      return;
    }
    // Use min(numTasks, pool size) chunks so each thread gets ~1 task.
    Int chunks = (Int)std::min((SizeT)numTasks, n);
    Pool().RunChunks(chunks, numTasks, body, ctx);
  }
}

#endif
