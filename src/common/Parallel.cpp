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
  namespace
  {
    // Nesting guard: a ParallelDo issued from inside a chunk body — worker
    // thread or the caller thread running a chunk — runs inline serially
    // instead of re-entering the pool. Top-level orchestration (e.g. the
    // ToString subtree fan-out) relies on this: subtree workers call
    // Newton/NTT code that itself uses ParallelDo.
    thread_local int tl_chunkDepth = 0;

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
      using Clock = std::chrono::steady_clock;
      const auto deadline = Clock::now() + std::chrono::microseconds(BIGMATH_POOL_SPIN_US);
      for (;;)
      {
        for (int i = 0; i < 32; ++i)
        {
          if (done()) return true;
          CpuRelax();
        }
        if (Clock::now() >= deadline) return done();
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

        uint32_t gen;
        {
          std::lock_guard<std::mutex> lk(m);
          curBody = body;
          curCtx = ctx;
          curTotal = total;
          curChunks = numChunks;
          gen = (uint32_t)(generation.load(std::memory_order_relaxed) + 1);
          pending.store(numChunks, std::memory_order_relaxed);
          claim.store((uint64_t)gen << 32, std::memory_order_relaxed);
          generation.store(gen, std::memory_order_release);
        }
        cv.notify_all();

        ClaimAndRun(gen, body, ctx, total, numChunks);

        if (!SpinUntil([this] { return pending.load(std::memory_order_acquire) == 0; }))
        {
          std::unique_lock<std::mutex> lk(m);
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
          Int c = (Int)(uint32_t)word;
          if (c >= chunks) return;
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
            std::lock_guard<std::mutex> lk(m);
            doneCv.notify_one();
          }
          word = claim.load(std::memory_order_acquire);
        }
      }

      void WorkerLoop()
      {
        uint32_t lastSeen = 0;
        for (;;)
        {
          // Stay hot for a short while after the previous job, then park.
          SpinUntil([this, lastSeen] {
            return generation.load(std::memory_order_acquire) != lastSeen ||
                   shutdown.load(std::memory_order_acquire);
          });

          void (*body)(Int, Int, void *) = nullptr;
          void *ctx = nullptr;
          Int total = 0;
          Int chunks = 0;
          {
            std::unique_lock<std::mutex> lk(m);
            cv.wait(lk, [this, lastSeen] {
              return generation.load(std::memory_order_acquire) != lastSeen ||
                     shutdown.load(std::memory_order_acquire);
            });
            if (shutdown.load(std::memory_order_acquire)) return;
            lastSeen = generation.load(std::memory_order_acquire);
            body = curBody;
            ctx = curCtx;
            total = curTotal;
            chunks = curChunks;
          }

          ClaimAndRun(lastSeen, body, ctx, total, chunks);
        }
      }

      SizeT numThreads = 1;
      std::vector<std::thread> workers;

      std::mutex dispatchMutex; // held by the one caller that owns the slot
      std::mutex m;             // guards the job description below + both cvs
      std::condition_variable cv;
      std::condition_variable doneCv;
      std::atomic<bool> shutdown{false};

      // Job description: written under `m` by the slot owner, read under `m`
      // by workers together with `generation`.
      void (*curBody)(Int, Int, void *) = nullptr;
      void *curCtx = nullptr;
      Int curTotal = 0;
      Int curChunks = 0;

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
