/**
 * BigMath: Internal thread pool for parallel NTT.
 *
 * Compiled only when BIGMATH_USE_THREADS=1. Lazy-initialized on first call.
 * Workers spin-wait on a condition variable for low dispatch latency.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/common/Parallel.h"

#if BIGMATH_USE_THREADS

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace BigMath
{
  namespace
  {
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
            workers.emplace_back([this, i] { WorkerLoop((Int)i); });
        }
      }

      ~ThreadPool()
      {
        {
          std::lock_guard<std::mutex> lk(m);
          shutdown = true;
          generation++;
        }
        cv.notify_all();
        for (auto &t : workers) t.join();
      }

      SizeT NumThreads() const { return numThreads; }

      // Dispatch `numChunks` parallel calls of body(start, end). The caller
      // thread runs chunk 0; workers 1..numChunks-1 run the others. Blocks
      // until all chunks complete.
      void RunChunks(Int numChunks, Int total, void (*body)(Int, Int, void *), void *ctx)
      {
        if (numChunks <= 1)
        {
          body(0, total, ctx);
          return;
        }
        if ((SizeT)numChunks > numThreads) numChunks = (Int)numThreads;

        // Publish work.
        {
          std::lock_guard<std::mutex> lk(m);
          curBody = body;
          curCtx = ctx;
          curTotal = total;
          curChunks = numChunks;
          remaining = numChunks - 1; // caller does chunk 0
          generation++;
        }
        cv.notify_all();

        // Caller runs chunk 0 itself.
        Int chunkSize = (total + numChunks - 1) / numChunks;
        Int s0 = 0;
        Int e0 = std::min(total, chunkSize);
        body(s0, e0, ctx);

        // Wait for workers.
        std::unique_lock<std::mutex> lk(m);
        doneCv.wait(lk, [this] { return remaining == 0; });
      }

    private:
      void WorkerLoop(Int workerId)
      {
        Long lastSeen = 0;
        for (;;)
        {
          void (*body)(Int, Int, void *) = nullptr;
          void *ctx = nullptr;
          Int total = 0;
          Int chunks = 0;
          {
            std::unique_lock<std::mutex> lk(m);
            cv.wait(lk, [this, &lastSeen] { return generation != lastSeen || shutdown; });
            if (shutdown) return;
            lastSeen = generation;
            body = curBody;
            ctx = curCtx;
            total = curTotal;
            chunks = curChunks;
          }

          if (workerId < chunks)
          {
            Int chunkSize = (total + chunks - 1) / chunks;
            Int s = workerId * chunkSize;
            Int e = std::min(total, s + chunkSize);
            if (s < e) body(s, e, ctx);
          }

          {
            std::lock_guard<std::mutex> lk(m);
            if (workerId < chunks)
              --remaining;
          }
          if (workerId < chunks)
            doneCv.notify_one();
        }
      }

      SizeT numThreads = 1;
      std::vector<std::thread> workers;
      std::mutex m;
      std::condition_variable cv;
      std::condition_variable doneCv;
      Long generation = 0;
      bool shutdown = false;

      void (*curBody)(Int, Int, void *) = nullptr;
      void *curCtx = nullptr;
      Int curTotal = 0;
      Int curChunks = 0;
      Int remaining = 0;
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
