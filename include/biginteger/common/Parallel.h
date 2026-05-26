/**
 * BigMath: Internal parallel-for helper for NTT-bound ops.
 *
 * Gated on BIGMATH_USE_THREADS=1. When unset (default), all calls reduce to
 * the serial body inline — zero overhead.
 *
 * Implementation uses a small persistent thread pool defined in
 * src/common/Parallel.cpp. The pool is created on first use (lazy init),
 * sized to min(hardware_concurrency(), BIGMATH_MAX_THREADS). Each worker
 * spins on a condition variable; dispatch latency is ~5µs per call.
 *
 * Public header stays free of <thread>/<mutex> includes so consumers don't
 * pick up pthread linkage unless they opt in via BIGMATH_USE_THREADS=1.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGMATH_PARALLEL
#define BIGMATH_PARALLEL

#include "Util.h"

namespace BigMath
{
#if BIGMATH_USE_THREADS

  // Number of worker threads in the pool, fixed for lifetime of process.
  // Workers + the calling thread cooperate, so effective parallelism is
  // NumThreads() in total (workers do chunks 1..N-1, caller does chunk 0).
  SizeT ParallelNumThreads();

  // Run `body(start, end)` over chunks of [0, total). The interval is split
  // into approximately NumThreads()-sized chunks; small `total` may use fewer
  // chunks (or stay serial) to avoid dispatch overhead exceeding work cost.
  //
  // `body` must be thread-safe with respect to other concurrent invocations
  // of itself — each chunk gets disjoint [start, end). Capture-by-reference
  // any shared state.
  //
  // Threshold: if `total < ParallelMinSize()`, runs serially in the caller.
  void ParallelFor(Int total, void (*body)(Int start, Int end, void *ctx), void *ctx);

  // Convenience overload for stateful lambdas. Compiles to the same
  // function-pointer dispatch via a trampoline.
  template <typename F>
  inline void ParallelFor(Int total, F &&body)
  {
    struct Ctx { F f; };
    Ctx ctx{std::forward<F>(body)};
    auto tramp = +[](Int s, Int e, void *p) {
      static_cast<Ctx *>(p)->f(s, e);
    };
    ParallelFor(total, tramp, &ctx);
  }

  // Cutoff below which ParallelFor runs serially. Tuned to roughly match
  // the per-dispatch overhead (~5µs * N threads).
  SizeT ParallelMinSize();

  // ParallelDo: dispatch numTasks chunks unconditionally (bypass MinSize).
  // For coarse-grained parallelism where each task is large enough that
  // dispatch overhead is irrelevant — e.g. running 3 independent NTTs in
  // parallel, one per prime.
  void ParallelDo(Int numTasks, void (*body)(Int start, Int end, void *ctx), void *ctx);

  template <typename F>
  inline void ParallelDo(Int numTasks, F &&body)
  {
    struct Ctx { F f; };
    Ctx ctx{std::forward<F>(body)};
    auto tramp = +[](Int s, Int e, void *p) {
      static_cast<Ctx *>(p)->f(s, e);
    };
    ParallelDo(numTasks, tramp, &ctx);
  }

#else

  // Stubs for the single-threaded build. Always returns 1 / runs serially.
  inline SizeT ParallelNumThreads() { return 1; }
  inline SizeT ParallelMinSize() { return 0; }

  template <typename F>
  inline void ParallelFor(Int total, F &&body)
  {
    if (total > 0)
      body(Int{0}, total);
  }

  template <typename F>
  inline void ParallelDo(Int numTasks, F &&body)
  {
    if (numTasks > 0)
      body(Int{0}, numTasks);
  }

#endif
}

#endif
