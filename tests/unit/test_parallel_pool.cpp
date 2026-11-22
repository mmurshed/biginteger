// Thread-pool contract tests (common/Parallel.h).
//
// The pool has a single work slot. These tests pin the three properties the
// rest of the library leans on: every index of a dispatch runs exactly once,
// nested dispatches complete (inline), and — the one that used to hang —
// distinct user threads may run pool-backed operations at the same time.

#include "unit_test_framework.h"

#include <atomic>
#include <chrono>
#include <cstdint>
#include <random>
#include <thread>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"
#include "biginteger/algorithms/Squaring.h"
#include "biginteger/common/Parallel.h"

using namespace BigMath;

REGISTER_TEST(ParallelPool, EveryIndexRunsExactlyOnce)
{
  for (Int tasks : {1, 2, 3, 6, 8, 13, 100})
  {
    for (int round = 0; round < 200; ++round)
    {
      std::vector<std::atomic<int>> hits((size_t)tasks);
      for (auto &h : hits)
        h.store(0);
      ParallelDo(tasks, [&hits](Int s, Int e) {
        for (Int i = s; i < e; ++i)
          hits[(size_t)i].fetch_add(1);
      });
      for (Int i = 0; i < tasks; ++i)
        ASSERT_EQ(hits[(size_t)i].load(), 1);
    }
  }
}

REGISTER_TEST(ParallelPool, NestedDispatchCompletes)
{
  std::atomic<long> sum{0};
  ParallelDo(6, [&sum](Int s, Int e) {
    for (Int i = s; i < e; ++i)
    {
      // A dispatch from inside a chunk must not re-enter the pool's slot.
      ParallelDo(5, [&sum, i](Int s2, Int e2) {
        for (Int j = s2; j < e2; ++j)
          sum.fetch_add((long)(i * 10 + j));
      });
    }
  });
  // sum over i in [0,6), j in [0,5) of (10 i + j) = 10*15*5 + 6*10 = 810
  ASSERT_EQ(sum.load(), 810L);
}

REGISTER_TEST(ParallelPool, ConcurrentCallersStayCorrect)
{
  // Four user threads multiply NTT-sized operands at the same time. Before
  // callers were serialized on the pool's work slot, a second dispatcher
  // overwrote the job in flight: lost chunks, a corrupted completion count,
  // and in practice a hang.
  const int threads = 4;
  const int rounds = 40;
  BaseT base = BigInteger::Base();
  std::vector<std::vector<DataT>> a((size_t)threads), b((size_t)threads), expect((size_t)threads);
  for (int t = 0; t < threads; ++t)
  {
    std::mt19937_64 gen(0xB16B00B5ULL + (uint64_t)t);
    a[(size_t)t].resize((size_t)(3000 + 37 * t));
    b[(size_t)t].resize((size_t)(2800 + 11 * t));
    for (auto &x : a[(size_t)t]) x = gen() & LimbMask;
    for (auto &x : b[(size_t)t]) x = gen() & LimbMask;
    expect[(size_t)t] = KaratsubaMultiplication::Multiply(a[(size_t)t], b[(size_t)t], base);
  }

  std::atomic<int> wrong{0};
  std::vector<std::thread> pool;
  for (int t = 0; t < threads; ++t)
  {
    pool.emplace_back([&, t] {
      for (int r = 0; r < rounds; ++r)
      {
        auto c = Multiply(a[(size_t)t], b[(size_t)t], base);
        if (c != expect[(size_t)t])
          wrong.fetch_add(1);
      }
    });
  }
  for (auto &th : pool)
    th.join();
  ASSERT_EQ(wrong.load(), 0);
}

#if BIGMATH_USE_THREADS

REGISTER_TEST(ParallelPool, ResultsIndependentOfTemperature)
{
  // The dispatcher picks Karatsuba or the NTT for mid-sized products
  // depending on whether the pool is awake. Both must give the same product,
  // whether calls arrive in a burst (pool warms up) or far apart (cold).
  BaseT base = BigInteger::Base();
  std::mt19937_64 gen(0x7E3BULL);
  for (SizeT n : {320, 400, 500, 600})
  {
    std::vector<DataT> a(n), b(n);
    for (auto &x : a) x = gen() & LimbMask;
    for (auto &x : b) x = gen() & LimbMask;
    auto expect = KaratsubaMultiplication::Multiply(a, b, base);
    auto expectSq = KaratsubaMultiplication::Multiply(a, a, base);
    for (int i = 0; i < 12; ++i) // burst
    {
      ASSERT_TRUE(Multiply(a, b, base) == expect);
      ASSERT_TRUE(Square(a, base) == expectSq);
    }
    for (int i = 0; i < 3; ++i) // spaced out
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(2));
      ASSERT_TRUE(Multiply(a, b, base) == expect);
      std::this_thread::sleep_for(std::chrono::milliseconds(2));
      ASSERT_TRUE(Square(a, base) == expectSq);
    }
  }
}

#endif


// Awake workers pick a job up without a lock, from a description the slot
// owner rewrites for every dispatch. Back-to-back jobs with different
// contexts and slice counts must never run a slice with another job's
// description: each job's context lives in a stack frame that is gone by the
// next dispatch, and carries the id the body expects.
REGISTER_TEST(ParallelPool, BackToBackJobsNeverMixDescriptions)
{
  std::mt19937 gen(123);
  for (int round = 0; round < 20000; ++round)
  {
    const Int tasks = 2 + (Int)(gen() % 9);
    struct Ctx
    {
      int id;
      Int tasks;
      std::atomic<int> hits[16];
      std::atomic<int> bad{0};
    } ctx;
    ctx.id = round;
    ctx.tasks = tasks;
    for (auto &h : ctx.hits)
      h.store(0);
    ParallelDo(tasks, [&ctx, round, tasks](Int s, Int e) {
      for (Int i = s; i < e; ++i)
      {
        if (i < 0 || i >= tasks || ctx.id != round || ctx.tasks != tasks)
          ctx.bad.fetch_add(1);
        else
          ctx.hits[(size_t)i].fetch_add(1);
      }
    });
    ASSERT_EQ(ctx.bad.load(), 0);
    for (Int i = 0; i < tasks; ++i)
      ASSERT_EQ(ctx.hits[(size_t)i].load(), 1);
  }
}

// No lost wake-ups. Each job has one slice per pool thread, and a slice does
// not return until every slice has been entered — so the job only completes
// if every worker shows up, wherever the dispatch caught it: spinning, on
// its way to park, parked, or just woken. The gaps between jobs straddle the
// spin window to hit each of those states.
REGISTER_TEST(ParallelPool, EveryWorkerAnswersEveryDispatch)
{
  const Int threads = (Int)ParallelNumThreads();
  if (threads < 2)
    return;
  using Clock = std::chrono::steady_clock;
  for (int gapUs : {0, 5, 20, 50, 80, 95, 100, 105, 120, 200, 500, 3000})
    for (int round = 0; round < 25; ++round)
    {
      const auto until = Clock::now() + std::chrono::microseconds(gapUs);
      while (Clock::now() < until)
      {
      }
      std::atomic<Int> entered{0};
      std::atomic<bool> timedOut{false};
      ParallelDo(threads, [&](Int s, Int e) {
        entered.fetch_add(e - s);
        const auto deadline = Clock::now() + std::chrono::seconds(20);
        while (entered.load() < threads)
        {
          if (Clock::now() > deadline)
          {
            timedOut.store(true);
            break;
          }
          std::this_thread::yield();
        }
      });
      ASSERT_TRUE(!timedOut.load());
      ASSERT_EQ(entered.load(), threads);
    }
}

#if BIGMATH_USE_THREADS
REGISTER_TEST(ParallelPool, ShareReadOnlyTransformPlans)
{
  const auto &held = NttCrt::GetPlan<NttCrt::F2, NttCrt::G2>(32768);
  constexpr int count = 4;
  std::array<const NttCrt::Plan<NttCrt::F2> *, count> plans{};
  std::vector<std::thread> threads;
  for (int i = 0; i < count; ++i)
  {
    threads.emplace_back([&, i] {
      plans[i] = &NttCrt::GetPlan<NttCrt::F2, NttCrt::G2>(32768);
    });
  }
  for (auto &thread : threads) thread.join();
  for (auto *plan : plans) ASSERT_TRUE(plan == &held);
}
#endif
