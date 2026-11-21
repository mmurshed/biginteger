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
    for (auto &x : a[(size_t)t]) x = gen();
    for (auto &x : b[(size_t)t]) x = gen();
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

namespace
{
  // Poll `cond` for up to two seconds (CI runners can be slow to schedule).
  template <typename Cond>
  bool Eventually(Cond &&cond)
  {
    auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(2);
    while (std::chrono::steady_clock::now() < deadline)
    {
      if (cond())
        return true;
      std::this_thread::sleep_for(std::chrono::microseconds(200));
    }
    return cond();
  }
}

REGISTER_TEST(ParallelPool, TemperatureTracksParkAndPrewarm)
{
  if (ParallelNumThreads() < 2)
    return; // no workers: the pool is never warm
  // Idle long enough and every worker parks.
  ASSERT_TRUE(Eventually([] { return !ParallelWarm(); }));
  // A prewarm wakes them without any work being dispatched...
  bool warmed = false;
  for (int attempt = 0; attempt < 200 && !warmed; ++attempt)
  {
    ParallelPrewarm();
    auto until = std::chrono::steady_clock::now() + std::chrono::microseconds(80);
    while (std::chrono::steady_clock::now() < until && !warmed)
      warmed = ParallelWarm();
  }
#if !defined(BIGMATH_POOL_SPIN_US) || BIGMATH_POOL_SPIN_US > 0
  ASSERT_TRUE(warmed);
#endif
  // ...and they park again once nothing arrives.
  ASSERT_TRUE(Eventually([] { return !ParallelWarm(); }));
}

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
    for (auto &x : a) x = gen();
    for (auto &x : b) x = gen();
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

