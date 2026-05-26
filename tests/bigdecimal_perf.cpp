// BigDecimal performance bench.
//
// Measures Add, Multiply, Divide, Parse, ToString at varying operand scales.
// Reports minimum time over a small iteration count.
//
// Build (standalone): from project root,
//   cmake --build build --target bigdecimal_perf
//   ./build/bigdecimal_perf

#include <chrono>
#include <cstdio>
#include <random>
#include <string>

#include "bigdecimal/BigDecimal.h"

using namespace BigMath;

namespace
{
  std::string RandDigits(int n, std::mt19937 &rng)
  {
    std::string s;
    s.reserve((size_t)n);
    std::uniform_int_distribution<int> d09(0, 9);
    std::uniform_int_distribution<int> d19(1, 9);
    s.push_back((char)('0' + d19(rng)));
    for (int i = 1; i < n; ++i) s.push_back((char)('0' + d09(rng)));
    return s;
  }

  // Build BigDecimal "int.frac" with given digit counts.
  BigDecimal MakeRandom(int intDigits, int fracDigits, std::mt19937 &rng)
  {
    std::string s = RandDigits(intDigits, rng);
    if (fracDigits > 0) { s.push_back('.'); s += RandDigits(fracDigits, rng); }
    return BigDecimal::FromString(s);
  }

  template <class F>
  double TimeMin(int iters, F &&f)
  {
    double best = 1e18;
    for (int i = 0; i < iters; ++i)
    {
      auto t0 = std::chrono::steady_clock::now();
      f();
      auto t1 = std::chrono::steady_clock::now();
      double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
      if (ms < best) best = ms;
    }
    return best;
  }
}

int main()
{
  std::mt19937 rng(42);

  std::printf("BigDecimal performance bench\n");
  std::printf("============================\n\n");

  struct Case { int intDigits, fracDigits, iters; };
  Case sizes[] = {
    { 100,   10,   200 },
    { 1000,  100,  100 },
    { 5000,  500,  20  },
    { 20000, 2000, 5   },
  };

  std::printf("%-22s %12s %12s %12s\n", "operand (int.frac)", "Add ms", "Mul ms", "Div ms (10dp)");

  for (auto const &c : sizes)
  {
    BigDecimal a = MakeRandom(c.intDigits, c.fracDigits, rng);
    BigDecimal b = MakeRandom(c.intDigits, c.fracDigits, rng);

    double tAdd = TimeMin(c.iters, [&]{ volatile auto r = (a + b); (void)r; });
    double tMul = TimeMin(c.iters, [&]{ volatile auto r = (a * b); (void)r; });
    double tDiv = TimeMin(std::max(1, c.iters / 5), [&]{
      volatile auto r = a.Divide(b, 10, RoundingMode::HALF_EVEN); (void)r;
    });

    char label[64];
    std::snprintf(label, sizeof(label), "%d.%d", c.intDigits, c.fracDigits);
    std::printf("%-22s %12.4f %12.4f %12.4f\n", label, tAdd, tMul, tDiv);
  }

  std::printf("\n");
  std::printf("Scale-mismatched Add (int=1000, fracA=10, fracB=1000):\n");
  {
    BigDecimal a = MakeRandom(1000, 10, rng);
    BigDecimal b = MakeRandom(1000, 1000, rng);
    double t = TimeMin(100, [&]{ volatile auto r = (a + b); (void)r; });
    std::printf("  %.4f ms\n", t);
  }

  std::printf("\nParse + ToString:\n");
  for (int n : {100, 1000, 10000}) {
    std::string s = RandDigits(n, rng);
    s.insert(s.size() / 2, ".");
    double tParse = TimeMin(50, [&]{
      volatile auto r = BigDecimal::FromString(s); (void)r;
    });
    BigDecimal v = BigDecimal::FromString(s);
    double tFmt = TimeMin(50, [&]{
      volatile auto r = v.ToPlainString(); (void)r;
    });
    std::printf("  %5d digits: parse %.4f ms, ToString %.4f ms\n", n, tParse, tFmt);
  }

  std::printf("\nDivide at varying target scales (operand=2000 int + 200 frac):\n");
  {
    BigDecimal a = MakeRandom(2000, 200, rng);
    BigDecimal b = MakeRandom(2000, 200, rng);
    for (int scale : {0, 10, 100, 1000, 5000}) {
      double t = TimeMin(3, [&]{
        volatile auto r = a.Divide(b, scale, RoundingMode::HALF_EVEN); (void)r;
      });
      std::printf("  scale=%-5d  %.4f ms\n", scale, t);
    }
  }

  return 0;
}
