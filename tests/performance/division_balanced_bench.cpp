// Near-balanced large-divisor division benchmark.
//
// Times the public dispatch path DivideAndRemainder against the individual
// NewtonDivision / BurnikelZieglerDivision strategies for near-balanced shapes
// (ratio 1.33–2) at large divisor sizes — the band where the dispatcher's
// strategy choice dominates wall-clock. FastDivision is intentionally excluded
// (quadratic; multi-second at these sizes).
//
// Build:
//   c++ -std=c++20 -O3 -march=native -Iinclude \
//       tests/performance/division_balanced_bench.cpp src/**/*.cpp -o division_balanced_bench
//
// Output is CSV: na_limbs,nb_limbs,ratio,dispatch_ms,newton_ms,bz_ms

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <random>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/Division.h"
#include "biginteger/algorithms/division/BurnikelZieglerDivision.h"
#include "biginteger/algorithms/division/NewtonDivision.h"

using namespace std;
using namespace BigMath;
using Clock = chrono::steady_clock;

static vector<DataT> RandomNumber(size_t n, uint64_t seed)
{
  mt19937_64 gen(seed);
  vector<DataT> v(n);
  for (auto &x : v)
    v[&x - v.data()] = (DataT)gen();
#if BIGMATH_LIMB_64
  v.back() |= 0x8000000000000000ULL;
#else
  for (auto &x : v) x &= 0xFFFFFFFFULL;
  v.back() |= 0x80000000ULL;
#endif
  return v;
}

template <class F>
static double BestMs(F &&fn, int reps)
{
  double best = numeric_limits<double>::max();
  for (int i = 0; i < reps; ++i)
  {
    auto t0 = Clock::now();
    auto r = fn();
    auto t1 = Clock::now();
    if (r.first.empty()) abort();
    best = min(best, chrono::duration<double, milli>(t1 - t0).count());
  }
  return best;
}

int main(int argc, char **argv)
{
  BaseT base = BigInteger::Base();
  printf("base=%d\n", (int)base);
  printf("na_limbs,nb_limbs,ratio,dispatch_ms,newton_ms,bz_ms\n");

  // ratio-2 (2n/n) shapes at realistic, non-power-of-2 divisor sizes — the
  // common case. BZ's recursive halving lands intermediate NTT multiplies just
  // over power-of-2 boundaries here, where its constant factor blows up 5–60×;
  // Newton pads once and stays flat. The two power-of-2 sizes are BZ's best
  // case (it ties Newton) and the 2^k+1 size is its worst (≈60×).
  vector<pair<size_t, size_t>> cases = {
      {200000, 100000},   // ratio 2, non-pow2
      {500000, 250000},   // ratio 2, non-pow2
      {1000000, 500000},  // ratio 2, non-pow2
      {2000000, 1000000}, // ratio 2, non-pow2
      {524288, 262144},   // ratio 2, pow2 (BZ best case — control)
      {524290, 262145},   // ratio 2, 2^18+1 (BZ worst case)
      {3000000, 1000000}, // ratio 3, non-pow2 (control: Newton pre-change too)
  };

  for (auto [na, nb] : cases)
  {
    auto a = RandomNumber(na, 0x1111 ^ na);
    auto b = RandomNumber(nb, 0x2222 ^ nb);
    int reps = na <= 300000 ? 3 : 1;

    double dispatch = BestMs([&]() { return DivideAndRemainder(a, b, base, true); }, reps);
    double newton = BestMs([&]() { return NewtonDivision::DivideAndRemainder(a, b, base, true); }, reps);
    double bz = BestMs([&]() { return BurnikelZieglerDivision::DivideAndRemainder(a, b, base, true); }, reps);

    printf("%zu,%zu,%.2f,%.3f,%.3f,%.3f\n", na, nb, (double)na / nb, dispatch, newton, bz);
    fflush(stdout);
  }
  return 0;
}
