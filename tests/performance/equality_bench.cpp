// Compare the former equality implementation with operator== on the same operands.
#include "biginteger/ops/Comparison.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

using namespace BigMath;

namespace
{
  using Predicate = bool (*)(const BigInteger &, const BigInteger &);

  double BatchNs(Predicate predicate, const BigInteger &a, const BigInteger &b, bool expected, int reps)
  {
    // Read the call target on every iteration so the compiler cannot hoist a
    // comparison of two unchanged inputs out of the timed loop.
    volatile Predicate compare = predicate;
    int matches = 0;
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < reps; ++i)
      matches += compare(a, b);
    const auto elapsed = std::chrono::steady_clock::now() - start;
    if (matches != (expected ? reps : 0))
      std::abort();
    return std::chrono::duration<double, std::nano>(elapsed).count();
  }

  void Measure(const char *label, const BigInteger &a, const BigInteger &b, bool expected)
  {
    const Predicate predicates[] = {
      [](const BigInteger &x, const BigInteger &y) { return x.CompareTo(y) == 0; },
      [](const BigInteger &x, const BigInteger &y) { return x == y; }
    };
    std::array<double, 9> samples[2];
    int reps[] = {1, 1};
    for (int p = 0; p < 2; ++p)
      while (BatchNs(predicates[p], a, b, expected, reps[p]) < 1'000'000)
        reps[p] *= 2;
    for (std::size_t i = 0; i < samples[0].size(); ++i)
      for (int step = 0; step < 2; ++step)
      {
        const auto p = (i + step) % 2;
        samples[p][i] = BatchNs(predicates[p], a, b, expected, reps[p]) / reps[p];
      }
    for (auto &times : samples)
      std::sort(times.begin(), times.end());
    std::printf("%u,%s,%.2f,%.2f\n", unsigned(a.size()), label, samples[0][4], samples[1][4]);
  }
}

int main()
{
  std::mt19937_64 random(193);
  std::puts("limbs,comparison,compare_to_ns,equality_ns");
  for (SizeT size : {1u, 4u, 1'000u, 1'000'000u})
  {
    std::vector<DataT> limbs(size);
    for (auto &limb : limbs)
      limb = random() & LimbMask;
    limbs.back() |= 2; // Keep the length when the benchmark flips bit zero.
    const BigInteger a(limbs, false);
    Measure("same object", a, a, true);
    Measure("equal copy", a, BigInteger(a), true);
    Measure("opposite sign", a, BigInteger(limbs, true), false);
    auto mismatch = [&](const char *label, SizeT index) {
      limbs[index] ^= 1;
      Measure(label, a, BigInteger(limbs, false), false);
      limbs[index] ^= 1;
    };
    mismatch("highest limb", size - 1);
    if (size > 1)
      mismatch("lowest limb", 0);
    if (size >= 4)
    {
      mismatch("second highest", size - 2);
      mismatch("second lowest", 1);
    }
    if (size > 4)
      mismatch("middle limb", size / 2);
    limbs.push_back(1);
    Measure("different length", a, BigInteger(limbs, false), false);
  }
}
