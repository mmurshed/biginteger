#include "biginteger/ops/Comparison.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>

using namespace BigMath;

namespace
{
  template <class A, class B>
  double Measure(Int (*function)(A const &, B const &), A const &a, B const &b, Int expected)
  {
    auto volatile compare = function;
    const auto batch = [&](int reps) {
      int sum = 0;
      const auto start = std::chrono::steady_clock::now();
      for (int i = 0; i < reps; ++i)
        sum += compare(a, b);
      const auto stop = std::chrono::steady_clock::now();
      if (sum != reps * expected)
        std::abort();
      return std::chrono::duration<double, std::nano>(stop - start).count();
    };
    int reps = 1;
    while (batch(reps) < 500'000 && reps < 1 << 20)
      reps *= 2;
    std::array<double, 7> samples;
    for (double &sample : samples)
      sample = batch(reps) / reps;
    std::sort(samples.begin(), samples.end());
    return samples[samples.size() / 2];
  }

  void Ordered(SizeT size, const char *shape, BigInteger const &a, BigInteger const &b, Int order)
  {
    using Operation = Int (*)(BigInteger const &, BigInteger const &);
    const std::array<std::pair<const char *, Operation>, 5> operations{{
        {"compare_to", [](auto const &x, auto const &y) -> Int { return x.CompareTo(y); }},
        {"greater_equal", [](auto const &x, auto const &y) -> Int { return x >= y; }},
        {"greater", [](auto const &x, auto const &y) -> Int { return x > y; }},
        {"less_equal", [](auto const &x, auto const &y) -> Int { return x <= y; }},
        {"less", [](auto const &x, auto const &y) -> Int { return x < y; }},
    }};
    const Int expected[] = {order, order >= 0, order > 0, order <= 0, order < 0};
    for (size_t i = 0; i < operations.size(); ++i)
    {
      const auto [name, operation] = operations[i];
      const auto ns = Measure(operation, a, b, expected[i]);
      std::printf("%s,%s,%u,%.3f\n", name, shape, size, ns);
    }
  }

  void Magnitude(SizeT size, const char *shape, std::vector<DataT> const &a,
                 std::vector<DataT> const &b, Int order)
  {
    using Vector = std::vector<DataT>;
    const auto ns = Measure<Vector, Vector>(
        [](auto const &x, auto const &y) { return Compare(x, y); }, a, b, order);
    std::printf("compare,%s,%u,%.3f\n", shape, size, ns);
  }
}

int main()
{
  std::mt19937_64 random(197);
  std::puts("operation,shape,limbs,ns");
  for (SizeT size : {1u, 4u, 8u, 16u, 64u, 256u, 4096u, 65536u, 1000000u})
  {
    std::vector<DataT> limbs(size);
    for (auto &limb : limbs)
      limb = random() & LimbMask;
    limbs.back() |= 2;
    const BigInteger a(limbs, false);
    const BigInteger copy(a);
    Ordered(size, "same_object", a, a, 0);
    Ordered(size, "equal_copy", a, copy, 0);
    Magnitude(size, "equal_copy", limbs, copy.GetInteger(), 0);
    Ordered(size, "opposite_sign", a, BigInteger(limbs, true), 1);

    auto mismatch = [&](const char *label, SizeT index) {
      auto other = limbs;
      other[index] ^= 1;
      const Int order = limbs[index] > other[index] ? 1 : -1;
      Ordered(size, label, a, BigInteger(other, false), order);
      Magnitude(size, label, limbs, other, order);
    };
    mismatch("highest_limb", size - 1);
    if (size > 1)
    {
      mismatch("second_highest", size - 2);
      mismatch("lowest_limb", 0);
    }
    if (size > 4)
      mismatch("middle_limb", size / 2);

    auto longer = limbs;
    longer.push_back(1);
    Ordered(size, "different_length", a, BigInteger(longer, false), -1);
    Magnitude(size, "different_length", limbs, longer, -1);
    longer.resize(size + 7, 0);
    longer[size] = 0;
    Magnitude(size, "padded_equal", limbs, longer, 0);

    const auto high = limbs.front();
    std::vector<DataT> padded(size, 0);
    padded.front() = high;
    const auto scalarNs = Measure<std::vector<DataT>, DataT>(
        [](auto const &x, auto const &y) { return Compare(x, y); }, padded, high, 0);
    std::printf("compare_scalar,padded_equal,%u,%.3f\n", size, scalarNs);
  }
}
