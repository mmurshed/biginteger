// Repeated same-operand benchmark for the prepared CRT NTT API.

#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/multiplication/NTTMultiplication.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/common/Util.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace BigMath;

namespace
{
  volatile SizeT sink = 0;

  std::vector<DataT> RandomNumber(SizeT limbs, std::mt19937_64 &gen)
  {
    std::vector<DataT> value(limbs);
    for (SizeT i = 0; i < limbs; ++i)
    {
#if BIGMATH_LIMB_64
      value[i] = (DataT)gen();
#else
      value[i] = (DataT)(gen() & 0xFFFFFFFFULL);
#endif
    }
    if (!value.empty())
      value.back() |= (DataT)1 << (LimbBits - 1);
    TrimZerosToOne(value);
    return value;
  }

  template <class F>
  double TimeMs(F &&fn)
  {
    auto start = std::chrono::steady_clock::now();
    SizeT check = fn();
    auto end = std::chrono::steady_clock::now();
    sink += check;
    return std::chrono::duration<double, std::milli>(end - start).count();
  }

  void CheckEqual(
      const std::vector<DataT> &actual,
      const std::vector<DataT> &expected,
      SizeT index)
  {
    if (Compare(actual, expected) != 0)
    {
      std::cerr << "prepared NTT mismatch at index " << index
                << " actual_limbs=" << actual.size()
                << " expected_limbs=" << expected.size() << '\n';
      std::abort();
    }
  }
}

int main(int argc, char **argv)
{
  SizeT preparedLimbs = argc > 1 ? (SizeT)std::strtoull(argv[1], nullptr, 10) : 100000;
  SizeT otherLimbs = argc > 2 ? (SizeT)std::strtoull(argv[2], nullptr, 10) : preparedLimbs;
  SizeT count = argc > 3 ? (SizeT)std::strtoull(argv[3], nullptr, 10) : 5;

  std::mt19937_64 gen(0xC0FFEE123456789ULL);
  std::vector<DataT> fixed = RandomNumber(preparedLimbs, gen);
  std::vector<std::vector<DataT>> others;
  others.reserve(count);
  for (SizeT i = 0; i < count; ++i)
    others.push_back(RandomNumber(otherLimbs, gen));

  std::vector<std::vector<DataT>> expected;
  expected.reserve(count);

  double normalMs = TimeMs([&]() {
    SizeT total = 0;
    for (SizeT i = 0; i < count; ++i)
    {
      auto product = NTTMultiplication::Multiply(fixed, others[i], BigInteger::Base());
      total += (SizeT)product.size();
      expected.push_back(std::move(product));
    }
    return total;
  });

  NTTMultiplication::PreparedOperand prepared;
  double prepareMs = TimeMs([&]() {
    prepared = NTTMultiplication::PrepareOperand(fixed, otherLimbs, BigInteger::Base());
    return (SizeT)prepared.n;
  });

  double preparedMs = TimeMs([&]() {
    SizeT total = 0;
    for (SizeT i = 0; i < count; ++i)
    {
      auto product = NTTMultiplication::Multiply(prepared, others[i]);
      CheckEqual(product, expected[i], i);
      total += (SizeT)product.size();
    }
    return total;
  });

  std::cout << "base=" << BigInteger::Base()
            << " prepared_limbs=" << preparedLimbs
            << " other_limbs=" << otherLimbs
            << " count=" << count
            << " transform_n=" << prepared.n << '\n';
  std::cout << std::fixed << std::setprecision(3)
            << "normal_ms=" << normalMs
            << " prepare_ms=" << prepareMs
            << " prepared_ms=" << preparedMs
            << " prepared_with_setup_ms=" << (prepareMs + preparedMs)
            << " steady_state_speedup=" << (normalMs / preparedMs)
            << " amortized_speedup=" << (normalMs / (prepareMs + preparedMs))
            << '\n';
  std::cerr << "checksum=" << sink << '\n';
  return 0;
}
