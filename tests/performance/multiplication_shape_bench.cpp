// Shape-focused multiplication benchmark for dispatcher tuning.

#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
#include "biginteger/algorithms/multiplication/NTTMultiplication.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/common/Util.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
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
    {
#if BIGMATH_LIMB_64
      value.back() |= 0x8000000000000000ULL;
#else
      value.back() |= 0x80000000ULL;
#endif
    }
    TrimZerosToOne(value);
    return value;
  }

  template <class F>
  double BestMs(F &&fn, int reps)
  {
    double best = std::numeric_limits<double>::max();
    for (int i = 0; i < reps; ++i)
    {
      auto start = std::chrono::steady_clock::now();
      SizeT check = fn();
      auto end = std::chrono::steady_clock::now();
      sink += check;
      best = std::min(best, std::chrono::duration<double, std::milli>(end - start).count());
    }
    return best;
  }

  int RepsFor(SizeT totalLimbs)
  {
    if (totalLimbs <= 256) return 9;
    if (totalLimbs <= 2048) return 7;
    if (totalLimbs <= 8192) return 5;
    if (totalLimbs <= 32768) return 3;
    return 1;
  }

  std::string Winner(
      double classic,
      double karatsuba,
      double ntt,
      double blockwise,
      double dispatch)
  {
    std::string name = "dispatch";
    double best = dispatch;
    if (!std::isnan(classic) && classic < best)
    {
      best = classic;
      name = "classic";
    }
    if (!std::isnan(karatsuba) && karatsuba < best)
    {
      best = karatsuba;
      name = "karatsuba";
    }
    if (!std::isnan(ntt) && ntt < best)
    {
      best = ntt;
      name = "ntt";
    }
    if (!std::isnan(blockwise) && blockwise < best)
    {
      best = blockwise;
      name = "blockwise";
    }
    return name;
  }

  struct Case
  {
    std::string shape;
    SizeT lhsLimbs;
    SizeT rhsLimbs;
  };

  void AddProductAt(
      std::vector<DataT> &out,
      std::vector<DataT> const &product,
      SizeT offset)
  {
    if (product.empty() || IsZero(product))
      return;
    if (out.size() < offset + product.size() + 1)
      out.resize(offset + product.size() + 1, 0);

    SizeT i = 0;
#if BIGMATH_LIMB_64
    ULong128 carry = 0;
    for (; i < product.size(); ++i)
    {
      ULong128 sum = (ULong128)out[offset + i] + product[i] + carry;
      out[offset + i] = (DataT)sum;
      carry = sum >> 64;
    }
    SizeT pos = offset + i;
    while (carry)
    {
      ULong128 sum = (ULong128)out[pos] + carry;
      out[pos] = (DataT)sum;
      carry = sum >> 64;
      ++pos;
      if (pos >= out.size() && carry)
        out.push_back(0);
    }
#else
    ULong carry = 0;
    for (; i < product.size(); ++i)
    {
      ULong sum = (ULong)out[offset + i] + product[i] + carry;
      out[offset + i] = (DataT)(sum & 0xFFFFFFFFULL);
      carry = sum >> 32;
    }
    SizeT pos = offset + i;
    while (carry)
    {
      ULong sum = (ULong)out[pos] + carry;
      out[pos] = (DataT)(sum & 0xFFFFFFFFULL);
      carry = sum >> 32;
      ++pos;
      if (pos >= out.size() && carry)
        out.push_back(0);
    }
#endif
  }

  std::vector<DataT> BlockwiseSkewMultiply(
      std::vector<DataT> const &a,
      std::vector<DataT> const &b,
      BaseT base)
  {
    std::vector<DataT> const *large = &a;
    std::vector<DataT> const *small = &b;
    if (a.size() < b.size())
      std::swap(large, small);

    SizeT blockLimbs = std::max((SizeT)256, (SizeT)small->size() * 4);
    std::vector<DataT> out(large->size() + small->size() + 2, 0);
    for (SizeT start = 0; start < large->size(); start += blockLimbs)
    {
      SizeT end = std::min((SizeT)large->size(), start + blockLimbs);
      std::vector<DataT> chunk(large->begin() + start, large->begin() + end);
      TrimZerosToOne(chunk);
      if (IsZero(chunk))
        continue;
      std::vector<DataT> product = Multiply(chunk, *small, base);
      AddProductAt(out, product, start);
    }
    TrimZerosToOne(out);
    return out;
  }

  void CheckEqual(
      std::vector<DataT> const &actual,
      std::vector<DataT> const &expected,
      std::string const &label)
  {
    if (Compare(actual, expected) != 0)
    {
      std::cerr << "multiplication result mismatch: " << label
                << " actual_limbs=" << actual.size()
                << " expected_limbs=" << expected.size() << '\n';
      std::abort();
    }
  }
}

int main(int argc, char **argv)
{
  std::vector<Case> cases;
  if (argc > 1)
  {
    for (int i = 1; i < argc; ++i)
    {
      SizeT n = (SizeT)std::strtoull(argv[i], nullptr, 10);
      cases.push_back({"n/n", n, n});
      cases.push_back({"2n/n", 2 * n, n});
      cases.push_back({"5n/n", 5 * n, n});
      cases.push_back({"10n/n", 10 * n, n});
      cases.push_back({"20n/n", 20 * n, n});
      cases.push_back({"50n/n", 50 * n, n});
    }
  }
  else
  {
    for (SizeT n : {64u, 96u, 128u, 256u, 512u, 1024u, 2048u, 4096u, 8192u})
      cases.push_back({"n/n", n, n});
    for (SizeT n : {64u, 128u, 256u, 512u, 1024u, 2048u})
    {
      cases.push_back({"5n/n", 5 * n, n});
      cases.push_back({"10n/n", 10 * n, n});
      cases.push_back({"20n/n", 20 * n, n});
    }
  }

  std::cout << "base=" << BigInteger::Base() << "\n";
  std::cout << "shape,lhs_limbs,rhs_limbs,total_limbs,classic_ms,karatsuba_ms,ntt_ms,blockwise_ms,dispatch_ms,winner\n";

  std::mt19937_64 gen(0xA11CEB16B00B5ULL);
  for (auto const &c : cases)
  {
    std::vector<DataT> a = RandomNumber(c.lhsLimbs, gen);
    std::vector<DataT> b = RandomNumber(c.rhsLimbs, gen);
    SizeT total = c.lhsLimbs + c.rhsLimbs;
    int reps = RepsFor(total);

    std::vector<DataT> expected = Multiply(a, b, BigInteger::Base());

    double classic = std::numeric_limits<double>::quiet_NaN();
    if (total <= 4096 || std::min(c.lhsLimbs, c.rhsLimbs) <= 128)
    {
      classic = BestMs([&]() {
        auto r = ClassicMultiplication::Multiply(a, b, BigInteger::Base());
        CheckEqual(r, expected, c.shape + " classic");
        return (SizeT)r.size();
      }, reps);
    }

    double karatsuba = std::numeric_limits<double>::quiet_NaN();
    if (total <= 8192 || std::min(c.lhsLimbs, c.rhsLimbs) <= 128)
    {
      karatsuba = BestMs([&]() {
        auto r = KaratsubaMultiplication::Multiply(a, b, BigInteger::Base());
        CheckEqual(r, expected, c.shape + " karatsuba");
        return (SizeT)r.size();
      }, reps);
    }

    double ntt = BestMs([&]() {
      auto r = NTTMultiplication::Multiply(a, b, BigInteger::Base());
      CheckEqual(r, expected, c.shape + " ntt");
      return (SizeT)r.size();
    }, std::max(1, reps / 2));

    double blockwise = std::numeric_limits<double>::quiet_NaN();
    if (std::max(c.lhsLimbs, c.rhsLimbs) >= 8 * std::min(c.lhsLimbs, c.rhsLimbs))
    {
      blockwise = BestMs([&]() {
        auto r = BlockwiseSkewMultiply(a, b, BigInteger::Base());
        CheckEqual(r, expected, c.shape + " blockwise");
        return (SizeT)r.size();
      }, std::max(1, reps / 2));
    }

    double dispatch = BestMs([&]() {
      auto r = Multiply(a, b, BigInteger::Base());
      CheckEqual(r, expected, c.shape + " dispatch");
      return (SizeT)r.size();
    }, reps);

    std::cout << c.shape << ','
              << c.lhsLimbs << ','
              << c.rhsLimbs << ','
              << total << ','
              << std::fixed << std::setprecision(4)
              << classic << ','
              << karatsuba << ','
              << ntt << ','
              << blockwise << ','
              << dispatch << ','
              << Winner(classic, karatsuba, ntt, blockwise, dispatch)
              << '\n';
  }

  std::cerr << "checksum=" << sink << '\n';
  return 0;
}
