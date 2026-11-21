// Single-limb division kernel (ClassicDivision::DivRem1) against plain
// 128-by-64 division, plus the vector entry points built on it.

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/division/ClassicDivision.h"

using namespace BigMath;

#if BIGMATH_LIMB_64

REGISTER_TEST(DivSingleLimb, KernelMatchesWideDivision)
{
  std::mt19937_64 gen(0xD1B1ULL);
  const DataT M = ~(DataT)0;
  // Divisors covering every normalization shift class and the extremes.
  std::vector<DataT> divisors{1, 2, 3, 10, 0xFFFFFFFFULL, (DataT)1 << 32, ((DataT)1 << 63) - 1, (DataT)1 << 63,
                              ((DataT)1 << 63) + 1, M - 1, M, 10000000000000000000ULL, 1000000000000000000ULL};
  for (int i = 0; i < 40; ++i)
    divisors.push_back((gen() >> (gen() % 64)) | 1);

  for (SizeT n : {1, 2, 3, 4, 7, 8, 9, 33, 100})
    for (DataT d : divisors)
      for (int kind = 0; kind < 3; ++kind)
      {
        std::vector<DataT> a(n);
        for (auto &x : a)
          x = kind == 0 ? gen() : (kind == 1 ? M : (gen() & 1 ? M : (DataT)0));
        // Reference.
        std::vector<DataT> refq(n);
        ULong128 rem = 0;
        for (SizeT j = n; j-- > 0;)
        {
          ULong128 cur = (rem << 64) | a[j];
          refq[j] = (DataT)(cur / d);
          rem = cur % d;
        }
        std::vector<DataT> q(n, 0x5555555555555555ULL);
        DataT r = ClassicDivision::DivRem1(q.data(), a.data(), n, d);
        ASSERT_TRUE(q == refq);
        ASSERT_EQ(r, (DataT)rem);
        // In place.
        std::vector<DataT> inplace = a;
        r = ClassicDivision::DivRem1(inplace.data(), inplace.data(), n, d);
        ASSERT_TRUE(inplace == refq);
        ASSERT_EQ(r, (DataT)rem);
      }
}

REGISTER_TEST(DivSingleLimb, VectorEntryPointsAgree)
{
  std::mt19937_64 gen(0x5CA1A2ULL);
  for (SizeT n : {1, 2, 5, 16, 65})
    for (DataT d : {(DataT)7, (DataT)1000000000000000000ULL, ~(DataT)0, (DataT)1 << 63})
    {
      std::vector<DataT> a(n);
      for (auto &x : a)
        x = gen();
      a.back() |= 1;
      auto qr = ClassicDivision::DivideAndRemainder(a, d, Base2_64);
      std::vector<DataT> q = ClassicDivision::Divide(a, d, Base2_64);
      std::vector<DataT> inplace = a;
      DataT r = ClassicDivision::DivModTo(inplace, d, Base2_64);

      // q * d + r == a
      std::vector<DataT> back(n + 1, 0);
      ULong128 carry = r;
      for (SizeT i = 0; i < n; ++i)
      {
        ULong128 t = (ULong128)(i < qr.first.size() ? qr.first[i] : 0) * d + carry;
        back[i] = (DataT)t;
        carry = t >> 64;
      }
      back[n] = (DataT)carry;
      back.pop_back();
      ASSERT_TRUE(back == a);
      ASSERT_EQ(qr.second[0], r);
      auto trim = [](std::vector<DataT> v) { while (v.size() > 1 && v.back() == 0) v.pop_back(); return v; };
      ASSERT_TRUE(trim(q) == trim(qr.first));
      ASSERT_TRUE(trim(inplace) == trim(qr.first));
    }
}

#endif
