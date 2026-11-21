// Newton division: the reciprocal's accuracy contract, and "no silent
// schoolbook fallback".
//
// ApproxReciprocal builds its result with Newton steps along a guarded
// precision schedule. The contract is 0 <= floor(B^(2n)/D) - R <= 3 at every
// size; an earlier exact-doubling schedule drifted by ~B^6 at n = 1024 and was
// patched with an extra full-size iteration. NewtonDivision::FallbackCount()
// counts the times a division gave up on the reciprocal path; the result is
// still right when that happens, only quadratic, so it has to be asserted.

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/Subtraction.h"
#include "biginteger/algorithms/division/BurnikelZieglerDivision.h"
#include "biginteger/algorithms/division/FastDivision.h"
#include "biginteger/algorithms/division/NewtonDivision.h"
#include "biginteger/common/Comparator.h"

using namespace BigMath;

namespace
{
  std::vector<DataT> Trimmed(std::vector<DataT> v)
  {
    while (v.size() > 1 && v.back() == 0)
      v.pop_back();
    return v;
  }

  // kind 0: random; 1: all ones; 2: B^n/2 + 1 (reciprocal just under 2·B^n);
  // 3: B^n/2 exactly; 4: random with a zeroed low half.
  std::vector<DataT> Divisor(SizeT n, int kind, std::mt19937_64 &gen)
  {
    std::vector<DataT> d(n);
    for (auto &x : d)
      x = (kind == 0 || kind == 4) ? gen() : (kind == 1 ? ~(DataT)0 : 0);
    if (kind == 4)
      for (SizeT i = 0; i < n / 2; ++i)
        d[i] = 0;
    d.back() |= (DataT)1 << 63;
    if (kind == 2)
      d[0] |= 1;
    return d;
  }
}

#if BIGMATH_LIMB_64

REGISTER_TEST(NewtonDivision, ReciprocalWithinThreeUlpsAndNeverOver)
{
  std::mt19937_64 gen(0x2EC1920CULL);
  const std::vector<DataT> three{3};
  for (SizeT n : {2, 3, 5, 16, 17, 31, 32, 33, 34, 63, 64, 65, 66, 100, 129, 130, 255, 256, 257, 258, 400, 511, 512, 513, 514,
                  639, 640, 641, 800, 1023, 1024, 1025, 1026, 1280, 1281, 1500, 2047, 2048, 2049, 3000, 4096, 5000})
  {
    for (int kind = 0; kind < 5; ++kind)
    {
      auto D = Divisor(n, kind, gen);
      // Exact floor(B^(2n) / D) by divide-and-conquer division.
      std::vector<DataT> num(2 * n + 1, 0);
      num[2 * n] = 1;
      auto exact = Trimmed(BurnikelZieglerDivision::DivideAndRemainder(num, D, Base2_64, false).first);
      auto R = Trimmed(NewtonDivision::ApproxReciprocal(D, false));
      ASSERT_TRUE(Compare(R, exact) <= 0);                          // never an overestimate
      ASSERT_TRUE(Compare(Add(R, three, Base2_64), exact) >= 0);    // within 3 ulps
    }
  }
}

REGISTER_TEST(NewtonDivision, MatchesReferenceWithoutFallback)
{
  std::mt19937_64 gen(0x9E77DEEDULL);
  const long before = NewtonDivision::FallbackCount();
  for (SizeT n : {2, 7, 33, 64, 65, 130, 257, 400, 513, 640, 641, 1000, 1024, 1025, 1281, 2048, 2049, 3001})
    for (double ratio : {1.02, 1.3, 1.99, 2.0, 2.01, 3.0, 4.5, 7.0})
      for (int kind = 0; kind < 5; ++kind)
      {
        SizeT na = (SizeT)(n * ratio + 0.5);
        if (na <= n)
          na = n + 1;
        auto b = Divisor(n, kind, gen);
        if (kind == 0)
          b.back() = gen() | 1; // un-normalized top limb as well
        std::vector<DataT> a(na);
        for (auto &x : a)
          x = kind == 1 ? ~(DataT)0 : gen();
        if (kind == 3)
        {
          // Dividend whose top limbs equal the divisor's: quotient digits B − 1.
          for (SizeT i = 0; i < n; ++i)
            a[na - n + i] = b[i];
          for (SizeT i = 0; i < na - n; ++i)
            a[i] = ~(DataT)0;
        }
        if (a.back() == 0)
          a.back() = 1;

        auto qr = NewtonDivision::DivideAndRemainder(a, b, Base2_64, true);
        auto q = Trimmed(qr.first), r = Trimmed(qr.second);
        ASSERT_TRUE(Compare(r, Trimmed(b)) < 0);
        ASSERT_TRUE(Trimmed(Add(Multiply(q, b, Base2_64), r, Base2_64)) == Trimmed(a));
        auto ref = BurnikelZieglerDivision::DivideAndRemainder(a, b, Base2_64, true);
        ASSERT_TRUE(q == Trimmed(ref.first));
        ASSERT_TRUE(r == Trimmed(ref.second));
      }
  ASSERT_EQ(NewtonDivision::FallbackCount(), before);
}

REGISTER_TEST(NewtonDivision, CachedDividerMatchesAndDoesNotFallBack)
{
  std::mt19937_64 gen(0xD171DE2ULL);
  const long before = NewtonDivision::FallbackCount();
  for (SizeT n : {64, 257, 1024, 2049})
  {
    auto b = Divisor(n, 0, gen);
    b.back() = gen() | 1;
    NewtonDivision::Divider divider(b, Base2_64);
    for (SizeT na : {n + 1, n + n / 2, 2 * n, 2 * n + 1, 3 * n + 5, 6 * n})
    {
      std::vector<DataT> a(na);
      for (auto &x : a)
        x = gen();
      if (a.back() == 0)
        a.back() = 1;
      auto qr = divider.DivideAndRemainder(a, true);
      auto ref = BurnikelZieglerDivision::DivideAndRemainder(a, b, Base2_64, true);
      ASSERT_TRUE(Trimmed(qr.first) == Trimmed(ref.first));
      ASSERT_TRUE(Trimmed(qr.second) == Trimmed(ref.second));
    }
  }
  ASSERT_EQ(NewtonDivision::FallbackCount(), before);
}

#endif
