// In-place divide-and-conquer division (BurnikelZieglerDivision, 64-bit
// limbs) against the schoolbook basecase and the division identity.
//
// The recursion splits the divisor into floor/ceil halves by pointer
// arithmetic, so the interesting shapes are odd lengths at every level,
// quotient blocks shorter than the divisor (the top-heavy step), dividends
// many blocks long, and operands that drive the "quotient digit one too
// large → add back" corrections and the high-bit (qh) paths.

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/division/BurnikelZieglerDivision.h"
#include "biginteger/algorithms/division/FastDivision.h"
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

  void Check(std::vector<DataT> const &a, std::vector<DataT> const &b)
  {
    auto qr = BurnikelZieglerDivision::DivideAndRemainder(a, b, Base2_64, true);
    auto q = Trimmed(qr.first);
    auto r = Trimmed(qr.second);
    ASSERT_TRUE(Compare(r, b) < 0);
    ASSERT_TRUE(Trimmed(Add(Multiply(q, b, Base2_64), r, Base2_64)) == Trimmed(a));
    auto ref = FastDivision::DivideAndRemainder(a, b, Base2_64, true);
    ASSERT_TRUE(q == Trimmed(ref.first));
    ASSERT_TRUE(r == Trimmed(ref.second));
    ASSERT_TRUE(Trimmed(BurnikelZieglerDivision::DivideAndRemainder(a, b, Base2_64, false).first) == q);
  }

  std::vector<DataT> Random(SizeT n, std::mt19937_64 &gen)
  {
    std::vector<DataT> v(n);
    for (auto &x : v)
      x = gen();
    if (v.back() == 0)
      v.back() = 1;
    return v;
  }
}

#if BIGMATH_LIMB_64

REGISTER_TEST(BzInPlace, ShapesAroundTheRecursion)
{
  std::mt19937_64 gen(0xB21A9EULL);
  const SizeT T = (SizeT)BIGMATH_DC_DIV_THRESHOLD;
  // Divisors from below the basecase floor to several recursion levels deep,
  // including odd lengths and 2^k ± 1.
  for (SizeT dn : {(SizeT)2, (SizeT)3, T - 1, T, T + 1, 2 * T - 1, 2 * T, 2 * T + 1, 3 * T + 2, (SizeT)255, (SizeT)256, (SizeT)257, (SizeT)389, (SizeT)513})
  {
    // Quotient lengths: 1, shorter than the floor, between floor and dn (the
    // top-heavy block), exactly dn, and several blocks with a ragged first one.
    for (SizeT qn : {(SizeT)1, (SizeT)5, T + 3, dn > 4 ? dn - 3 : dn, dn, dn + 1, 2 * dn, 2 * dn + T + 7, 3 * dn + 1})
    {
      auto b = Random(dn, gen);
      auto a = Random(dn + qn, gen);
      if (Compare(a, b) <= 0)
        a.push_back(1);
      Check(a, b);
      b.back() |= (DataT)1 << 63; // already normalized divisor
      if (Compare(a, b) <= 0)
        a.push_back(1);
      Check(a, b);
    }
  }
}

REGISTER_TEST(BzInPlace, StructuredOperands)
{
  const DataT M = ~(DataT)0;
  std::mt19937_64 gen(0x57BCULL);
  for (SizeT dn : {(SizeT)130, (SizeT)200, (SizeT)257, (SizeT)400})
    for (SizeT na : {dn + 1, dn + 100, 2 * dn, 2 * dn + 1, 3 * dn + 17})
    {
      // All ones: every quotient estimate is maximal.
      Check(std::vector<DataT>(na, M), std::vector<DataT>(dn, M));
      // Divisor = high half all ones, low half zero, and the reverse: the
      // low-part correction product is zero / maximal.
      std::vector<DataT> b(dn, 0);
      for (SizeT i = dn / 2; i < dn; ++i)
        b[i] = M;
      Check(std::vector<DataT>(na, M), b);
      b.assign(dn, M);
      for (SizeT i = dn / 2; i < dn; ++i)
        b[i] = 0;
      b.back() = (DataT)1 << 63;
      Check(std::vector<DataT>(na, M), b);
      // Dividend = divisor · B^k − 1: windows whose top limbs equal the
      // divisor's, quotient digits B − 1, and corrections at every level.
      std::vector<DataT> d = Random(dn, gen);
      d.back() |= (DataT)1 << 63;
      std::vector<DataT> a(na, M);
      for (SizeT i = 0; i < dn; ++i)
        a[na - dn + i] = d[i];
      if (Compare(a, d) > 0)
        Check(a, d);
      // Exact multiple of the divisor: zero remainder.
      auto qv = Random(na - dn, gen);
      Check(Trimmed(Multiply(qv, d, Base2_64)), d);
    }
}

#endif
