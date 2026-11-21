// Top-level DivideAndRemainder dispatch across the shape grid.
//
// The dispatcher chooses between schoolbook, in-place Burnikel–Ziegler,
// Newton and quotient-sized division from the divisor length and the
// quotient length. Whatever it chooses must be the right quotient and
// remainder — in particular on both sides of every edge of the rule — and
// Newton must never fall back to schoolbook on the way.

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Division.h"
#include "biginteger/algorithms/Multiplication.h"
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

  void Check(SizeT na, SizeT nb, int kind, std::mt19937_64 &gen)
  {
    std::vector<DataT> a(na), b(nb);
    for (auto &x : a) x = kind == 1 ? ~(DataT)0 : gen();
    for (auto &x : b) x = kind == 1 ? ~(DataT)0 : gen();
    if (kind == 2)
    {
      // Dividend's top limbs equal the divisor's: quotient digits B − 1.
      for (SizeT i = 0; i < nb; ++i) a[na - nb + i] = b[i];
      for (SizeT i = 0; i < na - nb; ++i) a[i] = ~(DataT)0;
    }
    if (a.back() == 0) a.back() = 1;
    if (b.back() == 0) b.back() = 1;

    auto qr = DivideAndRemainder(a, b, Base2_64, true);
    auto q = Trimmed(qr.first), r = Trimmed(qr.second);
    ASSERT_TRUE(Compare(r, Trimmed(b)) < 0);
    ASSERT_TRUE(Trimmed(Add(Multiply(q, b, Base2_64), r, Base2_64)) == Trimmed(a));
    ASSERT_TRUE(Trimmed(Divide(a, b, Base2_64)) == q);
  }
}

#if BIGMATH_LIMB_64

REGISTER_TEST(DivisionDispatch, ShapeGrid)
{
  std::mt19937_64 gen(0xD1F5A7C4ULL);
  const long before = NewtonDivision::FallbackCount();
  for (SizeT nb : {2, 7, 31, 95, 96, 97, 200, 500, 767, 768, 769, 1100})
    for (double ratio : {1.0, 1.05, 1.5, 1.99, 2.0, 2.01, 2.99, 3.0, 3.01, 6.0})
    {
      SizeT na = ratio == 1.0 ? nb + 1 : (SizeT)(nb * ratio + 0.5);
      if (na <= nb) na = nb + 1;
      for (int kind = 0; kind < 3; ++kind)
        Check(na, nb, kind, gen);
    }
  ASSERT_EQ(NewtonDivision::FallbackCount(), before);
}

REGISTER_TEST(DivisionDispatch, RuleEdges)
{
  std::mt19937_64 gen(0xED6ED1FULL);
  const long before = NewtonDivision::FallbackCount();
  const SizeT newtonB = DIV_NEWTON_MIN_B;      // ratio >= 2
  const SizeT near2B = DIV_NEWTON_NEAR2_MIN_B; // ratio >= 7/4
  const SizeT qMin = DIV_QSIZED_MIN_QUOTIENT;  // quotient limbs

  // Around the ratio-2 Newton floor, on both sides of ratio 2.
  for (SizeT nb : {newtonB - 1, newtonB, newtonB + 1})
    for (SizeT na : {2 * nb - 1, 2 * nb, 2 * nb + 1, 3 * nb})
      Check(na, nb, 0, gen);

  // Around the near-2 floor, on both sides of ratio 7/4.
  for (SizeT nb : {near2B - 1, near2B, near2B + 2})
    for (SizeT na : {nb + (3 * nb) / 4 - 1, nb + (3 * nb) / 4, nb + (3 * nb) / 4 + 1})
      Check(na, nb, 0, gen);

  // Around the quotient-sized entry: quotient just below / at / above the
  // threshold against a divisor comfortably longer than it, and against one
  // only just longer (the delta + guard < nb contract).
  for (SizeT delta : {qMin - 1, qMin, qMin + 1})
  {
    Check(2 * qMin + delta, 2 * qMin, 0, gen);
    Check(2 * qMin + delta, 2 * qMin, 2, gen);
    for (SizeT margin : {1, 7, 8, 9, 40})
      Check(delta + margin + delta, delta + margin, 0, gen);
  }
  ASSERT_EQ(NewtonDivision::FallbackCount(), before);
}

#endif
