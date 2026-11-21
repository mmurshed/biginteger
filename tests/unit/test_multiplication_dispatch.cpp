// Top-level Multiply dispatch across the shape grid.
//
// The dispatcher chooses between schoolbook, (sliced) Karatsuba and the CRT
// NTT from the operand lengths: total, shorter operand and la*lb. Whatever it
// chooses must be the right product, in particular along the rule's edges and
// for both operand orders. Reference: schoolbook.

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/Multiplication.h"

using namespace BigMath;

namespace
{
  std::vector<DataT> Limbs(SizeT n, int kind, std::mt19937_64 &gen)
  {
    std::vector<DataT> v(n);
    for (auto &x : v)
      x = kind == 0 ? gen() : ~(DataT)0;
    if (v.back() == 0)
      v.back() = 1;
    return v;
  }

  void Check(SizeT la, SizeT lb, std::mt19937_64 &gen)
  {
    for (int kind = 0; kind < 2; ++kind)
    {
      auto a = Limbs(la, kind, gen);
      auto b = Limbs(lb, kind, gen);
      auto expect = ClassicMultiplication::Multiply(a, b, Base2_64);
      ASSERT_TRUE(Multiply(a, b, Base2_64) == expect);
      ASSERT_TRUE(Multiply(b, a, Base2_64) == expect);
    }
  }
}

#if BIGMATH_LIMB_64

REGISTER_TEST(MultiplyDispatch, ShapeGridMatchesSchoolbook)
{
  std::mt19937_64 gen(0xD15BA7C4ULL);
  for (SizeT lb : {2, 8, 23, 24, 25, 28, 29, 32, 48, 63, 64, 65, 96, 128, 200, 310, 320})
    for (SizeT ratio : {1, 2, 3, 5, 10, 17})
      Check(lb * ratio, lb, gen);
}

REGISTER_TEST(MultiplyDispatch, NttRuleEdges)
{
  std::mt19937_64 gen(0xED6E0FA7EAULL);
  const SizeT minLimb = NTT_MIN_LIMB_THRESHOLD;
  const ULong area = NTT_AREA_THRESHOLD;
  const SizeT total = NTT_MULTIPLICATION_THRESHOLD;
  if (minLimb == 0 || area == 0)
    return;
  // Either side of the shorter-operand floor, at a length that satisfies the
  // area condition.
  for (SizeT lb : {minLimb - 1, minLimb, minLimb + 1})
  {
    SizeT la = (SizeT)(area / lb) + 1;
    Check(la - 2, lb, gen);
    Check(la, lb, gen);
    Check(la + 2, lb, gen);
  }
  // Either side of the area boundary for a mid-sized shorter operand, and the
  // balanced corner.
  for (SizeT lb : {(SizeT)100, (SizeT)160, (SizeT)256})
  {
    SizeT la = (SizeT)(area / lb);
    Check(la - 1, lb, gen);
    Check(la + 1, lb, gen);
  }
  Check(total / 2 - 1, total / 2 - 1, gen);
  Check(total / 2, total / 2, gen);
  Check(total / 2 + 1, total / 2 + 1, gen);
  // Long and thin: stays off the NTT however long it gets.
  Check(9000, minLimb - 1, gen);
  Check(9000, 24, gen);
  Check(9000, 25, gen);
}

#endif
