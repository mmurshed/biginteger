// Karatsuba on unbalanced operands.
//
// When one operand is at least twice the other, the Karatsuba split used to
// clamp to the short operand and peel it off the long one level by level —
// recursion depth la/lb instead of log2 — and for la > ~16*lb the parked
// temporaries ran past the 16n workspace (heap overflow; e.g. 1214 x 65
// limbs, a shape the dispatcher routes to Karatsuba). Such operands are now
// sliced into balanced pieces. These shapes must simply produce the right
// product; the overflow itself shows up under AddressSanitizer.

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"

using namespace BigMath;

namespace
{
  std::vector<DataT> Limbs(SizeT n, int kind, std::mt19937_64 &gen)
  {
    std::vector<DataT> v(n);
    for (auto &x : v)
      x = kind == 0 ? gen() : (kind == 1 ? ~(DataT)0 : ((gen() & 1) ? ~(DataT)0 : (DataT)1));
    if (v.back() == 0)
      v.back() = 1;
    return v;
  }

  void Check(SizeT la, SizeT lb, int kind, std::mt19937_64 &gen)
  {
    auto a = Limbs(la, kind, gen);
    auto b = Limbs(lb, kind, gen);
    auto expect = ClassicMultiplication::Multiply(a, b, Base2_64);
    ASSERT_TRUE(KaratsubaMultiplication::Multiply(a, b, Base2_64) == expect);
    ASSERT_TRUE(KaratsubaMultiplication::Multiply(b, a, Base2_64) == expect);
  }
}

REGISTER_TEST(KaratsubaUnbalanced, SkewedShapesMatchClassic)
{
  std::mt19937_64 gen(0x5CE3ED5EULL);
  // Short operand just above the leaf cutoff and well above it; long operand
  // from "barely unbalanced" to the far end of the Karatsuba dispatch band.
  for (SizeT lb : {25, 26, 33, 34, 40, 49, 65, 66, 70, 100, 129})
  {
    for (SizeT la : {2 * lb - 1, 2 * lb, 2 * lb + 1, 3 * lb, 3 * lb + 7, 5 * lb + 3, 10 * lb - 1, 17 * lb + 5})
    {
      if (la > 1300)
        la = 1300 - lb % 7;
      for (int kind = 0; kind < 3; ++kind)
        Check(la, lb, kind, gen);
    }
  }
}

REGISTER_TEST(KaratsubaUnbalanced, DispatchBandWorstCases)
{
  // The shapes that overflowed the workspace before the slicing path.
  std::mt19937_64 gen(0x0FF10ADULL);
  Check(1214, 65, 0, gen);
  Check(1214, 65, 1, gen);
  Check(1000, 70, 0, gen);
  Check(1100, 66, 2, gen);
  Check(1279, 33, 0, gen);
}
