/**
 * BigMath: Squaring dispatcher implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/Squaring.h"
#include "biginteger/algorithms/Multiplication.h"

namespace BigMath
{
  const SizeT NTT_SQUARE_THRESHOLD = BIGMATH_NTT_SQUARE_THRESHOLD;

  std::vector<DataT> Square(std::vector<DataT> const &a, BaseT base)
  {
    if (IsZero(a))
      return std::vector<DataT>{0};

    // Squaring fast paths are Base2_32-only. For other bases (incl. Base2_64
    // under LIMB_64=1) fall back to general multiplication. The squaring
    // family can be ported per-algo later if profile demands.
    if (base != Base2_32)
      return Multiply(a, a, base);

    if (a.size() == 1)
      return ClassicSquare::Square(a, base);

    if (a.size() < NTT_SQUARE_THRESHOLD)
      return KaratsubaSquare::Square(a, base);

    return NTTSquare::Square(a, base);
  }
}
