/**
 * BigMath: Squaring dispatcher implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/Squaring.h"

namespace BigMath
{
  const SizeT NTT_SQUARE_THRESHOLD = BIGMATH_NTT_SQUARE_THRESHOLD;

  std::vector<DataT> Square(std::vector<DataT> const &a, BaseT base)
  {
    if (IsZero(a))
      return std::vector<DataT>{0};

    if (a.size() == 1)
      return ClassicSquare::Square(a, base);

    if (a.size() < NTT_SQUARE_THRESHOLD)
      return KaratsubaSquare::Square(a, base);

    return NTTSquare::Square(a, base);
  }
}
