/**
 * BigMath: Multiplication dispatcher implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/Multiplication.h"

#include <algorithm>

namespace BigMath
{
  const SizeT CLASSIC_MULTIPLICATION_THRESHOLD = BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD;
  const SizeT TOOM3_MULTIPLICATION_THRESHOLD = BIGMATH_TOOM3_MULTIPLICATION_THRESHOLD;
  const SizeT NTT_MULTIPLICATION_THRESHOLD = BIGMATH_NTT_MULTIPLICATION_THRESHOLD;
  const SizeT TOOM3_SKEW_RATIO = BIGMATH_TOOM3_SKEW_RATIO;
  const SizeT CLASSIC_MIN_LIMB_THRESHOLD = BIGMATH_CLASSIC_MIN_LIMB_THRESHOLD;
  const SizeT CLASSIC_SKEW_MIN_LIMB_THRESHOLD = BIGMATH_CLASSIC_SKEW_MIN_LIMB_THRESHOLD;
  const SizeT CLASSIC_SKEW_RATIO = BIGMATH_CLASSIC_SKEW_RATIO;

  std::vector<DataT> Multiply(std::vector<DataT> const &a,
                              std::vector<DataT> const &b,
                              BaseT base)
  {
    if (IsZero(a) || IsZero(b))
      return std::vector<DataT>{0};

    if (b.size() == 1)
      return ClassicMultiplication::Multiply(a, b[0], base);
    if (a.size() == 1)
      return ClassicMultiplication::Multiply(b, a[0], base);

    SizeT size = a.size() + b.size();
    SizeT minSize = std::min(a.size(), b.size());
    SizeT maxSize = std::max(a.size(), b.size());

    if (size <= CLASSIC_MULTIPLICATION_THRESHOLD || minSize <= CLASSIC_MIN_LIMB_THRESHOLD)
      return ClassicMultiplication::Multiply(a, b, base);

    if (minSize <= CLASSIC_SKEW_MIN_LIMB_THRESHOLD && maxSize >= CLASSIC_SKEW_RATIO * minSize)
      return ClassicMultiplication::Multiply(a, b, base);

    if (size < TOOM3_MULTIPLICATION_THRESHOLD)
      return KaratsubaMultiplication::Multiply(a, b, base);

    if (size < NTT_MULTIPLICATION_THRESHOLD && maxSize < TOOM3_SKEW_RATIO * minSize)
      return ToomCookMultiplication::Multiply(a, b, base);

    if (size < NTT_MULTIPLICATION_THRESHOLD)
      return KaratsubaMultiplication::Multiply(a, b, base);

    return NTTMultiplication::Multiply(a, b, base);
  }

  std::vector<DataT> Multiply(std::vector<DataT> const &a,
                              DataT b,
                              BaseT base)
  {
    if (b == 0 || IsZero(a))
      return std::vector<DataT>{0};
    return ClassicMultiplication::Multiply(a, b, base);
  }
}
