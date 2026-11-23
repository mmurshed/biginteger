/**
 * BigMath: Multiplication dispatcher implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/Multiplication.h"

#include <algorithm>
#include "biginteger/algorithms/Squaring.h"
#include "biginteger/common/Parallel.h"

namespace BigMath
{
  const SizeT CLASSIC_MULTIPLICATION_THRESHOLD = BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD;
  const SizeT NTT_MULTIPLICATION_THRESHOLD = BIGMATH_NTT_MULTIPLICATION_THRESHOLD;
  const SizeT NTT_MIN_LIMB_THRESHOLD = BIGMATH_NTT_MIN_LIMB_THRESHOLD;
  const ULong NTT_AREA_THRESHOLD = BIGMATH_NTT_AREA_THRESHOLD;
  const ULong NTT_COLD_AREA_THRESHOLD = BIGMATH_NTT_COLD_AREA_THRESHOLD;
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

    // Below the measured small-square ceiling we can exploit symmetry before
    // general multiplication dispatch. CRT already recognizes aliases at scale.
    if (&a == &b && base == Base2_64 && a.size() <= BIGMATH_SMALL_SELF_SQUARE_MAX)
      return Square(a, base);

    SizeT size = a.size() + b.size();
    SizeT minSize = std::min(a.size(), b.size());
    SizeT maxSize = std::max(a.size(), b.size());

    if (detail::ClassicProductShape(minSize, maxSize))
      return ClassicMultiplication::Multiply(a, b, base);

    // NTT needs a long enough product AND a long enough shorter operand: for
    // a skewed product the alternative is sliced Karatsuba, whose cost grows
    // with max * K(min)/min, so a short operand keeps Karatsuba ahead well
    // past the balanced crossover. The area test la*lb captures that.
    const ULong area = (ULong)a.size() * (ULong)b.size();
    if (size >= NTT_MULTIPLICATION_THRESHOLD &&
        minSize >= NTT_MIN_LIMB_THRESHOLD &&
        area >= NTT_AREA_THRESHOLD)
    {
#if !BIGMATH_USE_THREADS && BIGMATH_NTT_SERIAL_BALANCED_MIN_LIMBS > 0
      if (base == Base2_64 && &a != &b && maxSize < BIGMATH_NTT_SERIAL_BALANCED_MIN_LIMBS &&
          maxSize - minSize <= 1)
        return KaratsubaMultiplication::Multiply(a, b, base);
#endif
      // Small NTT products only pay off while the pool is awake: a parked
      // pool makes them cost about the serial time, which loses to Karatsuba
      // below NTT_COLD_AREA. When cold, run Karatsuba; if these multiplies
      // are arriving in a burst the pool gets woken for the next one.
      if (area >= NTT_COLD_AREA_THRESHOLD || ParallelWarmOrPrewarm())
        return NTTMultiplication::Multiply(a, b, base);
      std::vector<DataT> product = (&a == &b && base == Base2_64)
          ? KaratsubaSquare::Square(a, base)
          : KaratsubaMultiplication::Multiply(a, b, base);
      ParallelNoteColdDone();
      return product;
    }

    return KaratsubaMultiplication::Multiply(a, b, base);
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
