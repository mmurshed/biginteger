/**
 * BigMath: Squaring dispatcher implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/Squaring.h"
#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"

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

#if BIGMATH_NTT_CRT
    // CRT NTT square: NttCrt::Multiply detects the self-operand and runs 3
    // forward transforms instead of 6; NEON Shoup + MFA apply. Measured
    // 1.6-6.4x over the single-prime Goldilocks NTTSquare at 2k-2M limbs
    // even before the self-operand skip (NTTSquare uses the 16-bit split at
    // 2x the transform length, no NEON, no MFA). NTTSquare remains the
    // -DBIGMATH_NTT_CRT=0 fallback and a test-exercised alternate.
    return NttCrt::Multiply(a, a, base);
#else
    return NTTSquare::Square(a, base);
#endif
  }
}
