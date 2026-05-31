/**
 * BigMath: Division dispatcher
 *
 * Dispatch order:
 *   1. NewtonDivision (blockwise handles arbitrary ratio via reciprocal cache),
 *      when any of these skew bands hold:
 *        - b ≥ NEWTON_MEDIUM_B    AND  a ≥ NEWTON_SKEW (3/1)          · b
 *        - b ≥ NEWTON_BALANCED_B  AND  a ≥ NEWTON_BALANCED (2/1)      · b
 *        - b ≥ NEWTON_HIGH_SKEW_B AND  a ≥ NEWTON_HIGH_SKEW (8/1)     · b
 *      The balanced (ratio ≥ 2) band starts higher (96k limbs) because BZ wins
 *      near-balanced below that; above it BZ degrades erratically (measured
 *      2×–4.5× slower at b ≥ 100k) while Newton stays smooth.
 *   2. Power-of-two base  AND  b > BZ_DIVISOR_THRESHOLD  AND  BZ band fits
 *      → BurnikelZieglerDivision    (balanced 2n/n recursion)
 *   3. else
 *      → FastDivision        (Knuth Algorithm D, hybrid-64-bit basecase)
 *   4. single-limb divisor inside the above → ClassicDivision
 *
 * Thresholds tunable at compile time via -DBIGMATH_*=N.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef DIVISION
#define DIVISION

#include <utility>
#include <vector>

#include "../BigInteger.h"
#include "division/BurnikelZieglerDivision.h"
#include "division/ClassicDivision.h"
#include "division/FastDivision.h"
#include "division/NewtonDivision.h"

namespace BigMath
{
#ifndef BIGMATH_NEWTON_MEDIUM_B
#define BIGMATH_NEWTON_MEDIUM_B 4096
#endif

#ifndef BIGMATH_BZ_DIVISOR_THRESHOLD
#define BIGMATH_BZ_DIVISOR_THRESHOLD 512
#endif

#ifndef BIGMATH_NEWTON_SKEW_NUMERATOR
#define BIGMATH_NEWTON_SKEW_NUMERATOR 3
#endif

#ifndef BIGMATH_NEWTON_SKEW_DENOMINATOR
#define BIGMATH_NEWTON_SKEW_DENOMINATOR 1
#endif

#ifndef BIGMATH_NEWTON_BALANCED_B
#define BIGMATH_NEWTON_BALANCED_B 98304
#endif

#ifndef BIGMATH_NEWTON_BALANCED_NUMERATOR
#define BIGMATH_NEWTON_BALANCED_NUMERATOR 2
#endif

#ifndef BIGMATH_NEWTON_BALANCED_DENOMINATOR
#define BIGMATH_NEWTON_BALANCED_DENOMINATOR 1
#endif

#ifndef BIGMATH_NEWTON_HIGH_SKEW_B
#define BIGMATH_NEWTON_HIGH_SKEW_B 2048
#endif

#ifndef BIGMATH_NEWTON_HIGH_SKEW_NUMERATOR
#define BIGMATH_NEWTON_HIGH_SKEW_NUMERATOR 8
#endif

#ifndef BIGMATH_NEWTON_HIGH_SKEW_DENOMINATOR
#define BIGMATH_NEWTON_HIGH_SKEW_DENOMINATOR 1
#endif

  extern const SizeT NEWTON_MEDIUM_B;
  extern const SizeT BZ_DIVISOR_THRESHOLD;
  extern const SizeT NEWTON_SKEW_NUMERATOR;
  extern const SizeT NEWTON_SKEW_DENOMINATOR;
  extern const SizeT NEWTON_BALANCED_B;
  extern const SizeT NEWTON_BALANCED_NUMERATOR;
  extern const SizeT NEWTON_BALANCED_DENOMINATOR;
  extern const SizeT NEWTON_HIGH_SKEW_B;
  extern const SizeT NEWTON_HIGH_SKEW_NUMERATOR;
  extern const SizeT NEWTON_HIGH_SKEW_DENOMINATOR;

  std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainder(
      std::vector<DataT> const &a,
      std::vector<DataT> const &b,
      BaseT base,
      bool computeRemainder = true);

  std::vector<DataT> Divide(std::vector<DataT> const &a,
                            std::vector<DataT> const &b,
                            BaseT base);

  std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainder(
      std::vector<DataT> const &a,
      DataT b,
      BaseT base);

  std::vector<DataT> Divide(std::vector<DataT> const &a,
                            DataT b,
                            BaseT base);
}

#endif
