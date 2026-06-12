/**
 * BigMath: Division dispatcher
 *
 * Dispatch order:
 *   1. NewtonDivision (blockwise handles arbitrary ratio via reciprocal cache),
 *      when any of these skew bands hold:
 *        - b ≥ NEWTON_MEDIUM_B    AND  a ≥ NEWTON_SKEW (5/2)          · b
 *        - b ≥ NEWTON_RATIO2_B    AND  a ≥ NEWTON_RATIO2 (8/5)        · b
 *        - b ≥ NEWTON_BALANCED_B  AND  a ≥ NEWTON_BALANCED (4/3)      · b
 *        - b ≥ NEWTON_HIGH_SKEW_B AND  a ≥ NEWTON_HIGH_SKEW (8/1)     · b
 *      The balanced (ratio ≥ 4/3) band starts at 24k limbs — the generic
 *      Newton/BZ crossover measured after the wraparound-Newton PRs (#85-#87);
 *      below it BZ wins near-balanced, above it BZ degrades erratically.
 *   2. QuotientSizedDivision when b ≥ NEWTON_BALANCED_B, a ≥ b + 64, and
 *      ratio < 4/3 — short-quotient shapes where cost should scale with the
 *      quotient, not the divisor (and where BZ blows up 7-128× on
 *      2^k+1-family divisor sizes).
 *   3. Power-of-two base  AND  b > BZ_DIVISOR_THRESHOLD  AND  BZ band fits
 *      → BurnikelZieglerDivision    (balanced 2n/n recursion)
 *   4. else
 *      → FastDivision        (Knuth Algorithm D, hybrid-64-bit basecase)
 *   5. single-limb divisor inside the above → ClassicDivision
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
#include "division/QuotientSizedDivision.h"

namespace BigMath
{
#ifndef BIGMATH_NEWTON_MEDIUM_B
#define BIGMATH_NEWTON_MEDIUM_B 2560
#endif

// Ratio-≥8/5 band. Sits between the medium (3/1) and balanced (4/3) bands:
// post-#85-87 Newton beats BZ at ratio 2 from ~6k limbs (measured 9.0 vs
// 9.9 ms at 6000, 13 vs 22 ms at 12000, 24 vs 49 ms at 24000) and ties from
// ratio ~1.6. 8/5 rather than 2/1 because digit-derived operands sit at limb
// ratios like 2.0000 ± 1 limb and BZ blows up 8-12× on non-pow2 divisor
// sizes right across that edge.
#ifndef BIGMATH_NEWTON_RATIO2_B
#define BIGMATH_NEWTON_RATIO2_B 6144
#endif
#ifndef BIGMATH_NEWTON_RATIO2_NUMERATOR
#define BIGMATH_NEWTON_RATIO2_NUMERATOR 8
#endif
#ifndef BIGMATH_NEWTON_RATIO2_DENOMINATOR
#define BIGMATH_NEWTON_RATIO2_DENOMINATOR 5
#endif

#ifndef BIGMATH_BZ_DIVISOR_THRESHOLD
#define BIGMATH_BZ_DIVISOR_THRESHOLD 512
#endif

#ifndef BIGMATH_NEWTON_SKEW_NUMERATOR
#define BIGMATH_NEWTON_SKEW_NUMERATOR 5
#endif

#ifndef BIGMATH_NEWTON_SKEW_DENOMINATOR
#define BIGMATH_NEWTON_SKEW_DENOMINATOR 2
#endif

#ifndef BIGMATH_NEWTON_BALANCED_B
#define BIGMATH_NEWTON_BALANCED_B 24576
#endif

#ifndef BIGMATH_NEWTON_BALANCED_NUMERATOR
#define BIGMATH_NEWTON_BALANCED_NUMERATOR 4
#endif

#ifndef BIGMATH_NEWTON_BALANCED_DENOMINATOR
#define BIGMATH_NEWTON_BALANCED_DENOMINATOR 3
#endif

#ifndef BIGMATH_QSIZED_MIN_DELTA
#define BIGMATH_QSIZED_MIN_DELTA 64
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
  extern const SizeT NEWTON_RATIO2_B;
  extern const SizeT NEWTON_RATIO2_NUMERATOR;
  extern const SizeT NEWTON_RATIO2_DENOMINATOR;
  extern const SizeT BZ_DIVISOR_THRESHOLD;
  extern const SizeT NEWTON_SKEW_NUMERATOR;
  extern const SizeT NEWTON_SKEW_DENOMINATOR;
  extern const SizeT NEWTON_BALANCED_B;
  extern const SizeT NEWTON_BALANCED_NUMERATOR;
  extern const SizeT NEWTON_BALANCED_DENOMINATOR;
  extern const SizeT QSIZED_MIN_DELTA;
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
