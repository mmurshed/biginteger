/**
 * BigMath: Division dispatcher
 *
 * Dispatch order:
 *   1. b ≥ NEWTON_MEDIUM_B  AND  a ≥ NEWTON_SKEW_NUMERATOR/DENOMINATOR · b
 *      → NewtonDivision      (blockwise handles arbitrary ratio via reciprocal cache)
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
#ifndef BIGMATH_NEWTON_LARGE_B
#define BIGMATH_NEWTON_LARGE_B 24576
#endif

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

#ifndef BIGMATH_NEWTON_HIGH_SKEW_B
#define BIGMATH_NEWTON_HIGH_SKEW_B 2048
#endif

#ifndef BIGMATH_NEWTON_HIGH_SKEW_NUMERATOR
#define BIGMATH_NEWTON_HIGH_SKEW_NUMERATOR 8
#endif

#ifndef BIGMATH_NEWTON_HIGH_SKEW_DENOMINATOR
#define BIGMATH_NEWTON_HIGH_SKEW_DENOMINATOR 1
#endif

  extern const SizeT NEWTON_LARGE_B;
  extern const SizeT NEWTON_MEDIUM_B;
  extern const SizeT BZ_DIVISOR_THRESHOLD;
  extern const SizeT NEWTON_SKEW_NUMERATOR;
  extern const SizeT NEWTON_SKEW_DENOMINATOR;
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
