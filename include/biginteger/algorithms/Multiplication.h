/**
 * BigMath: Multiplication dispatcher
 *
 * Picks the best strategy based on operand sizes:
 *   - single-limb operand                          → ClassicMultiplication (scalar)
 *   - sum ≤ CLASSIC_MULTIPLICATION_THRESHOLD OR
 *     min ≤ CLASSIC_MIN_LIMB_THRESHOLD             → ClassicMultiplication
 *   - sum < NTT_MULTIPLICATION_THRESHOLD           → KaratsubaMultiplication
 *   - otherwise                                    → NTTMultiplication
 *
 * Toom-Cook 3 is implemented and correctness-checked but excluded from
 * dispatch — NTT dominates above ~256 limbs and Karatsuba dominates below,
 * leaving no productive band. See MULTIPLICATION.md for details.
 *
 * Thresholds are tunable at compile time via -DBIGMATH_*=N.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef MULTIPLICATION
#define MULTIPLICATION

#include <vector>

#include "../algorithms/multiplication/ClassicMultiplication.h"
#include "../algorithms/multiplication/KaratsubaMultiplication.h"
#include "../algorithms/multiplication/ToomCookMultiplication.h"
#include "../algorithms/multiplication/NTTMultiplication.h"

namespace BigMath
{
#ifndef BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD
#define BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD 0
#endif

#ifndef BIGMATH_NTT_MULTIPLICATION_THRESHOLD
#define BIGMATH_NTT_MULTIPLICATION_THRESHOLD 4096
#endif

#ifndef BIGMATH_CLASSIC_MIN_LIMB_THRESHOLD
#define BIGMATH_CLASSIC_MIN_LIMB_THRESHOLD 0
#endif

  extern const SizeT CLASSIC_MULTIPLICATION_THRESHOLD;
  extern const SizeT NTT_MULTIPLICATION_THRESHOLD;
  extern const SizeT CLASSIC_MIN_LIMB_THRESHOLD;

  std::vector<DataT> Multiply(std::vector<DataT> const &a,
                              std::vector<DataT> const &b,
                              BaseT base);

  std::vector<DataT> Multiply(std::vector<DataT> const &a,
                              DataT b,
                              BaseT base);
}

#endif
