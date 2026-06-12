/**
 * BigMath: Multiplication dispatcher
 *
 * Picks the best strategy based on operand sizes:
 *   - single-limb operand                          → ClassicMultiplication (scalar)
 *   - sum ≤ CLASSIC_MULTIPLICATION_THRESHOLD OR
 *     min ≤ CLASSIC_MIN_LIMB_THRESHOLD             → ClassicMultiplication
 *   - sum < TOOM3_MULTIPLICATION_THRESHOLD         → KaratsubaMultiplication
 *   - sum < NTT_MULTIPLICATION_THRESHOLD and
 *     max < TOOM3_SKEW_RATIO · min                 → ToomCookMultiplication (Toom-3)
 *   - sum < NTT_MULTIPLICATION_THRESHOLD           → KaratsubaMultiplication
 *   - otherwise                                    → NTTMultiplication
 *
 * Toom-3's former dispatch window (total ≈ 2560-5120 limbs) was retired on
 * 2026-06-12: with TOOM3_MULTIPLICATION_THRESHOLD == NTT_MULTIPLICATION_THRESHOLD
 * the Toom-3 branch below is unreachable and CRT+NEON NTT covers the band
 * (see DispatchThresholds.h). ToomCookMultiplication stays as a
 * test-exercised alternate. Toom-5 has no productive band — it ties
 * Karatsuba below 256 per-operand limbs and degrades fast above, so it
 * stays excluded from dispatch.
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
// dispatch_tuner LIMB_64=1: Classic wins up to total 96 limbs (much larger
// schoolbook range than 32-bit, where Classic only wins at 1-limb operands).
#if BIGMATH_LIMB_64
#define BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD 96
#else
#define BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD 0
#endif
#endif

#ifndef BIGMATH_TOOM3_MULTIPLICATION_THRESHOLD
// Set equal to the NTT threshold since 2026-06-12, which retires the Toom-3
// dispatch window (CRT+NEON NTT beats it everywhere it used to win).
#define BIGMATH_TOOM3_MULTIPLICATION_THRESHOLD 1280
#endif

#ifndef BIGMATH_NTT_MULTIPLICATION_THRESHOLD
// NTT crossover. Dropped 5120 -> 1280 on 2026-06-12 after the NEON Shoup
// butterflies made CRT NTT ~2.5x faster in the small-mid band.
#define BIGMATH_NTT_MULTIPLICATION_THRESHOLD 1280
#endif

#ifndef BIGMATH_TOOM3_SKEW_RATIO
// Toom-3 is a balanced-product win. For 2:1+ skew inside the Toom window,
// Karatsuba is consistently faster on the measured AppleClang/native build.
#define BIGMATH_TOOM3_SKEW_RATIO 2
#endif

#ifndef BIGMATH_CLASSIC_MIN_LIMB_THRESHOLD
#define BIGMATH_CLASSIC_MIN_LIMB_THRESHOLD 0
#endif

#ifndef BIGMATH_CLASSIC_SKEW_MIN_LIMB_THRESHOLD
#define BIGMATH_CLASSIC_SKEW_MIN_LIMB_THRESHOLD 64
#endif

#ifndef BIGMATH_CLASSIC_SKEW_RATIO
#define BIGMATH_CLASSIC_SKEW_RATIO 10
#endif

  extern const SizeT CLASSIC_MULTIPLICATION_THRESHOLD;
  extern const SizeT TOOM3_MULTIPLICATION_THRESHOLD;
  extern const SizeT NTT_MULTIPLICATION_THRESHOLD;
  extern const SizeT TOOM3_SKEW_RATIO;
  extern const SizeT CLASSIC_MIN_LIMB_THRESHOLD;
  extern const SizeT CLASSIC_SKEW_MIN_LIMB_THRESHOLD;
  extern const SizeT CLASSIC_SKEW_RATIO;

  std::vector<DataT> Multiply(std::vector<DataT> const &a,
                              std::vector<DataT> const &b,
                              BaseT base);

  std::vector<DataT> Multiply(std::vector<DataT> const &a,
                              DataT b,
                              BaseT base);
}

#endif
