/**
 * BigMath: Multiplication dispatcher
 *
 * Picks the best strategy based on operand sizes (64-bit-limb defaults):
 *   - single-limb operand                          → ClassicMultiplication (scalar)
 *   - same operand, 2 through 319 limbs             → Square
 *   - sum ≤ CLASSIC_MULTIPLICATION_THRESHOLD (56) OR
 *     min ≤ CLASSIC_MIN_LIMB_THRESHOLD (24)        → ClassicMultiplication
 *   - sum ≥ NTT_MULTIPLICATION_THRESHOLD (640) AND
 *     min ≥ NTT_MIN_LIMB_THRESHOLD (64) AND
 *     la·lb ≥ NTT_AREA_THRESHOLD (96000) AND
 *     (la·lb ≥ NTT_COLD_AREA_THRESHOLD (262144)
 *      OR the thread pool is warm)                 → NTTMultiplication
 *   - cold NTT fallback for the same operand       → KaratsubaSquare
 *   - otherwise                                    → KaratsubaMultiplication
 *                                                    (slices unbalanced operands)
 *
 * The NTT condition is an area test rather than a plain total because the
 * alternative for a skewed product is sliced Karatsuba, whose cost follows
 * max · K(min)/min: a short operand keeps Karatsuba ahead far beyond the
 * balanced crossover (48 × 6144 limbs still favors it).
 *
 * Karatsuba selects the limb-based Toom splits inside its recursion.
 * The 64-bit defaults disable the legacy schoolbook skew guard
 * (CLASSIC_SKEW_MIN_LIMB_THRESHOLD = 0) because Karatsuba slices skewed inputs.
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
#include "../algorithms/multiplication/NTTMultiplication.h"

namespace BigMath
{
#ifndef BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD
// LIMB_64=1: Classic wins up to total 56 limbs (28x28); from 32x32 Karatsuba
// is faster. Must stay value-identical to build/DispatchThresholds.h.
#if BIGMATH_LIMB_64
#define BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD 56
#else
#define BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD 0
#endif
#endif

// The #ifndef fallbacks below must stay value-identical to
// build/DispatchThresholds.h (the canonical home, included via Constants.h);
// they only apply when that header is bypassed.
#ifndef BIGMATH_NTT_MULTIPLICATION_THRESHOLD
// Floor on a.size() + b.size() for NTT. 5120 -> 1280 on 2026-06-12, -> 640 on
// 2026-11-20 (together with the min-limb and area conditions below).
#if BIGMATH_LIMB_64
#define BIGMATH_NTT_MULTIPLICATION_THRESHOLD 640
#else
#define BIGMATH_NTT_MULTIPLICATION_THRESHOLD 1280
#endif
#endif

#ifndef BIGMATH_NTT_MIN_LIMB_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_NTT_MIN_LIMB_THRESHOLD 64
#else
#define BIGMATH_NTT_MIN_LIMB_THRESHOLD 0
#endif
#endif

#ifndef BIGMATH_NTT_AREA_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_NTT_AREA_THRESHOLD 96000
#else
#define BIGMATH_NTT_AREA_THRESHOLD 0
#endif
#endif

#ifndef BIGMATH_NTT_COLD_AREA_THRESHOLD
// Below this area the NTT is taken only while the thread pool is warm.
#if BIGMATH_LIMB_64
#define BIGMATH_NTT_COLD_AREA_THRESHOLD 262144
#else
#define BIGMATH_NTT_COLD_AREA_THRESHOLD 0
#endif
#endif

#ifndef BIGMATH_CLASSIC_MIN_LIMB_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_CLASSIC_MIN_LIMB_THRESHOLD 24
#else
#define BIGMATH_CLASSIC_MIN_LIMB_THRESHOLD 0
#endif
#endif

#ifndef BIGMATH_CLASSIC_SKEW_MIN_LIMB_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_CLASSIC_SKEW_MIN_LIMB_THRESHOLD 0
#else
#define BIGMATH_CLASSIC_SKEW_MIN_LIMB_THRESHOLD 64
#endif
#endif

#ifndef BIGMATH_CLASSIC_SKEW_RATIO
#define BIGMATH_CLASSIC_SKEW_RATIO 10
#endif

  extern const SizeT CLASSIC_MULTIPLICATION_THRESHOLD;
  extern const SizeT NTT_MULTIPLICATION_THRESHOLD;
  extern const SizeT NTT_MIN_LIMB_THRESHOLD;
  extern const ULong NTT_AREA_THRESHOLD;
  extern const ULong NTT_COLD_AREA_THRESHOLD;
  extern const SizeT CLASSIC_MIN_LIMB_THRESHOLD;
  extern const SizeT CLASSIC_SKEW_MIN_LIMB_THRESHOLD;
  extern const SizeT CLASSIC_SKEW_RATIO;

  namespace detail
  {
    inline bool ClassicProductShape(SizeT shorter, SizeT longer)
    {
      return (ULong)shorter + longer <= CLASSIC_MULTIPLICATION_THRESHOLD ||
             shorter <= CLASSIC_MIN_LIMB_THRESHOLD ||
             (shorter <= CLASSIC_SKEW_MIN_LIMB_THRESHOLD && longer >= (ULong)CLASSIC_SKEW_RATIO * shorter);
    }
  }

  std::vector<DataT> Multiply(std::vector<DataT> const &a,
                              std::vector<DataT> const &b,
                              BaseT base);

  std::vector<DataT> Multiply(std::vector<DataT> const &a,
                              DataT b,
                              BaseT base);
}

#endif
