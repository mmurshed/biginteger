#ifndef BIGMATH_DISPATCH_THRESHOLDS
#define BIGMATH_DISPATCH_THRESHOLDS

// Centralized dispatch threshold defaults.
//
// This header is intentionally kept as a build-time override point: the local
// dispatch tuner can regenerate it for a specific compiler / CPU combination,
// while the checked-in values remain the conservative defaults used by the
// repo. Each macro is only defined when the build has not supplied its own
// -DBIGMATH_* override.

#ifndef BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD 96
#else
#define BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD 0
#endif
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

#ifndef BIGMATH_TOOM3_MULTIPLICATION_THRESHOLD
#define BIGMATH_TOOM3_MULTIPLICATION_THRESHOLD 2560
#endif

#ifndef BIGMATH_TOOM3_SKEW_RATIO
#define BIGMATH_TOOM3_SKEW_RATIO 2
#endif

#ifndef BIGMATH_NTT_MULTIPLICATION_THRESHOLD
#define BIGMATH_NTT_MULTIPLICATION_THRESHOLD 5120
#endif

#ifndef BIGMATH_NTT_CRT_THRESHOLD
#define BIGMATH_NTT_CRT_THRESHOLD 5000
#endif

#ifndef BIGMATH_NTT_MFA_THRESHOLD
#define BIGMATH_NTT_MFA_THRESHOLD (1 << 24)
#endif

#ifndef BIGMATH_NTT_SQUARE_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_NTT_SQUARE_THRESHOLD 2048
#else
#define BIGMATH_NTT_SQUARE_THRESHOLD 512
#endif
#endif

#ifndef BIGMATH_BZ_DIVISOR_THRESHOLD
#define BIGMATH_BZ_DIVISOR_THRESHOLD 512
#endif

#ifndef BIGMATH_NEWTON_MEDIUM_B
#define BIGMATH_NEWTON_MEDIUM_B 2560
#endif

// Ratio-≥8/5 Newton band between the medium (3/1) and balanced (4/3) bands.
// 8/5 instead of a knife-edge 2/1: digit-derived operands land at limb ratios
// like 2.0000 ± 1 limb, and the BZ side of the edge blows up 8-12× on
// non-power-of-2 divisor sizes; Newton generic-ties BZ from ratio ~1.6 here.
#ifndef BIGMATH_NEWTON_RATIO2_B
#define BIGMATH_NEWTON_RATIO2_B 6144
#endif
#ifndef BIGMATH_NEWTON_RATIO2_NUMERATOR
#define BIGMATH_NEWTON_RATIO2_NUMERATOR 8
#endif
#ifndef BIGMATH_NEWTON_RATIO2_DENOMINATOR
#define BIGMATH_NEWTON_RATIO2_DENOMINATOR 5
#endif

// 5/2 rather than 3/1: digit-derived operands land at limb ratios like
// 3.0000 ± 1 limb, and BZ blows up ~7× on non-pow2 divisor sizes across
// that edge (15578×5193 limbs: BZ 69 ms vs Newton 9 ms). Newton wins from
// ratio ~2.5 at b ≥ 2560 (near-tie at 2600, clear by 4096).
#ifndef BIGMATH_NEWTON_SKEW_NUMERATOR
#define BIGMATH_NEWTON_SKEW_NUMERATOR 5
#endif

#ifndef BIGMATH_NEWTON_SKEW_DENOMINATOR
#define BIGMATH_NEWTON_SKEW_DENOMINATOR 2
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

#endif
