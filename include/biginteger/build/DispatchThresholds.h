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

// Post-NEON (PR #96) retune: CRT NTT with NEON Shoup butterflies beats
// Karatsuba from ~384 limbs per operand (sum 768) and Toom-3 everywhere,
// so the Toom-3 window is retired from dispatch (kept as a cross-check
// algorithm) and the NTT entry drops 5120 -> 1280 (sum 896-1024 also favors CRT but sits one bit_ceil tier below sum 1040-1250 where Karatsuba still edges it - 1280 keeps the dispatch cliff-free).
#ifndef BIGMATH_TOOM3_MULTIPLICATION_THRESHOLD
#define BIGMATH_TOOM3_MULTIPLICATION_THRESHOLD 1280
#endif

#ifndef BIGMATH_TOOM3_SKEW_RATIO
#define BIGMATH_TOOM3_SKEW_RATIO 2
#endif

#ifndef BIGMATH_NTT_MULTIPLICATION_THRESHOLD
#define BIGMATH_NTT_MULTIPLICATION_THRESHOLD 1280
#endif

// CRT+NEON beats single-prime Goldilocks at every measured size (512/op:
// 0.106 vs 0.176 ms), so the CRT gate drops below the NTT entry point —
// Goldilocks remains only as the non-aarch64-friendly fallback via
// -DBIGMATH_NTT_CRT=0.
#ifndef BIGMATH_NTT_CRT_THRESHOLD
#define BIGMATH_NTT_CRT_THRESHOLD 256
#endif

// Post-PR-#107 retune (2026-06-12): the row-chunked ParallelDo(6) fused MFA
// stages flipped the old 2^24 break-even — the whole-transform non-MFA path
// idles cores (6-unit forward, 3-unit inverse) while the fused path keeps
// 6-12 units busy. Warm-state sweep (quiet M1 Max, best-of-3 interleaved):
// n=2^22 mul 3.0× faster, n=2^23 2.2×, div 100M÷20M digits 1.37→0.74× vs
// GMP, 200M÷40M 1.47→1.01×. 2^18 measured ≈ wash vs 2^20 (≤4% div); 2^20
// keeps the gate out of the latency-sensitive sub-ms band.
#ifndef BIGMATH_NTT_MFA_THRESHOLD
#define BIGMATH_NTT_MFA_THRESHOLD (1 << 20)
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
#define BIGMATH_NEWTON_MEDIUM_B 1280
#endif

// Newton floors re-swept 2026-06-12 (smallskew_div_plan.md S1) AFTER fixing
// the BZ odd-size FastDivision fallback (operand padding in
// BurnikelZieglerDivision.h). That fix changed every crossover: BZ used to
// collapse to O(n·Δ) on odd divisor sizes (1.3-10× losses), which had pushed
// the Newton floors artificially low as insurance. With padded BZ
// well-behaved on all sizes the floors are data-driven again — the Newton
// frontier is now (640, 8/1), (1024, 7/2), (1280, 14/5), (1792, 5/2),
// (4096, 8/5), (24576, 4/3); ratio-8/5@4096 re-confirmed unchanged.
// Ratio-≥8/5 Newton band between the mid (5/2) and balanced (4/3) bands.
// 8/5 instead of a knife-edge 2/1: digit-derived operands land at limb ratios
// like 2.0000 ± 1 limb; Newton near-ties padded BZ from ratio ~1.6 at 4096
// (at ratio 1.6 BZ wins 18-22% at b=2048-3072, Newton from 4096).
#ifndef BIGMATH_NEWTON_RATIO2_B
#define BIGMATH_NEWTON_RATIO2_B 4096
#endif
#ifndef BIGMATH_NEWTON_RATIO2_NUMERATOR
#define BIGMATH_NEWTON_RATIO2_NUMERATOR 8
#endif
#ifndef BIGMATH_NEWTON_RATIO2_DENOMINATOR
#define BIGMATH_NEWTON_RATIO2_DENOMINATOR 5
#endif

// 14/5 rather than 3/1: digit-derived operands land at limb ratios like
// 3.0000 ± 1 limb; 2.8 keeps that edge on the Newton side (Newton is still
// 2.3× better than even padded BZ at 15579×5193). Not 5/2: at exactly
// ratio 2.5, padded BZ beats Newton up to b≈1792 (06-12 probe: 7-34% at
// 1024-1543), so 2.5 has its own band below with a higher floor. Floor
// 1280: BZ wins ratio 2.8-3.0 at 1037-1163, Newton from 1291.
#ifndef BIGMATH_NEWTON_SKEW_NUMERATOR
#define BIGMATH_NEWTON_SKEW_NUMERATOR 14
#endif

#ifndef BIGMATH_NEWTON_SKEW_DENOMINATOR
#define BIGMATH_NEWTON_SKEW_DENOMINATOR 5
#endif

// Ratio ≥ 7/2 from 1024: Newton wins ratio 4-6 across all of [1024, 1280)
// (3-45%, growing with ratio); 7/2 rather than 4/1 keeps the ratio-4.0000
// ± 1 limb knife-edge (40k÷10k-digit class) on the Newton side, and 3.5
// itself is a measured tie. Below 3.2, BZ wins — that segment belongs to
// the 14/5@1280 band above.
#ifndef BIGMATH_NEWTON_RATIO35_B
#define BIGMATH_NEWTON_RATIO35_B 1024
#endif
#ifndef BIGMATH_NEWTON_RATIO35_NUMERATOR
#define BIGMATH_NEWTON_RATIO35_NUMERATOR 7
#endif
#ifndef BIGMATH_NEWTON_RATIO35_DENOMINATOR
#define BIGMATH_NEWTON_RATIO35_DENOMINATOR 2
#endif

// Exactly-2.5 shapes (digit-derived 2.5000 ± 1 limb): Newton from 1792
// (14% win at 1791, 6% at 2049, 29% at 2600); below, padded BZ wins.
#ifndef BIGMATH_NEWTON_MID_B
#define BIGMATH_NEWTON_MID_B 1792
#endif
#ifndef BIGMATH_NEWTON_MID_NUMERATOR
#define BIGMATH_NEWTON_MID_NUMERATOR 5
#endif
#ifndef BIGMATH_NEWTON_MID_DENOMINATOR
#define BIGMATH_NEWTON_MID_DENOMINATOR 2
#endif

// Thin-quotient extension of the quotient-sized band below the balanced
// floor: delta <= b/8 at b >= 8192. Generic wins are modest (8192 limbs at
// ratio 1.1: 3.8 vs BZ 5.5 ms) but the 2^k+1-family BZ pathology there is
// 3-14x (16385 limbs ratio 1.1: 11.9 vs 66.4 ms).
#ifndef BIGMATH_QSIZED_SMALL_B
#define BIGMATH_QSIZED_SMALL_B 8192
#endif
#ifndef BIGMATH_QSIZED_SMALL_DELTA_DIV
#define BIGMATH_QSIZED_SMALL_DELTA_DIV 8
#endif

// 640 (06-12 probe): ratio 8 ties BZ at 640 and wins 16% at 704; ratio 10
// wins 12% at 640. 576 stays BZ (16% better) and Newton degrades sharply
// at b ≈ 520, so the floor must not go below 640.
#ifndef BIGMATH_NEWTON_HIGH_SKEW_B
#define BIGMATH_NEWTON_HIGH_SKEW_B 640
#endif

#ifndef BIGMATH_NEWTON_HIGH_SKEW_NUMERATOR
#define BIGMATH_NEWTON_HIGH_SKEW_NUMERATOR 8
#endif

#ifndef BIGMATH_NEWTON_HIGH_SKEW_DENOMINATOR
#define BIGMATH_NEWTON_HIGH_SKEW_DENOMINATOR 1
#endif

#endif
