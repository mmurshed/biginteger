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

// Newton frontier re-swept 2026-06-12 a second time after the BZ basecase
// retune (BIGMATH_BZ_RECURSION_THRESHOLD 512 -> 128) made BZ another
// 1.2-1.5× faster across its whole band: every Newton floor moves up.
// Frontier (floor, min ratio): (896, 8/1), (1280, 7/2), (1792, 14/5),
// (2560, 5/2), (4096, 2/1), (8192, 8/5), (131072, 4/3) — the balanced
// band macros live further down this file. The old exact-ratio knife-edge
// hazard (BZ collapsing on the wrong side of 2.0000/3.0000 ± 1 limb) died
// with the odd-size padding fix — both sides of every knife are now
// well-behaved, so the cuts below track measured crossovers only.
#ifndef BIGMATH_NEWTON_MEDIUM_B
#define BIGMATH_NEWTON_MEDIUM_B 1792
#endif

// Ratio ≥ 2 from 4096: Newton edges BZ at 3584 (tie) and wins 8-14% from
// 4096-6144. Below, BZ.
#ifndef BIGMATH_NEWTON_RATIO20_B
#define BIGMATH_NEWTON_RATIO20_B 4096
#endif
#ifndef BIGMATH_NEWTON_RATIO20_NUMERATOR
#define BIGMATH_NEWTON_RATIO20_NUMERATOR 2
#endif
#ifndef BIGMATH_NEWTON_RATIO20_DENOMINATOR
#define BIGMATH_NEWTON_RATIO20_DENOMINATOR 1
#endif

// Ratio ≥ 8/5 from 8192: BZ wins 1.6-band shapes at 4096-7168 (7-19%),
// Newton from ~8192-9216.
#ifndef BIGMATH_NEWTON_RATIO2_B
#define BIGMATH_NEWTON_RATIO2_B 8192
#endif
#ifndef BIGMATH_NEWTON_RATIO2_NUMERATOR
#define BIGMATH_NEWTON_RATIO2_NUMERATOR 8
#endif
#ifndef BIGMATH_NEWTON_RATIO2_DENOMINATOR
#define BIGMATH_NEWTON_RATIO2_DENOMINATOR 5
#endif

// 14/5 keeps the ratio-3.0000 ± 1 limb digit-derived edge inside one band
// (both sides Newton at the floor). Newton from 1792 at ratio 2.8-3.0
// (17-21%); 1543 is a tie, below BZ wins.
#ifndef BIGMATH_NEWTON_SKEW_NUMERATOR
#define BIGMATH_NEWTON_SKEW_NUMERATOR 14
#endif

#ifndef BIGMATH_NEWTON_SKEW_DENOMINATOR
#define BIGMATH_NEWTON_SKEW_DENOMINATOR 5
#endif

// Ratio ≥ 7/2 from 1280: Newton wins ratio 4-5 at 1280-1408 (15-41%);
// 1024 is BZ. 7/2 keeps the ratio-4.0000 ± 1 limb class on the Newton side.
#ifndef BIGMATH_NEWTON_RATIO35_B
#define BIGMATH_NEWTON_RATIO35_B 1280
#endif
#ifndef BIGMATH_NEWTON_RATIO35_NUMERATOR
#define BIGMATH_NEWTON_RATIO35_NUMERATOR 7
#endif
#ifndef BIGMATH_NEWTON_RATIO35_DENOMINATOR
#define BIGMATH_NEWTON_RATIO35_DENOMINATOR 2
#endif

// Exactly-2.5 shapes (digit-derived 2.5000 ± 1 limb): Newton from ~2560
// (2816: 18%); 2304 and below, BZ.
#ifndef BIGMATH_NEWTON_MID_B
#define BIGMATH_NEWTON_MID_B 2560
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

// 896 (06-12 post-basecase-retune probe): Newton wins ratio 8 from 896
// (12%) and 1024 (21%); 768 stays BZ (7%) and Newton degrades sharply at
// b ≈ 520, so the floor must not go below ~832.
#ifndef BIGMATH_NEWTON_HIGH_SKEW_B
#define BIGMATH_NEWTON_HIGH_SKEW_B 896
#endif

#ifndef BIGMATH_NEWTON_HIGH_SKEW_NUMERATOR
#define BIGMATH_NEWTON_HIGH_SKEW_NUMERATOR 8
#endif

#ifndef BIGMATH_NEWTON_HIGH_SKEW_DENOMINATOR
#define BIGMATH_NEWTON_HIGH_SKEW_DENOMINATOR 1
#endif

// Balanced band (ratio >= 4/3): Newton from 131072 limbs — raised from
// 24576 on 2026-06-12 after the BZ odd-size padding fix + 128-limb basecase
// retune made padded BZ win ratio 1.4-1.5 through ~98k limbs.
#ifndef BIGMATH_NEWTON_BALANCED_B
#define BIGMATH_NEWTON_BALANCED_B 131072
#endif
#ifndef BIGMATH_NEWTON_BALANCED_NUMERATOR
#define BIGMATH_NEWTON_BALANCED_NUMERATOR 4
#endif
#ifndef BIGMATH_NEWTON_BALANCED_DENOMINATOR
#define BIGMATH_NEWTON_BALANCED_DENOMINATOR 3
#endif

// Quotient-sized band: b >= 24576 (decoupled from the balanced floor),
// a >= b + BIGMATH_QSIZED_MIN_DELTA, ratio < 4/3.
#ifndef BIGMATH_QSIZED_MAIN_B
#define BIGMATH_QSIZED_MAIN_B 24576
#endif
#ifndef BIGMATH_QSIZED_MIN_DELTA
#define BIGMATH_QSIZED_MIN_DELTA 64
#endif

// Burnikel-Ziegler recursion basecase (retuned 512 -> 128 on 2026-06-12:
// Karatsuba-backed multiplies beat 512-limb Knuth-D basecase calls).
#ifndef BIGMATH_BZ_RECURSION_THRESHOLD
#define BIGMATH_BZ_RECURSION_THRESHOLD 128
#endif

// Newton division: minimum combined size for the cyclic (mod 2^k-1)
// wrapped-remainder multiply.
#ifndef BIGMATH_CYCLIC_NTT_THRESHOLD
#define BIGMATH_CYCLIC_NTT_THRESHOLD 1280
#endif

#endif
