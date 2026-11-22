#ifndef BIGMATH_DISPATCH_THRESHOLDS
#define BIGMATH_DISPATCH_THRESHOLDS

// Centralized dispatch threshold defaults.
//
// This header is intentionally kept as a build-time override point: the local
// dispatch tuner can emit a profile for a specific compiler / CPU combination,
// while the checked-in values remain the conservative defaults used by the
// repo. Each macro is only defined when the build has not supplied its own
// -DBIGMATH_* override.

// Re-swept twice on 2026-11-20. 96 -> 64: the Karatsuba leaf cutoff had moved
// 48 -> 32 without this following it, so balanced 33-48-limb products stayed
// on schoolbook while Karatsuba was already 4-27% faster there. 64 -> 56 with
// the flag-chain limb kernels and the leaf cutoff at 24: Karatsuba now wins
// from 32x32 (519 vs 556 ns), schoolbook still wins at 28x28 (433 vs 452 ns).
#ifndef BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD 56
#else
#define BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD 0
#endif
#endif

// Multiplication dispatch, re-swept 2026-11-20 over a (min, max) shape grid
// after the limb kernels, the Karatsuba slicing path and the NTT kernel /
// threading / truncated-length work (M1 Max):
//
//   - Schoolbook is best exactly when the shorter operand is at or below the
//     Karatsuba leaf cutoff (24): Karatsuba would hand the whole product to
//     the same kernel anyway. For min 32-64 sliced Karatsuba now beats
//     schoolbook by 7-25% at every ratio, so the old "min <= 64 and
//     max >= 10*min -> classic" skew guard is retired (threshold 0).
//   - NTT beats Karatsuba when the shorter operand has at least 64 limbs and
//     the schoolbook area la*lb reaches ~96k limb^2: balanced from 320 limbs
//     per operand including the total-size floor (was 640), 128 x 768, 64 x ~1500.
//     Below min 64 sliced Karatsuba wins at every measured length (48 x 6144: 139 vs 146 us).
//   - BIGMATH_NTT_MULTIPLICATION_THRESHOLD stays as a floor on the total.
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

// The Toom-3 window is [TOOM3, NTT) on the total; equal values keep it
// retired (CRT NTT beats Toom-3 everywhere it used to win).
#ifndef BIGMATH_TOOM3_MULTIPLICATION_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_TOOM3_MULTIPLICATION_THRESHOLD 640
#else
#define BIGMATH_TOOM3_MULTIPLICATION_THRESHOLD 1280
#endif
#endif

#ifndef BIGMATH_TOOM3_SKEW_RATIO
#define BIGMATH_TOOM3_SKEW_RATIO 2
#endif

#ifndef BIGMATH_NTT_MULTIPLICATION_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_NTT_MULTIPLICATION_THRESHOLD 640
#else
#define BIGMATH_NTT_MULTIPLICATION_THRESHOLD 1280
#endif
#endif

// NTT additionally requires min(la, lb) >= NTT_MIN_LIMB and
// la * lb >= NTT_AREA (limb^2). 0 / 0 reduces the rule to the total alone.
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

// The area threshold above was measured in a tight loop, where the thread
// pool stays awake between multiplies. An isolated multiply finds the workers
// parked and pays about the serial cost — at 320-450 limbs per operand ~52 us
// against Karatsuba's 24-43 us (and 24 us for the same NTT with a warm pool).
// So between NTT_AREA and NTT_COLD_AREA the dispatcher takes the NTT only
// when ParallelWarm() says the pool is awake; otherwise it runs Karatsuba and
// pokes the pool awake, so a second multiply right behind this one goes to
// the NTT. From NTT_COLD_AREA up (cold NTT 53 us vs Karatsuba 50 us at
// 512 x 512) the NTT wins regardless. Serial builds are always "cold".
#ifndef BIGMATH_NTT_COLD_AREA_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_NTT_COLD_AREA_THRESHOLD 262144
#else
#define BIGMATH_NTT_COLD_AREA_THRESHOLD 0
#endif
#endif

// CRT+NEON beats single-prime Goldilocks at every measured size (512/op:
// 0.106 vs 0.176 ms), so the CRT gate drops below the NTT entry point —
// Goldilocks remains only as the non-aarch64-friendly fallback via
// -DBIGMATH_NTT_CRT=0.
#ifndef BIGMATH_NTT_CRT_THRESHOLD
#define BIGMATH_NTT_CRT_THRESHOLD 256
#endif

// Wider CRT primes for products beyond the standard residue limits.
#ifndef BIGMATH_NTT_WIDE
#define BIGMATH_NTT_WIDE 1
#endif
#ifndef BIGMATH_NTT_WIDE_MAX_COEFFS
#define BIGMATH_NTT_WIDE_MAX_COEFFS 1811939328ULL
#endif
#ifndef BIGMATH_NTT_WIDE_MAX_TERMS
#define BIGMATH_NTT_WIDE_MAX_TERMS 836075521ULL
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

// On balanced 64-bit products, the warm NTT entry starts at 320 limbs.
// Below it an identical operand can use symmetric squaring directly.
#ifndef BIGMATH_SMALL_SELF_SQUARE_MAX
#define BIGMATH_SMALL_SELF_SQUARE_MAX 319
#endif

// Retuned 2048 -> 640 on 2026-06-12 when the square dispatch moved from the
// single-prime Goldilocks NTTSquare to the CRT+NEON NttCrt self-multiply, and
// 640 -> 352 on 2026-11-20 after the NTT kernel/threading work: KaratsubaSquare
// ties the CRT self-multiply at 352 limbs (19.5 vs 19.4 us) and loses 1.2x at
// 384, 1.8x at 512.
// Toom-Cook splits inside the Karatsuba recursion (64-bit limbs; swept
// 2026-11-21 on Apple M1 Max, header-only probes of KaratsubaMultiplication /
// KaratsubaSquare). A level whose shorter operand has at least this many
// limbs splits in three (Toom-3) instead of two:
//   products: neutral at 144-160 limbs, -3.5% at 168, -7% at 176-224,
//             -9% at 384, -7% at 448-512, -10% at 768 (256 flat);
//   squares (threshold 160; 128 was +2.4% at 128 limbs in a public-API
//             A/B): neutral at 96-160, -4.5% at 192, -6% at 320, -9% at 384,
//             -8.5% at 512, -9% at 608.
// Toom-4 is implemented but off: below the NTT entry it was -2.5% at 320
// limbs, +5% at 384 and flat at 448-512 against Toom-3, and balanced
// products past 512 limbs go to the NTT. Platforms whose NTT starts later
// may want it (try ~500).
#ifndef BIGMATH_MUL_TOOM33_THRESHOLD
#define BIGMATH_MUL_TOOM33_THRESHOLD 168
#endif
#ifndef BIGMATH_MUL_TOOM44_THRESHOLD
#define BIGMATH_MUL_TOOM44_THRESHOLD 1000000
#endif
#ifndef BIGMATH_SQR_TOOM3_THRESHOLD
#define BIGMATH_SQR_TOOM3_THRESHOLD 160
#endif

#ifndef BIGMATH_NTT_SQUARE_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_NTT_SQUARE_THRESHOLD 352
#else
#define BIGMATH_NTT_SQUARE_THRESHOLD 512
#endif
#endif

// Same pool-temperature split as BIGMATH_NTT_COLD_AREA_THRESHOLD: between the
// two square thresholds the CRT self-multiply is taken only while the pool is
// warm (cold it costs 41-47 us against KaratsubaSquare's 19.5 us at 352 limbs,
// 36 us at 512, 43 us at 576); from here up it wins cold as well (640 limbs:
// 47 vs 51 us).
#ifndef BIGMATH_NTT_SQUARE_COLD_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_NTT_SQUARE_COLD_THRESHOLD 608
#else
#define BIGMATH_NTT_SQUARE_COLD_THRESHOLD 0
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

// ─── 64-bit-limb division dispatch (re-swept 2026-11-21) ─────────────────────
// Measured on a (divisor limbs 32..262144) x (ratio 1.02..10) grid after the
// raw-buffer basecase, the in-place Burnikel-Ziegler and the Newton reciprocal
// rework, each algorithm timed directly (M1 Max). With nb = divisor limbs,
// delta = dividend limbs - nb (the quotient has delta + 1 limbs):
//
//   nb < DIV_BZ_MIN_B                     -> FastDivision (schoolbook)
//   delta >= nb (ratio >= 2):
//     nb >= DIV_NEWTON_MIN_B, or
//     nb >= DIV_NEWTON_SKEW_MIN_B and ratio >= DIV_NEWTON_SKEW_RATIO
//                                         -> Newton
//     otherwise                           -> Burnikel-Ziegler
//   delta < nb (quotient shorter than the divisor):
//     delta >= DIV_QSIZED_MIN_QUOTIENT    -> quotient-sized (Newton on the tops)
//     nb >= DIV_NEWTON_NEAR2_MIN_B and ratio >= 7/4
//                                         -> Newton
//     otherwise                           -> Burnikel-Ziegler
//
// In-place BZ ties schoolbook below ~96 limbs (it bottoms out in the same
// kernel) and beats it from there, at every ratio. For a short quotient BZ
// already has the quotient-sized cost shape (recursive divide of the tops plus
// one delta x nb product), so the old thin-quotient band is gone; quotient-
// sized division proper wins once the quotient is long enough (~3k limbs) for
// a Newton tops-divide to beat BZ's recursion (delta 3072: 0.78-0.87x of the
// better of BZ/Newton; delta 2703: 1.11x). At ratio ~1.75 it trades places
// with Newton within +-20% depending on whether nb is a power of two. The legacy 32-in-64 layout
// keeps the NEWTON_* / QSIZED_* / BZ_DIVISOR_THRESHOLD bands above.
#ifndef BIGMATH_DIV_BZ_MIN_B
#define BIGMATH_DIV_BZ_MIN_B 96
#endif
#ifndef BIGMATH_DIV_NEWTON_MIN_B
#define BIGMATH_DIV_NEWTON_MIN_B 1536
#endif
#ifndef BIGMATH_DIV_NEWTON_SKEW_MIN_B
#define BIGMATH_DIV_NEWTON_SKEW_MIN_B 768
#endif
#ifndef BIGMATH_DIV_NEWTON_SKEW_RATIO
#define BIGMATH_DIV_NEWTON_SKEW_RATIO 3
#endif
#ifndef BIGMATH_DIV_NEWTON_NEAR2_MIN_B
#define BIGMATH_DIV_NEWTON_NEAR2_MIN_B 3072
#endif
#ifndef BIGMATH_DIV_QSIZED_MIN_QUOTIENT
#define BIGMATH_DIV_QSIZED_MIN_QUOTIENT 3072
#endif

#endif
