/**
 * BigMath: Multi-prime CRT NTT multiplication.
 *
 * Alternate to Goldilocks NTT (NTTMultiplication.h). Three NTT-friendly
 * primes (2013265921, 469762049, 1811939329) with 32-bit coefficient
 * splitting: each 64-bit input limb → 2 coefficients (vs Goldilocks' 4).
 * Halves the NTT length at the cost of 3 parallel transforms + CRT
 * reconstruction.
 *
 * Prime selection: all three primes have 2-adic order ≥ 2^26, so the CRT
 * NTT length ceiling is 2^26 ≈ 67M coefficients per transform. Longer
 * products are recombined from several residues of that length
 * (MultiplyAboveTransform): up to seven transform lengths and a shorter
 * operand of about 1.79 billion digits. Wider CRT primes extend the
 * residue path to about 8.05 billion balanced digits before overlap-add. Earlier triple (998244353, 985661441,
 * 754974721) had P2 capped at 2^22, silently breaking correctness for
 * operand pairs above ~20M digits.
 *
 * Trade-off vs Goldilocks:
 *   - 3 transforms × NTT(N/2) work each = ~1.5× more transform layers,
 *     BUT each layer runs on 32-bit modular ops (single MUL + simple
 *     reduce) rather than ULong128 Goldilocks mul+reduce.
 *   - Coefficients fit in 32-bit lanes → better cache density.
 *   - Combined modulus = p1·p2·p3 ≈ 2^91, comfortably holds 32×32 product
 *     summed over up to 2^26 ≈ 67 M coefficients.
 *
 * Gated on BIGMATH_NTT_CRT (default ON, see below); selected at compile
 * time in NTTMultiplication's dispatch.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

// Pull in the platform/threshold configuration BEFORE the #ifndef fallbacks
// below: a TU that includes this header first would otherwise lock in the
// local defaults and silently diverge from TUs that saw a tuned profile
// (inline functions => ODR violation).
#include "../../build/PlatformConfig.h"
#include "../../build/DispatchThresholds.h"

// Default: CRT NTT enabled via size-gated dispatch (BIGMATH_NTT_CRT_THRESHOLD,
// 256 limbs sum — effectively always once NTT is selected). Below the gate,
// Goldilocks NTT runs as before — zero-cost fall-through. Opt out via
// -DBIGMATH_NTT_CRT=0.
#ifndef BIGMATH_NTT_CRT
#define BIGMATH_NTT_CRT 1
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

// Default-on Matrix Fourier Algorithm (Bailey 6-step) for very large NTT
// transforms. Recursive 2D decomposition N = N1*N2 with sub-FFTs sized to
// fit L1. Targets the bandwidth-bound regime where N exceeds L2 (~32 MB on
// M1 Max ≈ 2^23 UInt entries per prime). Threshold below = transform
// length, NOT operand limb count. Opt out via -DBIGMATH_NTT_MFA=0.
// NEON acceleration of the radix-8 butterfly layers (aarch64 only). The
// 30-bit CRT primes fit 32-bit lanes; twiddle multiplies use Shoup's method
// with a precomputed floor(w*2^31/P) companion, so the quotient estimate is
// one doubling-multiply-high (vqdmulhq_s32) per four lanes. Opt out via
// -DBIGMATH_NEON=0.
//
// The interleaved (twiddle, shoup) tables keep the strided gathers at one
// cache line per pair, so NEON wins or ties at every transform size
// (separate tables doubled twiddle traffic above L2 and lost 15-30% at
// n = 2^21 before the interleave). BIGMATH_NEON_NTT_MAX remains as an
// escape hatch; default is effectively unlimited.
#ifndef BIGMATH_NEON
#if defined(__aarch64__)
#define BIGMATH_NEON 1
#else
#define BIGMATH_NEON 0
#endif
#endif

#ifndef BIGMATH_NEON_NTT_MAX
#define BIGMATH_NEON_NTT_MAX (1 << 30)
#endif

#ifndef BIGMATH_NTT_MFA
#define BIGMATH_NTT_MFA 1
#endif

// Below this NTT length, MFA leaves run the existing radix-8 chain in-place.
// Picked so sub-FFTs of 2^13 = 8192 UInt entries (32 KB) fit M1 L1 (128 KB)
// with twiddles + workspace.
#ifndef BIGMATH_NTT_MFA_LEAF
#define BIGMATH_NTT_MFA_LEAF (1 << 13)
#endif

// MFA fires when transform length n reaches this. Canonical value lives in
// build/DispatchThresholds.h (2^20 since the 2026-06-12 post-PR-#107
// retune); this fallback only applies if that header isn't in the TU.
// History: 2^21 (too early pre-NEON) -> 2^24 (MFA lost through ~2^23 with
// whole-transform ParallelDo batching) -> 2^20: the row-chunked
// ParallelDo(6) fused stages from PR #107 keep 6-12 work units busy where
// the non-MFA path idles cores, flipping the break-even (warm-state
// n=2^22 mul 3.0x, n=2^23 2.2x vs non-MFA).
#ifndef BIGMATH_NTT_MFA_THRESHOLD
#define BIGMATH_NTT_MFA_THRESHOLD (1 << 20)
#endif

// Fuse the MFA transpose into the adjacent row-FFT pass. The plain MFA writes a
// full transposed copy to scratch and immediately reads it back for the row
// FFTs — two full round-trips through DRAM per axis. Since large MFA is memory-
// bandwidth bound on Apple silicon, that round-trip dominates. Fusing tiles the
// transpose into an L2-resident buffer, runs the row FFT there, and writes the
// result once: 8n -> 4n bytes of traffic per forward/inverse transform. The
// per-row FFT math is identical to the unfused path, so results are bit-exact.
// Only the single-level case (n1,n2 <= LEAF) is fused; that covers every
// in-practice MFA length (2^24..2^26 all factor to sub-FFTs <= 2^13). Measured
// ~11% wall-clock at 4M-10M limbs (40M-100M digits) on M1 Max. Opt out with
// -DBIGMATH_NTT_MFA_FUSE=0.
#ifndef BIGMATH_NTT_MFA_FUSE
#define BIGMATH_NTT_MFA_FUSE 1
#endif

// Rows transposed per tile block. tile = TILE * (maxSubFFT + PAD) UInt must
// stay L2-resident: 32 * 8208 * 4B = 1.0 MB. Re-swept 2026-11-21 together with
// the row padding below (was 16): at n >= 2^23 on M1 Max, 32 beats 16 by
// 10-12% and 64 by 3-5% single-threaded; with 8 threads it beats 16 by 5-6%
// and ties 64. Below 2^23 the three are within 7% of each other: 32 ahead
// single-threaded, level with 8 threads.
#ifndef BIGMATH_NTT_MFA_FUSE_TILE
#define BIGMATH_NTT_MFA_FUSE_TILE 32
#endif

// Extra UInt entries between consecutive tile rows (one 64-byte line). The
// gather writes one element into every tile row per source column, so the
// rows' write cursors advance together. With rows packed back to back they
// sit exactly one sub-FFT length apart; from a sub-FFT of 4096 entries that
// is 16 KB = one way of the 8-way 128 KB L1D on Apple P-cores, so all TILE
// cursors map to the same cache set and evict each other on every store.
// One line of padding per row walks the cursors across TILE distinct sets.
// Measured on M1 Max: single-thread per-point cost at n = 2^22 / 2^23 / 2^24
// was 73 / 114 / 131 ns unpadded and 73 / 89 / 104 ns padded (TILE = 16).
// Keep it a multiple of 16 entries (whole cache lines) so every row keeps the
// tile's alignment.
#ifndef BIGMATH_NTT_MFA_FUSE_TILE_PAD
#define BIGMATH_NTT_MFA_FUSE_TILE_PAD 16
#endif

// Output limbs from which the Garner/carry finalize pass is split across the
// thread pool. Below it one NEON pass over the residues costs less than a
// dispatch (≈ 2.5 ns per coefficient serial).
#ifndef BIGMATH_NTT_FINALIZE_PARALLEL_MIN
#define BIGMATH_NTT_FINALIZE_PARALLEL_MIN 4096
#endif

// Transform length from which the whole-transform (non-MFA) path dispatches
// its forward/inverse units to the thread pool. Below it the transforms run
// on the calling thread: at n = 1024 six threads buy about 1.3x, which is not
// worth waking (and then spinning) the pool for a ~25 us multiply.
#ifndef BIGMATH_NTT_PARALLEL_MIN
#define BIGMATH_NTT_PARALLEL_MIN 2048
#endif

// Transform length from which the three inverse transforms are split into
// six units (two half-length inverses per prime, then a shared top layer).
#ifndef BIGMATH_NTT_SPLIT_INVERSE_MIN
#define BIGMATH_NTT_SPLIT_INVERSE_MIN 4096
#endif

// Truncated transform length for products that land just above a power of
// two. A product of `c` coefficients normally needs a cyclic transform of
// length bit_ceil(c); with H = bit_ceil(c) / 2 and c only slightly above H
// that wastes almost half the transform. The H + L scheme instead computes
// the product mod (x^H + 1) (weighted, length H) and mod (x^L − 1) (length L,
// a power of two ≥ c − H) and recombines them — transform work H + L instead
// of 2H.
//
// Used when c ≤ H + (H >> BIGMATH_NTT_TRUNC_BAND_SHIFT) and
// BIGMATH_NTT_TRUNC_MIN ≤ H ≤ BIGMATH_NTT_TRUNC_MAX. Measured (M1 Max,
// threaded, vs the full 2H transform): 1.3-2.2x for c ≤ 1.25·H across
// H = 2^11..2^20, 1.0-1.35x for 1.25·H < c ≤ 1.5·H. Below H = 2048 the full
// transform is as fast (H = 1024 runs under the serial gate); above
// H = 2^20 the fused-MFA full transform wins when threaded (single-threaded,
// see BIGMATH_NTT_TRUNC_SERIAL_MAX below).
#ifndef BIGMATH_NTT_TRUNC
#define BIGMATH_NTT_TRUNC 1
#endif
#ifndef BIGMATH_NTT_TRUNC_MIN
#define BIGMATH_NTT_TRUNC_MIN 2048
#endif
#ifndef BIGMATH_NTT_TRUNC_MAX
#define BIGMATH_NTT_TRUNC_MAX (1 << 20)
#endif
#ifndef BIGMATH_NTT_TRUNC_BAND_SHIFT
#define BIGMATH_NTT_TRUNC_BAND_SHIFT 1
#endif

// Above BIGMATH_NTT_TRUNC_MAX the H + L scheme continues only when the
// process has a single thread (a -DBIGMATH_USE_THREADS=0 build, or a pool of
// one), up to this H, with the length-H and length-L transforms run through
// the fused MFA stages (CrtTruncSerial). The threaded build keeps the full 2H
// fused-MFA transform there: its tile-parallel stages scale across
// cores, the truncated pipeline has no such split, and with 8 threads the
// full transform wins. 2^25 is the largest H the 2^26 length ceiling admits.
// Set to BIGMATH_NTT_TRUNC_MAX or lower to turn the path off.
//
// Its band covers the whole octave (shift 0: c ≤ 2H, where L reaches H and
// the transform work equals the full transform's). Measured single-threaded
// on M1 Max against the full 2H fused-MFA transform, H = 2^21..2^24:
// 0.45-0.55x for c ≤ 1.25·H, 0.61-0.67x for 1.25·H < c ≤ 1.5·H, and
// 0.82-0.88x for 1.5·H < c ≤ 2H — two length-H transforms and five live
// planes cost less than one length-2H transform and nine. At H = 2^20 the
// whole-array path below BIGMATH_NTT_TRUNC_MAX is still 10-14% ahead.
#ifndef BIGMATH_NTT_TRUNC_SERIAL_MAX
#define BIGMATH_NTT_TRUNC_SERIAL_MAX (1 << 25)
#endif
#ifndef BIGMATH_NTT_TRUNC_SERIAL_BAND_SHIFT
#define BIGMATH_NTT_TRUNC_SERIAL_BAND_SHIFT 0
#endif

#ifndef NTT_MULTIPLICATION_CRT
#define NTT_MULTIPLICATION_CRT

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <vector>

#if BIGMATH_NEON
#include <arm_neon.h>
#endif

#include "../../common/DecimalConversionScope.h"
#include "../../common/Parallel.h"
#include "../../common/Util.h"
#include "ClassicMultiplication.h"
#include "../Addition.h"

namespace BigMath
{
  namespace NttCrt
  {
    // ─── Modular field over a single 30-bit NTT prime ────────────────────────

    template <UInt P>
    struct ModField
    {
      static_assert(P > 0, "prime must be non-zero");

      static constexpr UInt Prime = P;

      static inline UInt Add(UInt a, UInt b)
      {
        UInt s = a + b;
        if (s >= P) s -= P;
        return s;
      }

      static inline UInt Sub(UInt a, UInt b)
      {
        // Both inputs in [0, P); diff in (-P, P); add P if negative.
        return a >= b ? a - b : a + P - b;
      }

      static inline UInt Mul(UInt a, UInt b)
      {
        return (UInt)(((ULong)a * b) % P);
      }

      static UInt Power(UInt base, UInt exp)
      {
        UInt r = 1;
        base %= P;
        while (exp > 0)
        {
          if (exp & 1) r = Mul(r, base);
          base = Mul(base, base);
          exp >>= 1;
        }
        return r;
      }

      static UInt Inv(UInt a) { return Power(a, P - 2); }
    };

    // Three NTT-friendly primes used for the CRT scheme. Each P-1 is divisible
    // by a large power of 2; the minimum across all three is the maximum NTT
    // length the CRT path supports.
    //   P1 = 15·2^27 + 1 → supports n up to 2^27
    //   P2 =  7·2^26 + 1 → supports n up to 2^26
    //   P3 = 27·2^26 + 1 → supports n up to 2^26
    // CRT length ceiling = 2^26 = 67 108 864 coefficients ≈ 640M-digit operands.
    constexpr UInt P1 = 2013265921u;
    constexpr UInt P2 = 469762049u;
    constexpr UInt P3 = 1811939329u;
    constexpr UInt G1 = 31;
    constexpr UInt G2 = 3;
    constexpr UInt G3 = 13;

    using F1 = ModField<P1>;
    using F2 = ModField<P2>;
    using F3 = ModField<P3>;

    // Garner returns a convolution coefficient modulo p1·p2·p3 ≈ 2^90.47, so
    // the coefficient itself must stay below that. With unsigned 32-bit
    // inputs a sum of n products is at most n·(2^32 − 1)^2. With inputs
    // shifted to [−2^31, 2^31) it is at most n·2^62 in magnitude and the
    // centered residue is exact for twice as many terms. n is the coefficient
    // count of the shorter operand.
    constexpr ULong128 CrtModulus = (ULong128)P1 * P2 * P3;
    constexpr ULong CrtMaxTerms =
        (ULong)((CrtModulus - 1) / ((ULong128)0xFFFFFFFFu * 0xFFFFFFFFu)); // 92 897 280
    constexpr ULong CrtMaxBiasedTerms = (ULong)(((CrtModulus - 1) / 2) >> 62); // 185 794 560
    static_assert(CrtMaxTerms >= ((ULong)1 << 26), "one transform of 2^26 unsigned coefficients must fit");

    // ─── NTT core templated by ModField ──────────────────────────────────────

    template <typename F>
    struct Plan
    {
      std::vector<UInt> forwardRoots;
      std::vector<UInt> inverseRoots;
#if BIGMATH_NEON
      // Interleaved twiddle/Shoup pairs: il[2k] = roots[k],
      // il[2k+1] = ShoupCompanion(roots[k]) = floor(roots[k] * 2^31 / P). One cache line serves both the
      // twiddle and its Shoup companion — with separate arrays the strided
      // gathers above L2 doubled twiddle traffic and erased the NEON win at
      // n >= 2^21.
      std::vector<UInt> forwardIl;
      std::vector<UInt> inverseIl;
#endif
      UInt invSize = 1;
    };

#if BIGMATH_NEON
    // Shoup companion of a fixed multiplier w < P: floor(w * 2^31 / P), which
    // is < 2^31 and therefore a non-negative int32 — see NeonShoupMul.
    template <typename F>
    static inline UInt ShoupCompanion(UInt w)
    {
      return (UInt)(((ULong)w << 31) / F::Prime);
    }

    template <typename F>
    inline std::vector<UInt> BuildInterleaved(const std::vector<UInt> &roots)
    {
      std::vector<UInt> il(roots.size() * 2);
      const UInt *rp = roots.data();
      UInt *ip = il.data();
      Int count = (Int)roots.size();
      auto body = [rp, ip](Int start, Int end) {
        for (Int k = start; k < end; ++k)
        {
          ip[2 * k] = rp[k];
          ip[2 * k + 1] = ShoupCompanion<F>(rp[k]);
        }
      };
      if ((SizeT)count >= ParallelMinSize()) ParallelFor(count, body);
      else body(0, count);
      return il;
    }
#endif

    template <typename F, UInt G>
    inline std::vector<UInt> BuildRoots(Int n, bool invert)
    {
      UInt root = F::Power(G, (F::Prime - 1) / (UInt)n);
      if (invert) root = F::Inv(root);
      Int half = n / 2;
      std::vector<UInt> r(half);
      // Parallel chunk fill: each thread does one Power(root, chunk_start)
      // then sequential Mul within its chunk. Eliminates the 16.5M-Mul
      // serial dependency at N=33M (≈5-13% of 20M-100M mul time per row).
      UInt *rPtr = r.data();
      if ((SizeT)half >= ParallelMinSize())
      {
        ParallelFor(half, [rPtr, root](Int start, Int end) {
          if (end <= start) return;
          UInt cur = (start == 0) ? (UInt)1 : F::Power(root, (UInt)start);
          rPtr[start] = cur;
          for (Int k = start + 1; k < end; ++k)
          {
            cur = F::Mul(cur, root);
            rPtr[k] = cur;
          }
        });
      }
      else
      {
        rPtr[0] = 1;
        for (Int k = 1; k < half; ++k) rPtr[k] = F::Mul(rPtr[k - 1], root);
      }
      return r;
    }

    template <typename F, UInt G>
    inline Plan<F> BuildPlan(Int n)
    {
      Plan<F> p;
      p.forwardRoots = BuildRoots<F, G>(n, false);
      p.inverseRoots = BuildRoots<F, G>(n, true);
#if BIGMATH_NEON
      p.forwardIl = BuildInterleaved<F>(p.forwardRoots);
      p.inverseIl = BuildInterleaved<F>(p.inverseRoots);
#endif
      p.invSize = F::Inv((UInt)n);
      return p;
    }

    // Every key is a transform length (a power of two). The builder type
    // identifies the prime and kind of table. The hot lookup is thread-local.
    // Workers share immutable storage; the last owner releases it on exit.
    template <typename Table, typename Builder>
    inline const Table &CachedTransformTable(Int n, Builder build)
    {
#if BIGMATH_USE_THREADS
      struct Slot
      {
        std::mutex mutex;
        std::weak_ptr<const Table> published;
      };
      constexpr size_t Slots = std::numeric_limits<UInt>::digits;
      static std::array<Slot, Slots> slots;
      static thread_local std::array<std::shared_ptr<const Table>, Slots> owners;
      const auto index = std::countr_zero((UInt)n);
      auto &owner = owners[index];
      if (!owner)
      {
        auto &slot = slots[index];
        std::lock_guard lock(slot.mutex);
        owner = slot.published.lock();
        if (!owner)
        {
          owner = std::make_shared<const Table>(build());
          slot.published = owner;
        }
      }
      return *owner;
#else
      static thread_local std::unordered_map<Int, Table> tables;
      auto it = tables.find(n);
      if (it != tables.end()) return it->second;
      return tables.emplace(n, build()).first->second;
#endif
    }

    template <typename F, UInt G>
    inline const Plan<F> &GetPlan(Int n)
    {
      return CachedTransformTable<Plan<F>>(n, [n] { return BuildPlan<F, G>(n); });
    }

    // Always-serial Forward/Inverse. Parallelism in the CRT path is coarse-
    // grained: NttCrt::Multiply dispatches the 6 transforms (3 forward of a,
    // 3 forward of b) and 3 inverses as batched ParallelFor work units, one
    // per prime. Layer-level parallelism inside a single NTT was tried (see
    // git history) but the per-layer dispatch overhead (~30µs × log2(n)) at
    // sub-megabyte NTT sizes exceeded the parallel win. Coarse-grained
    // parallelism amortizes ~3 dispatches per Multiply across 3-6 work units.
    //
    // Butterflies are radix-4-fused: each layer-pair issues 4 muls + 8 add/subs
    // per group of 4 elements, with one load+store pair instead of two —
    // halves memory traffic vs separate radix-2 layers. Mirrors NTTCore.
    template <typename F>
    inline void ForwardRadix2Layer(UInt *a, Int n, Int len, const UInt *roots)
    {
      Int halflen = len >> 1;
      Int stride = n / len;
      for (Int i = 0; i < n; i += len)
      {
        for (Int j = 0; j < halflen; ++j)
        {
          UInt u = a[i + j];
          UInt v = a[i + j + halflen];
          a[i + j] = F::Add(u, v);
          a[i + j + halflen] = F::Mul(F::Sub(u, v), roots[j * stride]);
        }
      }
    }

    template <typename F>
    inline void ForwardRadix4Layer(UInt *a, Int n, Int outerLen, const UInt *roots)
    {
      Int halflen = outerLen >> 1;
      Int qlen = outerLen >> 2;
      Int stride = n / outerLen;
      Int omega4_off = n / 4;
      for (Int i = 0; i < n; i += outerLen)
      {
        for (Int j = 0; j < qlen; ++j)
        {
          UInt x0 = a[i + j];
          UInt x1 = a[i + j + qlen];
          UInt x2 = a[i + j + halflen];
          UInt x3 = a[i + j + halflen + qlen];

          UInt g  = roots[j * stride];
          UInt gw = roots[j * stride + omega4_off];
          UInt g2 = roots[2 * j * stride];

          UInt t0 = F::Add(x0, x2);
          UInt t1 = F::Add(x1, x3);
          UInt t2 = F::Mul(F::Sub(x0, x2), g);
          UInt t3 = F::Mul(F::Sub(x1, x3), gw);

          a[i + j]                  = F::Add(t0, t1);
          a[i + j + qlen]           = F::Mul(F::Sub(t0, t1), g2);
          a[i + j + halflen]        = F::Add(t2, t3);
          a[i + j + halflen + qlen] = F::Mul(F::Sub(t2, t3), g2);
        }
      }
    }

    template <typename F>
    inline void InverseRadix2Layer(UInt *a, Int n, Int len, const UInt *roots)
    {
      Int halflen = len >> 1;
      Int stride = n / len;
      for (Int i = 0; i < n; i += len)
      {
        for (Int j = 0; j < halflen; ++j)
        {
          UInt u = a[i + j];
          UInt v = F::Mul(a[i + j + halflen], roots[j * stride]);
          a[i + j] = F::Add(u, v);
          a[i + j + halflen] = F::Sub(u, v);
        }
      }
    }

    template <typename F>
    inline void InverseRadix4Layer(UInt *a, Int n, Int outerLen, const UInt *roots)
    {
      Int halflen = outerLen >> 1;
      Int qlen = outerLen >> 2;
      Int stride = n / outerLen;
      Int omega4_off = n / 4;
      for (Int i = 0; i < n; i += outerLen)
      {
        for (Int j = 0; j < qlen; ++j)
        {
          UInt x0 = a[i + j];
          UInt x1 = a[i + j + qlen];
          UInt x2 = a[i + j + halflen];
          UInt x3 = a[i + j + halflen + qlen];

          UInt g  = roots[j * stride];
          UInt gw = roots[j * stride + omega4_off];
          UInt g2 = roots[2 * j * stride];

          UInt y1 = F::Mul(x1, g2);
          UInt y3 = F::Mul(x3, g2);
          UInt t0 = F::Add(x0, y1);
          UInt t1 = F::Sub(x0, y1);
          UInt t2 = F::Mul(F::Add(x2, y3), g);
          UInt t3 = F::Mul(F::Sub(x2, y3), gw);

          a[i + j]                  = F::Add(t0, t2);
          a[i + j + halflen]        = F::Sub(t0, t2);
          a[i + j + qlen]           = F::Add(t1, t3);
          a[i + j + halflen + qlen] = F::Sub(t1, t3);
        }
      }
    }

    // Radix-8 fused DIF: three layers (len, len/2, len/4) in one load/store.
#if BIGMATH_NEON
    // ─── NEON butterfly primitives ──────────────────────────────────────────
    // All values in [0, P), P < 2^31.
    //
    // ShoupMul computes (w*x) mod P from the precomputed companion
    // ws = floor(w * 2^31 / P). Both ws and x are < 2^31, i.e. non-negative
    // as int32, so the doubling multiply-high vqdmulhq_s32(ws, x) =
    // floor(2*ws*x / 2^32) = floor(ws*x / 2^31) cannot saturate and yields
    // the quotient estimate for all four lanes in ONE instruction (the
    // unsigned 2^32 form needed two widening multiplies, two narrowing
    // shifts and a combine). Since ws*x/2^31 lies in (w*x/P - 1, w*x/P], the
    // estimate q is floor(w*x/P) or one less, so t = w*x - q*P (mod 2^32) is
    // the true value in [0, 2P) — the same bound as classic Shoup — and one
    // conditional subtract finishes it. Bit-exact with (w*x) % P.
    //
    // Conditional subtract / add use an unsigned min instead of
    // compare + mask: for t in [0, 2P), t - P wraps above t exactly when
    // t < P, so min(t, t - P) selects the reduced value in two ops.
    template <typename F>
    static inline uint32x4_t NeonShoupMul(uint32x4_t w, uint32x4_t ws, uint32x4_t x)
    {
      const uint32x4_t vP = vdupq_n_u32(F::Prime);
      uint32x4_t q = vreinterpretq_u32_s32(
          vqdmulhq_s32(vreinterpretq_s32_u32(ws), vreinterpretq_s32_u32(x)));
      uint32x4_t t = vmlsq_u32(vmulq_u32(w, x), q, vP);
      return vminq_u32(t, vsubq_u32(t, vP));
    }

    template <typename F>
    static inline uint32x4_t NeonAdd(uint32x4_t a, uint32x4_t b)
    {
      const uint32x4_t vP = vdupq_n_u32(F::Prime);
      uint32x4_t sum = vaddq_u32(a, b); // < 2P < 2^32
      return vminq_u32(sum, vsubq_u32(sum, vP));
    }

    template <typename F>
    static inline uint32x4_t NeonSub(uint32x4_t a, uint32x4_t b)
    {
      const uint32x4_t vP = vdupq_n_u32(F::Prime);
      // a >= b: d in [0, P) and d + P > d. a < b: d wrapped above 2^32 - P
      // and d + P wraps back to the reduced value below it.
      uint32x4_t d = vsubq_u32(a, b);
      return vminq_u32(d, vaddq_u32(d, vP));
    }

    // Strided gather of 4 consecutive-j (twiddle, shoup) pairs from the
    // interleaved table: one 64-bit load per lane brings both halves on the
    // same cache line; two unzips split them into the (w, w') vectors.
    struct NeonPair
    {
      uint32x4_t w;
      uint32x4_t ws;
    };

    // Loads use vld1_u32 on the UInt elements directly — never a
    // reinterpret_cast<const uint64_t*> over the uint32_t table, which would
    // be a strict-aliasing violation (works under today's AppleClang, but
    // this header must not depend on a compiler shrugging at UB; see the
    // PR #103 post-mortem for what silent miscompilation costs here).
    static inline NeonPair NeonGatherPair(const UInt *il, Int base, Int step)
    {
      uint32x2_t d0 = vld1_u32(il + 2 * base);
      uint32x2_t d1 = vld1_u32(il + 2 * (base + step));
      uint32x2_t d2 = vld1_u32(il + 2 * (base + 2 * step));
      uint32x2_t d3 = vld1_u32(il + 2 * (base + 3 * step));
      uint32x4_t a = vcombine_u32(d0, d1);
      uint32x4_t b = vcombine_u32(d2, d3);
      NeonPair r;
      r.w = vuzp1q_u32(a, b);
      r.ws = vuzp2q_u32(a, b);
      return r;
    }

    // ─── NEON tail layers (outerLen == 4 / 2, the log n % 3 remainder) ──────
    // At outerLen == 4 the twiddles collapse to constants: g = g2 = roots[0]
    // = 1 and gw = roots[n/4] — one broadcast Shoup pair serves the whole
    // pass, and the groups-of-4 are contiguous, so vld4q/vst4q de-interleave
    // does the lane work.
    template <typename F>
    inline void ForwardRadix4LayerNeon(UInt *a, Int n, const UInt *roots)
    {
      const UInt gw = roots[n / 4];
      const uint32x4_t vgw = vdupq_n_u32(gw);
      const uint32x4_t vgws = vdupq_n_u32(ShoupCompanion<F>(gw));
      Int i = 0;
      for (; i + 16 <= n; i += 16)
      {
        uint32x4x4_t x = vld4q_u32(a + i);
        uint32x4_t t0 = NeonAdd<F>(x.val[0], x.val[2]);
        uint32x4_t t1 = NeonAdd<F>(x.val[1], x.val[3]);
        uint32x4_t t2 = NeonSub<F>(x.val[0], x.val[2]);
        uint32x4_t t3 = NeonShoupMul<F>(vgw, vgws, NeonSub<F>(x.val[1], x.val[3]));
        uint32x4x4_t y;
        y.val[0] = NeonAdd<F>(t0, t1);
        y.val[1] = NeonSub<F>(t0, t1);
        y.val[2] = NeonAdd<F>(t2, t3);
        y.val[3] = NeonSub<F>(t2, t3);
        vst4q_u32(a + i, y);
      }
      for (; i < n; i += 4)
      {
        UInt x0 = a[i], x1 = a[i + 1], x2 = a[i + 2], x3 = a[i + 3];
        UInt t0 = F::Add(x0, x2);
        UInt t1 = F::Add(x1, x3);
        UInt t2 = F::Sub(x0, x2);
        UInt t3 = F::Mul(F::Sub(x1, x3), gw);
        a[i] = F::Add(t0, t1);
        a[i + 1] = F::Sub(t0, t1);
        a[i + 2] = F::Add(t2, t3);
        a[i + 3] = F::Sub(t2, t3);
      }
    }

    template <typename F>
    inline void InverseRadix4LayerNeon(UInt *a, Int n, const UInt *roots)
    {
      const UInt gw = roots[n / 4];
      const uint32x4_t vgw = vdupq_n_u32(gw);
      const uint32x4_t vgws = vdupq_n_u32(ShoupCompanion<F>(gw));
      Int i = 0;
      for (; i + 16 <= n; i += 16)
      {
        uint32x4x4_t x = vld4q_u32(a + i);
        uint32x4_t t0 = NeonAdd<F>(x.val[0], x.val[1]);
        uint32x4_t t1 = NeonSub<F>(x.val[0], x.val[1]);
        uint32x4_t t2 = NeonAdd<F>(x.val[2], x.val[3]);
        uint32x4_t t3 = NeonShoupMul<F>(vgw, vgws, NeonSub<F>(x.val[2], x.val[3]));
        uint32x4x4_t y;
        y.val[0] = NeonAdd<F>(t0, t2);
        y.val[2] = NeonSub<F>(t0, t2);
        y.val[1] = NeonAdd<F>(t1, t3);
        y.val[3] = NeonSub<F>(t1, t3);
        vst4q_u32(a + i, y);
      }
      for (; i < n; i += 4)
      {
        UInt x0 = a[i], x1 = a[i + 1], x2 = a[i + 2], x3 = a[i + 3];
        UInt t0 = F::Add(x0, x1);
        UInt t1 = F::Sub(x0, x1);
        UInt t2 = F::Add(x2, x3);
        UInt t3 = F::Mul(F::Sub(x2, x3), gw);
        a[i] = F::Add(t0, t2);
        a[i + 2] = F::Sub(t0, t2);
        a[i + 1] = F::Add(t1, t3);
        a[i + 3] = F::Sub(t1, t3);
      }
    }

    template <typename F>
    inline void InverseRadix2LayerNeon(UInt *a, Int n)
    {
      // len == 2: roots[0] == 1, pure add/sub on adjacent pairs.
      Int i = 0;
      for (; i + 8 <= n; i += 8)
      {
        uint32x4x2_t x = vld2q_u32(a + i);
        uint32x4x2_t y;
        y.val[0] = NeonAdd<F>(x.val[0], x.val[1]);
        y.val[1] = NeonSub<F>(x.val[0], x.val[1]);
        vst2q_u32(a + i, y);
      }
      for (; i < n; i += 2)
      {
        UInt u = a[i], v = a[i + 1];
        a[i] = F::Add(u, v);
        a[i + 1] = F::Sub(u, v);
      }
    }

    // ─── Fused NEON tails for the log n % 3 != 2 size classes ───────────────
    // The radix-8 chain leaves len == 16 (log n % 3 == 1) or len == 8
    // (log n % 3 == 0) for its last step. At those lengths the radix-8
    // j-loop has only 2 or 1 iterations per block, so its 4-lane kernel
    // cannot run and the layer used to fall back to scalar (a*b) % P — three
    // to four of every transform's layers, on two thirds of all sizes
    // (measured 0.60-0.63 ns per butterfly vs 0.43 for the all-NEON class).
    //
    // These tails cover the same layers with 4-lane steps instead. One step
    // pairs the 4-element quarters of each 16-element block (its twiddles
    // depend only on j = 0..3, so they are loop-invariant vectors); a 4x4
    // transpose then turns "quarter" into "lane" for the len-4 radix-4 step,
    // exactly as ForwardRadix4LayerNeon does through vld4q/vst4q. All
    // arithmetic is exact mod P, so the output is bit-identical to the
    // radix-2 chain regardless of how the layers are grouped.

    // In-place 4x4 transpose of four rows (an involution).
    static inline void NeonTranspose4(uint32x4_t &r0, uint32x4_t &r1, uint32x4_t &r2, uint32x4_t &r3)
    {
      uint64x2_t t0 = vreinterpretq_u64_u32(vtrn1q_u32(r0, r1));
      uint64x2_t t1 = vreinterpretq_u64_u32(vtrn2q_u32(r0, r1));
      uint64x2_t t2 = vreinterpretq_u64_u32(vtrn1q_u32(r2, r3));
      uint64x2_t t3 = vreinterpretq_u64_u32(vtrn2q_u32(r2, r3));
      r0 = vreinterpretq_u32_u64(vtrn1q_u64(t0, t2));
      r1 = vreinterpretq_u32_u64(vtrn1q_u64(t1, t3));
      r2 = vreinterpretq_u32_u64(vtrn2q_u64(t0, t2));
      r3 = vreinterpretq_u32_u64(vtrn2q_u64(t1, t3));
    }

    // Loop-invariant (twiddle, Shoup companion) vector pair for j = 0..3.
    template <typename F>
    static inline NeonPair NeonConstPair(const UInt *roots, Int base, Int step)
    {
      UInt w[4], ws[4];
      for (Int j = 0; j < 4; ++j)
      {
        w[j] = roots[base + j * step];
        ws[j] = ShoupCompanion<F>(w[j]);
      }
      NeonPair r;
      r.w = vld1q_u32(w);
      r.ws = vld1q_u32(ws);
      return r;
    }

    // Last four DIF layers (len 16, 8, 4, 2) in one pass. Requires n >= 16.
    template <typename F>
    inline void ForwardTail16Neon(UInt *a, Int n, const UInt *roots)
    {
      const Int stride = n / 16;
      const NeonPair g = NeonConstPair<F>(roots, 0, stride);
      const NeonPair gw = NeonConstPair<F>(roots, n / 4, stride);
      const NeonPair g2 = NeonConstPair<F>(roots, 0, 2 * stride);
      const UInt gwT = roots[n / 4];
      const uint32x4_t vgwT = vdupq_n_u32(gwT);
      const uint32x4_t vgwTs = vdupq_n_u32(ShoupCompanion<F>(gwT));
      for (Int i = 0; i < n; i += 16)
      {
        uint32x4_t x0 = vld1q_u32(a + i);
        uint32x4_t x1 = vld1q_u32(a + i + 4);
        uint32x4_t x2 = vld1q_u32(a + i + 8);
        uint32x4_t x3 = vld1q_u32(a + i + 12);

        // Radix-4 step at outerLen == 16 (lanes are j).
        uint32x4_t t0 = NeonAdd<F>(x0, x2);
        uint32x4_t t1 = NeonAdd<F>(x1, x3);
        uint32x4_t t2 = NeonShoupMul<F>(g.w, g.ws, NeonSub<F>(x0, x2));
        uint32x4_t t3 = NeonShoupMul<F>(gw.w, gw.ws, NeonSub<F>(x1, x3));
        uint32x4_t y0 = NeonAdd<F>(t0, t1);
        uint32x4_t y1 = NeonShoupMul<F>(g2.w, g2.ws, NeonSub<F>(t0, t1));
        uint32x4_t y2 = NeonAdd<F>(t2, t3);
        uint32x4_t y3 = NeonShoupMul<F>(g2.w, g2.ws, NeonSub<F>(t2, t3));

        // Each y_k is one contiguous len-4 group; transpose so vector index
        // is the element within its group, then the len-4 radix-4 step.
        NeonTranspose4(y0, y1, y2, y3);
        uint32x4_t u0 = NeonAdd<F>(y0, y2);
        uint32x4_t u1 = NeonAdd<F>(y1, y3);
        uint32x4_t u2 = NeonSub<F>(y0, y2);
        uint32x4_t u3 = NeonShoupMul<F>(vgwT, vgwTs, NeonSub<F>(y1, y3));
        uint32x4x4_t z;
        z.val[0] = NeonAdd<F>(u0, u1);
        z.val[1] = NeonSub<F>(u0, u1);
        z.val[2] = NeonAdd<F>(u2, u3);
        z.val[3] = NeonSub<F>(u2, u3);
        vst4q_u32(a + i, z); // interleaving store == transpose back
      }
    }

    // Last three DIF layers (len 8, 4, 2) in one pass. Requires n >= 16.
    template <typename F>
    inline void ForwardTail8Neon(UInt *a, Int n, const UInt *roots)
    {
      const NeonPair w = NeonConstPair<F>(roots, 0, n / 8);
      const UInt gwT = roots[n / 4];
      const uint32x4_t vgwT = vdupq_n_u32(gwT);
      const uint32x4_t vgwTs = vdupq_n_u32(ShoupCompanion<F>(gwT));
      for (Int i = 0; i < n; i += 16)
      {
        uint32x4_t lo0 = vld1q_u32(a + i);
        uint32x4_t hi0 = vld1q_u32(a + i + 4);
        uint32x4_t lo1 = vld1q_u32(a + i + 8);
        uint32x4_t hi1 = vld1q_u32(a + i + 12);

        // Radix-2 step at len == 8 on two adjacent blocks (lanes are j).
        uint32x4_t y0 = NeonAdd<F>(lo0, hi0);
        uint32x4_t y1 = NeonShoupMul<F>(w.w, w.ws, NeonSub<F>(lo0, hi0));
        uint32x4_t y2 = NeonAdd<F>(lo1, hi1);
        uint32x4_t y3 = NeonShoupMul<F>(w.w, w.ws, NeonSub<F>(lo1, hi1));

        NeonTranspose4(y0, y1, y2, y3);
        uint32x4_t u0 = NeonAdd<F>(y0, y2);
        uint32x4_t u1 = NeonAdd<F>(y1, y3);
        uint32x4_t u2 = NeonSub<F>(y0, y2);
        uint32x4_t u3 = NeonShoupMul<F>(vgwT, vgwTs, NeonSub<F>(y1, y3));
        uint32x4x4_t z;
        z.val[0] = NeonAdd<F>(u0, u1);
        z.val[1] = NeonSub<F>(u0, u1);
        z.val[2] = NeonAdd<F>(u2, u3);
        z.val[3] = NeonSub<F>(u2, u3);
        vst4q_u32(a + i, z);
      }
    }

    // First four DIT layers (len 2, 4, 8, 16) in one pass. Requires n >= 16.
    template <typename F>
    inline void InverseTail16Neon(UInt *a, Int n, const UInt *roots)
    {
      const Int stride = n / 16;
      const NeonPair g = NeonConstPair<F>(roots, 0, stride);
      const NeonPair gw = NeonConstPair<F>(roots, n / 4, stride);
      const NeonPair g2 = NeonConstPair<F>(roots, 0, 2 * stride);
      const UInt gwT = roots[n / 4];
      const uint32x4_t vgwT = vdupq_n_u32(gwT);
      const uint32x4_t vgwTs = vdupq_n_u32(ShoupCompanion<F>(gwT));
      for (Int i = 0; i < n; i += 16)
      {
        // De-interleaving load: x.val[e] holds element e of each len-4 group.
        uint32x4x4_t x = vld4q_u32(a + i);
        uint32x4_t t0 = NeonAdd<F>(x.val[0], x.val[1]);
        uint32x4_t t1 = NeonSub<F>(x.val[0], x.val[1]);
        uint32x4_t t2 = NeonAdd<F>(x.val[2], x.val[3]);
        uint32x4_t t3 = NeonShoupMul<F>(vgwT, vgwTs, NeonSub<F>(x.val[2], x.val[3]));
        uint32x4_t r0 = NeonAdd<F>(t0, t2);
        uint32x4_t r2 = NeonSub<F>(t0, t2);
        uint32x4_t r1 = NeonAdd<F>(t1, t3);
        uint32x4_t r3 = NeonSub<F>(t1, t3);

        // Back to one vector per group, then the radix-4 step at outerLen == 16.
        NeonTranspose4(r0, r1, r2, r3);
        uint32x4_t p1 = NeonShoupMul<F>(g2.w, g2.ws, r1);
        uint32x4_t p3 = NeonShoupMul<F>(g2.w, g2.ws, r3);
        uint32x4_t q0 = NeonAdd<F>(r0, p1);
        uint32x4_t q1 = NeonSub<F>(r0, p1);
        uint32x4_t q2 = NeonShoupMul<F>(g.w, g.ws, NeonAdd<F>(r2, p3));
        uint32x4_t q3 = NeonShoupMul<F>(gw.w, gw.ws, NeonSub<F>(r2, p3));
        vst1q_u32(a + i, NeonAdd<F>(q0, q2));
        vst1q_u32(a + i + 8, NeonSub<F>(q0, q2));
        vst1q_u32(a + i + 4, NeonAdd<F>(q1, q3));
        vst1q_u32(a + i + 12, NeonSub<F>(q1, q3));
      }
    }

    // First three DIT layers (len 2, 4, 8) in one pass. Requires n >= 16.
    template <typename F>
    inline void InverseTail8Neon(UInt *a, Int n, const UInt *roots)
    {
      const NeonPair w = NeonConstPair<F>(roots, 0, n / 8);
      const UInt gwT = roots[n / 4];
      const uint32x4_t vgwT = vdupq_n_u32(gwT);
      const uint32x4_t vgwTs = vdupq_n_u32(ShoupCompanion<F>(gwT));
      for (Int i = 0; i < n; i += 16)
      {
        uint32x4x4_t x = vld4q_u32(a + i);
        uint32x4_t t0 = NeonAdd<F>(x.val[0], x.val[1]);
        uint32x4_t t1 = NeonSub<F>(x.val[0], x.val[1]);
        uint32x4_t t2 = NeonAdd<F>(x.val[2], x.val[3]);
        uint32x4_t t3 = NeonShoupMul<F>(vgwT, vgwTs, NeonSub<F>(x.val[2], x.val[3]));
        uint32x4_t r0 = NeonAdd<F>(t0, t2);
        uint32x4_t r2 = NeonSub<F>(t0, t2);
        uint32x4_t r1 = NeonAdd<F>(t1, t3);
        uint32x4_t r3 = NeonSub<F>(t1, t3);

        // r0/r1 and r2/r3 become the low/high halves of two len-8 blocks.
        NeonTranspose4(r0, r1, r2, r3);
        uint32x4_t v0 = NeonShoupMul<F>(w.w, w.ws, r1);
        uint32x4_t v1 = NeonShoupMul<F>(w.w, w.ws, r3);
        vst1q_u32(a + i, NeonAdd<F>(r0, v0));
        vst1q_u32(a + i + 4, NeonSub<F>(r0, v0));
        vst1q_u32(a + i + 8, NeonAdd<F>(r2, v1));
        vst1q_u32(a + i + 12, NeonSub<F>(r2, v1));
      }
    }

    // Constant-multiplier Shoup pass for the inverse 1/n scaling.
    template <typename F>
    inline void NeonScale(UInt *a, Int n, UInt c)
    {
      const uint32x4_t vc = vdupq_n_u32(c);
      const uint32x4_t vcs = vdupq_n_u32(ShoupCompanion<F>(c));
      Int i = 0;
      for (; i + 4 <= n; i += 4)
        vst1q_u32(a + i, NeonShoupMul<F>(vc, vcs, vld1q_u32(a + i)));
      for (; i < n; ++i)
        a[i] = F::Mul(a[i], c);
    }

    // NEON 4-lane j-loop of the radix-8 DIF forward layer. Bit-exact with the
    // scalar layer below (Shoup and (a*b)%P agree exactly for inputs < P).
    template <typename F>
    inline void ForwardRadix8LayerNeon(UInt *a, Int n, Int outerLen,
                                       const UInt *roots, const UInt *il)
    {
      Int halflen = outerLen >> 1;
      Int qlen4 = outerLen >> 2;
      Int qlen8 = outerLen >> 3;
      Int stride = n / outerLen;
      Int stride4 = stride << 2;
      Int omega4_off = n / 4;
      for (Int i = 0; i < n; i += outerLen)
      {
        Int j = 0;
        for (; j + 4 <= qlen8; j += 4)
        {
          uint32x4_t x0 = vld1q_u32(a + i + j);
          uint32x4_t x1 = vld1q_u32(a + i + j + qlen8);
          uint32x4_t x2 = vld1q_u32(a + i + j + qlen4);
          uint32x4_t x3 = vld1q_u32(a + i + j + qlen4 + qlen8);
          uint32x4_t x4 = vld1q_u32(a + i + j + halflen);
          uint32x4_t x5 = vld1q_u32(a + i + j + halflen + qlen8);
          uint32x4_t x6 = vld1q_u32(a + i + j + halflen + qlen4);
          uint32x4_t x7 = vld1q_u32(a + i + j + halflen + qlen4 + qlen8);

          NeonPair g_a  = NeonGatherPair(il, j * stride, stride);
          NeonPair gw_a = NeonGatherPair(il, j * stride + omega4_off, stride);
          NeonPair g2_a = NeonGatherPair(il, 2 * j * stride, 2 * stride);
          NeonPair g_b  = NeonGatherPair(il, (j + qlen8) * stride, stride);
          NeonPair gw_b = NeonGatherPair(il, (j + qlen8) * stride + omega4_off, stride);
          NeonPair g2_b = NeonGatherPair(il, 2 * (j + qlen8) * stride, 2 * stride);
          NeonPair g3   = NeonGatherPair(il, j * stride4, stride4);

          uint32x4_t t0 = NeonAdd<F>(x0, x4);
          uint32x4_t t1 = NeonAdd<F>(x2, x6);
          uint32x4_t t2 = NeonShoupMul<F>(g_a.w, g_a.ws, NeonSub<F>(x0, x4));
          uint32x4_t t3 = NeonShoupMul<F>(gw_a.w, gw_a.ws, NeonSub<F>(x2, x6));
          uint32x4_t y0 = NeonAdd<F>(t0, t1);
          uint32x4_t y2 = NeonShoupMul<F>(g2_a.w, g2_a.ws, NeonSub<F>(t0, t1));
          uint32x4_t y4 = NeonAdd<F>(t2, t3);
          uint32x4_t y6 = NeonShoupMul<F>(g2_a.w, g2_a.ws, NeonSub<F>(t2, t3));

          uint32x4_t u0 = NeonAdd<F>(x1, x5);
          uint32x4_t u1 = NeonAdd<F>(x3, x7);
          uint32x4_t u2 = NeonShoupMul<F>(g_b.w, g_b.ws, NeonSub<F>(x1, x5));
          uint32x4_t u3 = NeonShoupMul<F>(gw_b.w, gw_b.ws, NeonSub<F>(x3, x7));
          uint32x4_t y1 = NeonAdd<F>(u0, u1);
          uint32x4_t y3 = NeonShoupMul<F>(g2_b.w, g2_b.ws, NeonSub<F>(u0, u1));
          uint32x4_t y5 = NeonAdd<F>(u2, u3);
          uint32x4_t y7 = NeonShoupMul<F>(g2_b.w, g2_b.ws, NeonSub<F>(u2, u3));

          uint32x4_t z0 = NeonAdd<F>(y0, y1);
          uint32x4_t z1 = NeonShoupMul<F>(g3.w, g3.ws, NeonSub<F>(y0, y1));
          uint32x4_t z2 = NeonAdd<F>(y2, y3);
          uint32x4_t z3 = NeonShoupMul<F>(g3.w, g3.ws, NeonSub<F>(y2, y3));
          uint32x4_t z4 = NeonAdd<F>(y4, y5);
          uint32x4_t z5 = NeonShoupMul<F>(g3.w, g3.ws, NeonSub<F>(y4, y5));
          uint32x4_t z6 = NeonAdd<F>(y6, y7);
          uint32x4_t z7 = NeonShoupMul<F>(g3.w, g3.ws, NeonSub<F>(y6, y7));

          vst1q_u32(a + i + j, z0);
          vst1q_u32(a + i + j + qlen8, z1);
          vst1q_u32(a + i + j + qlen4, z2);
          vst1q_u32(a + i + j + qlen4 + qlen8, z3);
          vst1q_u32(a + i + j + halflen, z4);
          vst1q_u32(a + i + j + halflen + qlen8, z5);
          vst1q_u32(a + i + j + halflen + qlen4, z6);
          vst1q_u32(a + i + j + halflen + qlen4 + qlen8, z7);
        }
        // Scalar tail (qlen8 % 4 remainder).
        for (; j < qlen8; ++j)
        {
          UInt x0 = a[i + j];
          UInt x1 = a[i + j + qlen8];
          UInt x2 = a[i + j + qlen4];
          UInt x3 = a[i + j + qlen4 + qlen8];
          UInt x4 = a[i + j + halflen];
          UInt x5 = a[i + j + halflen + qlen8];
          UInt x6 = a[i + j + halflen + qlen4];
          UInt x7 = a[i + j + halflen + qlen4 + qlen8];

          UInt g_a  = roots[j * stride];
          UInt gw_a = roots[j * stride + omega4_off];
          UInt g2_a = roots[2 * j * stride];
          UInt g_b  = roots[(j + qlen8) * stride];
          UInt gw_b = roots[(j + qlen8) * stride + omega4_off];
          UInt g2_b = roots[2 * (j + qlen8) * stride];
          UInt g3   = roots[j * stride4];

          UInt t0 = F::Add(x0, x4);
          UInt t1 = F::Add(x2, x6);
          UInt t2 = F::Mul(F::Sub(x0, x4), g_a);
          UInt t3 = F::Mul(F::Sub(x2, x6), gw_a);
          UInt y0 = F::Add(t0, t1);
          UInt y2 = F::Mul(F::Sub(t0, t1), g2_a);
          UInt y4 = F::Add(t2, t3);
          UInt y6 = F::Mul(F::Sub(t2, t3), g2_a);

          UInt u0 = F::Add(x1, x5);
          UInt u1 = F::Add(x3, x7);
          UInt u2 = F::Mul(F::Sub(x1, x5), g_b);
          UInt u3 = F::Mul(F::Sub(x3, x7), gw_b);
          UInt y1 = F::Add(u0, u1);
          UInt y3 = F::Mul(F::Sub(u0, u1), g2_b);
          UInt y5 = F::Add(u2, u3);
          UInt y7 = F::Mul(F::Sub(u2, u3), g2_b);

          UInt z0 = F::Add(y0, y1);
          UInt z1 = F::Mul(F::Sub(y0, y1), g3);
          UInt z2 = F::Add(y2, y3);
          UInt z3 = F::Mul(F::Sub(y2, y3), g3);
          UInt z4 = F::Add(y4, y5);
          UInt z5 = F::Mul(F::Sub(y4, y5), g3);
          UInt z6 = F::Add(y6, y7);
          UInt z7 = F::Mul(F::Sub(y6, y7), g3);

          a[i + j]                          = z0;
          a[i + j + qlen8]                  = z1;
          a[i + j + qlen4]                  = z2;
          a[i + j + qlen4 + qlen8]          = z3;
          a[i + j + halflen]                = z4;
          a[i + j + halflen + qlen8]        = z5;
          a[i + j + halflen + qlen4]        = z6;
          a[i + j + halflen + qlen4 + qlen8] = z7;
        }
      }
    }

    // NEON 4-lane j-loop of the radix-8 DIT inverse layer.
    template <typename F>
    inline void InverseRadix8LayerNeon(UInt *a, Int n, Int outerLen,
                                       const UInt *roots, const UInt *il)
    {
      Int halflen = outerLen >> 1;
      Int qlen4 = outerLen >> 2;
      Int qlen8 = outerLen >> 3;
      Int stride = n / outerLen;
      Int stride4 = stride << 2;
      for (Int i = 0; i < n; i += outerLen)
      {
        Int j = 0;
        for (; j + 4 <= qlen8; j += 4)
        {
          uint32x4_t x0 = vld1q_u32(a + i + j);
          uint32x4_t x1 = vld1q_u32(a + i + j + qlen8);
          uint32x4_t x2 = vld1q_u32(a + i + j + qlen4);
          uint32x4_t x3 = vld1q_u32(a + i + j + qlen4 + qlen8);
          uint32x4_t x4 = vld1q_u32(a + i + j + halflen);
          uint32x4_t x5 = vld1q_u32(a + i + j + halflen + qlen8);
          uint32x4_t x6 = vld1q_u32(a + i + j + halflen + qlen4);
          uint32x4_t x7 = vld1q_u32(a + i + j + halflen + qlen4 + qlen8);

          NeonPair g_a = NeonGatherPair(il, j * stride4, stride4);

          uint32x4_t v01 = NeonShoupMul<F>(g_a.w, g_a.ws, x1);
          uint32x4_t y0 = NeonAdd<F>(x0, v01);
          uint32x4_t y1 = NeonSub<F>(x0, v01);
          uint32x4_t v23 = NeonShoupMul<F>(g_a.w, g_a.ws, x3);
          uint32x4_t y2 = NeonAdd<F>(x2, v23);
          uint32x4_t y3 = NeonSub<F>(x2, v23);
          uint32x4_t v45 = NeonShoupMul<F>(g_a.w, g_a.ws, x5);
          uint32x4_t y4 = NeonAdd<F>(x4, v45);
          uint32x4_t y5 = NeonSub<F>(x4, v45);
          uint32x4_t v67 = NeonShoupMul<F>(g_a.w, g_a.ws, x7);
          uint32x4_t y6 = NeonAdd<F>(x6, v67);
          uint32x4_t y7 = NeonSub<F>(x6, v67);

          NeonPair g_b1 = NeonGatherPair(il, 2 * j * stride, 2 * stride);
          NeonPair g_b2 = NeonGatherPair(il, 2 * (j + qlen8) * stride, 2 * stride);

          uint32x4_t w02 = NeonShoupMul<F>(g_b1.w, g_b1.ws, y2);
          uint32x4_t z0 = NeonAdd<F>(y0, w02);
          uint32x4_t z2 = NeonSub<F>(y0, w02);
          uint32x4_t w13 = NeonShoupMul<F>(g_b2.w, g_b2.ws, y3);
          uint32x4_t z1 = NeonAdd<F>(y1, w13);
          uint32x4_t z3 = NeonSub<F>(y1, w13);
          uint32x4_t w46 = NeonShoupMul<F>(g_b1.w, g_b1.ws, y6);
          uint32x4_t z4 = NeonAdd<F>(y4, w46);
          uint32x4_t z6 = NeonSub<F>(y4, w46);
          uint32x4_t w57 = NeonShoupMul<F>(g_b2.w, g_b2.ws, y7);
          uint32x4_t z5 = NeonAdd<F>(y5, w57);
          uint32x4_t z7 = NeonSub<F>(y5, w57);

          NeonPair g_c0 = NeonGatherPair(il, j * stride, stride);
          NeonPair g_c1 = NeonGatherPair(il, (j + qlen8) * stride, stride);
          NeonPair g_c2 = NeonGatherPair(il, (j + qlen4) * stride, stride);
          NeonPair g_c3 = NeonGatherPair(il, (j + qlen4 + qlen8) * stride, stride);

          uint32x4_t w04 = NeonShoupMul<F>(g_c0.w, g_c0.ws, z4);
          uint32x4_t w15 = NeonShoupMul<F>(g_c1.w, g_c1.ws, z5);
          uint32x4_t w26 = NeonShoupMul<F>(g_c2.w, g_c2.ws, z6);
          uint32x4_t w37 = NeonShoupMul<F>(g_c3.w, g_c3.ws, z7);

          vst1q_u32(a + i + j, NeonAdd<F>(z0, w04));
          vst1q_u32(a + i + j + halflen, NeonSub<F>(z0, w04));
          vst1q_u32(a + i + j + qlen8, NeonAdd<F>(z1, w15));
          vst1q_u32(a + i + j + halflen + qlen8, NeonSub<F>(z1, w15));
          vst1q_u32(a + i + j + qlen4, NeonAdd<F>(z2, w26));
          vst1q_u32(a + i + j + halflen + qlen4, NeonSub<F>(z2, w26));
          vst1q_u32(a + i + j + qlen4 + qlen8, NeonAdd<F>(z3, w37));
          vst1q_u32(a + i + j + halflen + qlen4 + qlen8, NeonSub<F>(z3, w37));
        }
        // Scalar tail.
        for (; j < qlen8; ++j)
        {
          UInt x0 = a[i + j];
          UInt x1 = a[i + j + qlen8];
          UInt x2 = a[i + j + qlen4];
          UInt x3 = a[i + j + qlen4 + qlen8];
          UInt x4 = a[i + j + halflen];
          UInt x5 = a[i + j + halflen + qlen8];
          UInt x6 = a[i + j + halflen + qlen4];
          UInt x7 = a[i + j + halflen + qlen4 + qlen8];

          UInt g_a = roots[j * stride4];
          UInt v01 = F::Mul(x1, g_a);
          UInt y0 = F::Add(x0, v01);
          UInt y1 = F::Sub(x0, v01);
          UInt v23 = F::Mul(x3, g_a);
          UInt y2 = F::Add(x2, v23);
          UInt y3 = F::Sub(x2, v23);
          UInt v45 = F::Mul(x5, g_a);
          UInt y4 = F::Add(x4, v45);
          UInt y5 = F::Sub(x4, v45);
          UInt v67 = F::Mul(x7, g_a);
          UInt y6 = F::Add(x6, v67);
          UInt y7 = F::Sub(x6, v67);

          UInt g_b1 = roots[2 * j * stride];
          UInt g_b2 = roots[2 * (j + qlen8) * stride];
          UInt w02 = F::Mul(y2, g_b1);
          UInt z0 = F::Add(y0, w02);
          UInt z2 = F::Sub(y0, w02);
          UInt w13 = F::Mul(y3, g_b2);
          UInt z1 = F::Add(y1, w13);
          UInt z3 = F::Sub(y1, w13);
          UInt w46 = F::Mul(y6, g_b1);
          UInt z4 = F::Add(y4, w46);
          UInt z6 = F::Sub(y4, w46);
          UInt w57 = F::Mul(y7, g_b2);
          UInt z5 = F::Add(y5, w57);
          UInt z7 = F::Sub(y5, w57);

          UInt g_c0 = roots[j * stride];
          UInt g_c1 = roots[(j + qlen8) * stride];
          UInt g_c2 = roots[(j + qlen4) * stride];
          UInt g_c3 = roots[(j + qlen4 + qlen8) * stride];
          UInt w04 = F::Mul(z4, g_c0);
          UInt w15 = F::Mul(z5, g_c1);
          UInt w26 = F::Mul(z6, g_c2);
          UInt w37 = F::Mul(z7, g_c3);

          a[i + j]                          = F::Add(z0, w04);
          a[i + j + halflen]                = F::Sub(z0, w04);
          a[i + j + qlen8]                  = F::Add(z1, w15);
          a[i + j + halflen + qlen8]        = F::Sub(z1, w15);
          a[i + j + qlen4]                  = F::Add(z2, w26);
          a[i + j + halflen + qlen4]        = F::Sub(z2, w26);
          a[i + j + qlen4 + qlen8]          = F::Add(z3, w37);
          a[i + j + halflen + qlen4 + qlen8] = F::Sub(z3, w37);
        }
      }
    }
#endif // BIGMATH_NEON

    template <typename F>
    inline void ForwardRadix8Layer(UInt *a, Int n, Int outerLen, const UInt *roots)
    {
      Int halflen = outerLen >> 1;
      Int qlen4 = outerLen >> 2;
      Int qlen8 = outerLen >> 3;
      Int stride = n / outerLen;
      Int stride4 = stride << 2;
      Int omega4_off = n / 4;
      for (Int i = 0; i < n; i += outerLen)
      {
        for (Int j = 0; j < qlen8; ++j)
        {
          UInt x0 = a[i + j];
          UInt x1 = a[i + j + qlen8];
          UInt x2 = a[i + j + qlen4];
          UInt x3 = a[i + j + qlen4 + qlen8];
          UInt x4 = a[i + j + halflen];
          UInt x5 = a[i + j + halflen + qlen8];
          UInt x6 = a[i + j + halflen + qlen4];
          UInt x7 = a[i + j + halflen + qlen4 + qlen8];

          UInt g_a  = roots[j * stride];
          UInt gw_a = roots[j * stride + omega4_off];
          UInt g2_a = roots[2 * j * stride];
          UInt g_b  = roots[(j + qlen8) * stride];
          UInt gw_b = roots[(j + qlen8) * stride + omega4_off];
          UInt g2_b = roots[2 * (j + qlen8) * stride];
          UInt g3   = roots[j * stride4];

          UInt t0 = F::Add(x0, x4);
          UInt t1 = F::Add(x2, x6);
          UInt t2 = F::Mul(F::Sub(x0, x4), g_a);
          UInt t3 = F::Mul(F::Sub(x2, x6), gw_a);
          UInt y0 = F::Add(t0, t1);
          UInt y2 = F::Mul(F::Sub(t0, t1), g2_a);
          UInt y4 = F::Add(t2, t3);
          UInt y6 = F::Mul(F::Sub(t2, t3), g2_a);

          UInt u0 = F::Add(x1, x5);
          UInt u1 = F::Add(x3, x7);
          UInt u2 = F::Mul(F::Sub(x1, x5), g_b);
          UInt u3 = F::Mul(F::Sub(x3, x7), gw_b);
          UInt y1 = F::Add(u0, u1);
          UInt y3 = F::Mul(F::Sub(u0, u1), g2_b);
          UInt y5 = F::Add(u2, u3);
          UInt y7 = F::Mul(F::Sub(u2, u3), g2_b);

          UInt z0 = F::Add(y0, y1);
          UInt z1 = F::Mul(F::Sub(y0, y1), g3);
          UInt z2 = F::Add(y2, y3);
          UInt z3 = F::Mul(F::Sub(y2, y3), g3);
          UInt z4 = F::Add(y4, y5);
          UInt z5 = F::Mul(F::Sub(y4, y5), g3);
          UInt z6 = F::Add(y6, y7);
          UInt z7 = F::Mul(F::Sub(y6, y7), g3);

          a[i + j]                          = z0;
          a[i + j + qlen8]                  = z1;
          a[i + j + qlen4]                  = z2;
          a[i + j + qlen4 + qlen8]          = z3;
          a[i + j + halflen]                = z4;
          a[i + j + halflen + qlen8]        = z5;
          a[i + j + halflen + qlen4]        = z6;
          a[i + j + halflen + qlen4 + qlen8] = z7;
        }
      }
    }

    // Radix-8 fused DIT: inverse of ForwardRadix8Layer.
    template <typename F>
    inline void InverseRadix8Layer(UInt *a, Int n, Int outerLen, const UInt *roots)
    {
      Int halflen = outerLen >> 1;
      Int qlen4 = outerLen >> 2;
      Int qlen8 = outerLen >> 3;
      Int stride = n / outerLen;
      Int stride4 = stride << 2;
      for (Int i = 0; i < n; i += outerLen)
      {
        for (Int j = 0; j < qlen8; ++j)
        {
          UInt x0 = a[i + j];
          UInt x1 = a[i + j + qlen8];
          UInt x2 = a[i + j + qlen4];
          UInt x3 = a[i + j + qlen4 + qlen8];
          UInt x4 = a[i + j + halflen];
          UInt x5 = a[i + j + halflen + qlen8];
          UInt x6 = a[i + j + halflen + qlen4];
          UInt x7 = a[i + j + halflen + qlen4 + qlen8];

          // Layer A (innermost, twiddle g_a).
          UInt g_a = roots[j * stride4];
          UInt v01 = F::Mul(x1, g_a);
          UInt y0 = F::Add(x0, v01);
          UInt y1 = F::Sub(x0, v01);
          UInt v23 = F::Mul(x3, g_a);
          UInt y2 = F::Add(x2, v23);
          UInt y3 = F::Sub(x2, v23);
          UInt v45 = F::Mul(x5, g_a);
          UInt y4 = F::Add(x4, v45);
          UInt y5 = F::Sub(x4, v45);
          UInt v67 = F::Mul(x7, g_a);
          UInt y6 = F::Add(x6, v67);
          UInt y7 = F::Sub(x6, v67);

          // Layer B (middle).
          UInt g_b1 = roots[2 * j * stride];
          UInt g_b2 = roots[2 * (j + qlen8) * stride];
          UInt w02 = F::Mul(y2, g_b1);
          UInt z0 = F::Add(y0, w02);
          UInt z2 = F::Sub(y0, w02);
          UInt w13 = F::Mul(y3, g_b2);
          UInt z1 = F::Add(y1, w13);
          UInt z3 = F::Sub(y1, w13);
          UInt w46 = F::Mul(y6, g_b1);
          UInt z4 = F::Add(y4, w46);
          UInt z6 = F::Sub(y4, w46);
          UInt w57 = F::Mul(y7, g_b2);
          UInt z5 = F::Add(y5, w57);
          UInt z7 = F::Sub(y5, w57);

          // Layer C (outermost).
          UInt g_c0 = roots[j * stride];
          UInt g_c1 = roots[(j + qlen8) * stride];
          UInt g_c2 = roots[(j + qlen4) * stride];
          UInt g_c3 = roots[(j + qlen4 + qlen8) * stride];
          UInt w04 = F::Mul(z4, g_c0);
          UInt w15 = F::Mul(z5, g_c1);
          UInt w26 = F::Mul(z6, g_c2);
          UInt w37 = F::Mul(z7, g_c3);

          a[i + j]                          = F::Add(z0, w04);
          a[i + j + halflen]                = F::Sub(z0, w04);
          a[i + j + qlen8]                  = F::Add(z1, w15);
          a[i + j + halflen + qlen8]        = F::Sub(z1, w15);
          a[i + j + qlen4]                  = F::Add(z2, w26);
          a[i + j + halflen + qlen4]        = F::Sub(z2, w26);
          a[i + j + qlen4 + qlen8]          = F::Add(z3, w37);
          a[i + j + halflen + qlen4 + qlen8] = F::Sub(z3, w37);
        }
      }
    }

    // Pointer-based leaf forward — radix-8/4/2 fused chain, in-place, no scaling.
    // Used by both the vector wrapper below and by MFA sub-FFT calls.
    template <typename F>
    inline void ForwardPtr(UInt *aPtr, Int n, const Plan<F> &plan)
    {
      if (n <= 1) return;
      const UInt *roots = plan.forwardRoots.data();
#if BIGMATH_NEON
      const UInt *fil = plan.forwardIl.data();
      if (n >= 16 && n <= BIGMATH_NEON_NTT_MAX)
      {
        // All-NEON chain: radix-8 layers while the j-loop has >= 4 lanes
        // (len >= 32), then the tail kernel for whatever is left.
        Int l = n;
        while (l >= 32)
        {
          ForwardRadix8LayerNeon<F>(aPtr, n, l, roots, fil);
          l >>= 3;
        }
        if (l == 16)
          ForwardTail16Neon<F>(aPtr, n, roots);
        else if (l == 8)
          ForwardTail8Neon<F>(aPtr, n, roots);
        else
          ForwardRadix4LayerNeon<F>(aPtr, n, roots); // l == 4
        return;
      }
#endif
      Int len = n;
      while (len >= 8)
      {
#if BIGMATH_NEON
        if ((len >> 3) >= 4 && n <= BIGMATH_NEON_NTT_MAX)
          ForwardRadix8LayerNeon<F>(aPtr, n, len, roots, fil);
        else
#endif
        ForwardRadix8Layer<F>(aPtr, n, len, roots);
        len >>= 3;
      }
      if (len == 4)
      {
#if BIGMATH_NEON
        ForwardRadix4LayerNeon<F>(aPtr, n, roots);
#else
        ForwardRadix4Layer<F>(aPtr, n, 4, roots);
#endif
      }
      else if (len == 2)
        ForwardRadix2Layer<F>(aPtr, n, 2, roots);
    }

    // Pointer-based leaf inverse. `scale` controls whether the final 1/n
    // multiplication is applied — MFA sub-inverses leave it on so the per-leaf
    // 1/leaf_n scalings compose to 1/N total at the top level.
    template <typename F>
    inline void InversePtr(UInt *aPtr, Int n, const Plan<F> &plan, bool scale)
    {
      const UInt *roots = plan.inverseRoots.data();
#if BIGMATH_NEON
      const UInt *iil = plan.inverseIl.data();
      if (n >= 16 && n <= BIGMATH_NEON_NTT_MAX)
      {
        // All-NEON chain, mirror of ForwardPtr: tail kernel first, then
        // radix-8 layers (always len >= 32 here).
        Int rem = __builtin_ctz((unsigned)n) % 3;
        Int l;
        if (rem == 2)
        {
          InverseRadix4LayerNeon<F>(aPtr, n, roots);
          l = 32;
        }
        else if (rem == 1)
        {
          InverseTail16Neon<F>(aPtr, n, roots);
          l = 128;
        }
        else
        {
          InverseTail8Neon<F>(aPtr, n, roots);
          l = 64;
        }
        while (l <= n)
        {
          InverseRadix8LayerNeon<F>(aPtr, n, l, roots, iil);
          l <<= 3;
        }
      }
      else
#endif
      if (n >= 2)
      {
        Int logn = __builtin_ctz((unsigned)n);
        Int rem = logn % 3;
        Int len;
        if (rem == 1)
        {
#if BIGMATH_NEON
          InverseRadix2LayerNeon<F>(aPtr, n);
#else
          InverseRadix2Layer<F>(aPtr, n, 2, roots);
#endif
          len = 16;
        }
        else if (rem == 2)
        {
#if BIGMATH_NEON
          InverseRadix4LayerNeon<F>(aPtr, n, roots);
#else
          InverseRadix4Layer<F>(aPtr, n, 4, roots);
#endif
          len = 32;
        }
        else
        {
          len = 8;
        }
        while (len <= n)
        {
#if BIGMATH_NEON
          if ((len >> 3) >= 4 && n <= BIGMATH_NEON_NTT_MAX)
            InverseRadix8LayerNeon<F>(aPtr, n, len, roots, iil);
          else
#endif
          InverseRadix8Layer<F>(aPtr, n, len, roots);
          len <<= 3;
        }
      }

      if (scale)
      {
        UInt invSize = plan.invSize;
#if BIGMATH_NEON
        NeonScale<F>(aPtr, n, invSize);
#else
        for (Int i = 0; i < n; ++i) aPtr[i] = F::Mul(aPtr[i], invSize);
#endif
      }
    }

    // ─── Split inverse: two half-length inverses + a shared top layer ───────
    // A DIT inverse of length n is its first log2(n) − 1 layers applied to
    // each half independently — exactly a length-n/2 inverse with the n/2
    // plan (ω_{n/2}^k = ω_n^{2k}) — followed by the single len == n layer
    // that combines the halves. The halves are independent work units, and
    // the top layer splits freely over j, so one prime's inverse occupies two
    // threads instead of one. The 1/n scaling rides the top layer's pass.

    // Top DIT layer (len == n) fused with the 1/n scaling, for j in [j0, j1).
    template <typename F>
    inline void InverseTopLayerScaled(UInt *a, Int n, const Plan<F> &plan, Int j0, Int j1)
    {
      const Int h = n >> 1;
      const UInt *roots = plan.inverseRoots.data();
      const UInt invN = plan.invSize;
      Int j = j0;
#if BIGMATH_NEON
      const UInt *il = plan.inverseIl.data();
      const uint32x4_t vc = vdupq_n_u32(invN);
      const uint32x4_t vcs = vdupq_n_u32(ShoupCompanion<F>(invN));
      for (; j + 4 <= j1; j += 4)
      {
        uint32x4x2_t tw = vld2q_u32(il + 2 * j); // val[0] = twiddle, val[1] = companion
        uint32x4_t u = vld1q_u32(a + j);
        uint32x4_t v = NeonShoupMul<F>(tw.val[0], tw.val[1], vld1q_u32(a + j + h));
        vst1q_u32(a + j, NeonShoupMul<F>(vc, vcs, NeonAdd<F>(u, v)));
        vst1q_u32(a + j + h, NeonShoupMul<F>(vc, vcs, NeonSub<F>(u, v)));
      }
#endif
      for (; j < j1; ++j)
      {
        UInt u = a[j];
        UInt v = F::Mul(a[j + h], roots[j]);
        a[j] = F::Mul(F::Add(u, v), invN);
        a[j + h] = F::Mul(F::Sub(u, v), invN);
      }
    }

    // Serial form of the split inverse; equals InversePtr(a, n, plan, true).
    // `halfPlan` must be the plan for n / 2. Requires n >= 4.
    template <typename F>
    inline void InverseSplitPtr(UInt *a, Int n, const Plan<F> &plan, const Plan<F> &halfPlan)
    {
      const Int h = n >> 1;
      InversePtr<F>(a, h, halfPlan, /*scale=*/ false);
      InversePtr<F>(a + h, h, halfPlan, /*scale=*/ false);
      InverseTopLayerScaled<F>(a, n, plan, 0, h);
    }

    template <typename F>
    inline void Forward(std::vector<UInt> &a, const Plan<F> &plan)
    {
      ForwardPtr<F>(a.data(), (Int)a.size(), plan);
    }

    template <typename F>
    inline void Inverse(std::vector<UInt> &a, const Plan<F> &plan)
    {
      InversePtr<F>(a.data(), (Int)a.size(), plan, /*scale=*/ true);
    }

    // ─── Garner CRT reconstruction ───────────────────────────────────────────
    //
    // Given residues r1, r2, r3 with r_i = x mod p_i, recover x in
    // [0, p1·p2·p3). Output as (low64, mid64, high64) for downstream carry
    // propagation; total bit-width ≈ 90, fits two 64-bit words.

    struct InvTable
    {
      UInt p1_inv_mod_p2;  // p1^-1 mod p2
      UInt p1_inv_mod_p3;  // p1^-1 mod p3
      UInt p2_inv_mod_p3;  // p2^-1 mod p3
      ULong p1_long = P1;
      ULong128 p1p2 = (ULong128)P1 * P2;
    };

    inline const InvTable &GarnerInverses()
    {
      static const InvTable inv = []() {
        InvTable t;
        t.p1_inv_mod_p2 = F2::Inv(P1 % P2);
        t.p1_inv_mod_p3 = F3::Inv(P1 % P3);
        t.p2_inv_mod_p3 = F3::Inv(P2 % P3);
        return t;
      }();
      return inv;
    }

    // Reconstruct full convolution coefficient from three residues.
    // Returns the value as ULong128 (fits in ~90 bits comfortably).
    inline ULong128 Garner(UInt r1, UInt r2, UInt r3, const InvTable &inv)
    {
      // u1 = r1                                  (mod p1)
      // u2 = (r2 - r1) * p1^-1                   (mod p2)
      // u3 = ((r3 - r1) * p1^-1 - u2) * p2^-1    (mod p3)
      // x  = u1 + p1*u2 + p1*p2*u3
      UInt u1 = r1;

      UInt diff2 = F2::Sub(r2 % P2, r1 % P2);
      UInt u2 = F2::Mul(diff2, inv.p1_inv_mod_p2);

      UInt diff3 = F3::Sub(r3 % P3, r1 % P3);
      UInt tmp = F3::Mul(diff3, inv.p1_inv_mod_p3);
      UInt u3 = F3::Mul(F3::Sub(tmp, u2 % P3), inv.p2_inv_mod_p3);

      return (ULong128)u1 + (ULong128)P1 * u2 + inv.p1p2 * u3;
    }

    inline SizeT CoeffsPerLimb(BaseT base)
    {
      return base == Base2_64 ? 2u : 1u;
    }

    inline void PackOperand(const std::vector<DataT> &v,
                            BaseT base,
                            std::vector<UInt> &dst1,
                            std::vector<UInt> &dst2,
                            std::vector<UInt> &dst3)
    {
      if (base == Base2_64)
      {
        for (SizeT i = 0; i < v.size(); ++i)
        {
          SizeT j = i * 2;
          UInt lo = (UInt)(v[i] & 0xFFFFFFFFULL);
          UInt hi = (UInt)((v[i] >> 32) & 0xFFFFFFFFULL);
          dst1[j]     = lo % P1; dst1[j + 1] = hi % P1;
          dst2[j]     = lo % P2; dst2[j + 1] = hi % P2;
          dst3[j]     = lo % P3; dst3[j + 1] = hi % P3;
        }
      }
      else
      {
        // Base2_32: each limb is already a 32-bit value.
        for (SizeT i = 0; i < v.size(); ++i)
        {
          UInt vv = (UInt)v[i];
          dst1[i] = vv % P1;
          dst2[i] = vv % P2;
          dst3[i] = vv % P3;
        }
      }
    }

    // Residues of one operand modulo one prime: the packed 32-bit coefficients
    // in dst[0 .. used), zero padding in dst[used .. n). One plane at a time
    // so each (operand, prime) forward unit can pack its own input.
    template <typename F>
    inline void PackPlaneLimbs(const DataT *v, SizeT limbs, BaseT base, UInt *dst, Int n)
    {
      SizeT used;
      if (base == Base2_64)
      {
        for (SizeT i = 0; i < limbs; ++i)
        {
          dst[2 * i] = (UInt)(v[i] & 0xFFFFFFFFULL) % F::Prime;
          dst[2 * i + 1] = (UInt)((v[i] >> 32) & 0xFFFFFFFFULL) % F::Prime;
        }
        used = 2 * limbs;
      }
      else
      {
        for (SizeT i = 0; i < limbs; ++i)
          dst[i] = (UInt)v[i] % F::Prime;
        used = limbs;
      }
      if ((SizeT)n > used)
        std::memset(dst + used, 0, ((SizeT)n - used) * sizeof(UInt));
    }

    template <typename F>
    inline void PackPlane(const std::vector<DataT> &v, BaseT base, UInt *dst, Int n)
    {
      PackPlaneLimbs<F>(v.data(), (SizeT)v.size(), base, dst, n);
    }

    // ─── Garner + carry for 64-bit limbs (NEON digits, chunk-parallel) ──────
    //
    // The finalize pass turns the three residue planes into limbs: Garner
    // gives each coefficient's value x_i < p1·p2·p3 ≈ 2^90.5, and the product
    // is sum_i x_i · 2^(32 i). It used to be one scalar loop on the calling
    // thread — eight `% P` reductions per coefficient and a serial carry
    // chain — which after the transforms went parallel was 17-33% of a
    // threaded multiply's wall-clock (27 ms of 100 ms at n = 2^22).
    //
    // Two changes:
    //   * The mixed-radix digits u2, u3 are computed four coefficients at a
    //     time with constant-multiplier Shoup products (NeonGarnerDigits);
    //     only the 96-bit assembly x = r1 + p1·u2 + p1·p2·u3 and the carry
    //     stay scalar.
    //   * The carry chain is cut into per-thread limb ranges. Each range
    //     runs with carry-in 0 and reports its carry-out; because the sum is
    //     linear, adding carry-out_k at the first limb of range k+1 (with an
    //     ordinary ripple) reconstructs the exact result.

    static_assert((ULong)5 * P2 > P1, "Garner: r2 + 5*P2 - r1 must stay positive");
    static_assert((ULong)6 * P2 < ((ULong)1 << 32), "Garner: r2 + 5*P2 - r1 must fit 32 bits");
    static_assert((ULong)3 * P2 < ((ULong)1 << 31), "Garner: reduced difference must be < 2^31 for NeonShoupMul");
    static_assert(P1 < (ULong)2 * P3, "Garner: r1 mod P3 is one conditional subtract");
    static_assert(P2 < P3, "Garner: u2 is already reduced mod P3");
    static_assert((ULong128)P1 * P2 < ((ULong128)1 << 63), "Garner: r1 + p1*u2 fits 64 bits");

#if BIGMATH_NEON
    struct GarnerNeonConsts
    {
      uint32x4_t fiveP2, threeP2, p3;
      uint32x4_t inv12, inv12s, inv13, inv13s, inv23, inv23s;

      explicit GarnerNeonConsts(const InvTable &inv)
          : fiveP2(vdupq_n_u32(5u * P2)), threeP2(vdupq_n_u32(3u * P2)), p3(vdupq_n_u32(P3)),
            inv12(vdupq_n_u32(inv.p1_inv_mod_p2)), inv12s(vdupq_n_u32(ShoupCompanion<F2>(inv.p1_inv_mod_p2))),
            inv13(vdupq_n_u32(inv.p1_inv_mod_p3)), inv13s(vdupq_n_u32(ShoupCompanion<F3>(inv.p1_inv_mod_p3))),
            inv23(vdupq_n_u32(inv.p2_inv_mod_p3)), inv23s(vdupq_n_u32(ShoupCompanion<F3>(inv.p2_inv_mod_p3)))
      {
      }
    };

    // Mixed-radix digits for four coefficients:
    //   u2 = (r2 − r1) · p1^-1                (mod p2)
    //   u3 = ((r3 − r1) · p1^-1 − u2) · p2^-1 (mod p3)
    // Bit-identical to the scalar Garner() digits.
    static inline void NeonGarnerDigits(uint32x4_t r1, uint32x4_t r2, uint32x4_t r3,
                                        const GarnerNeonConsts &k, uint32x4_t &u2, uint32x4_t &u3)
    {
      // r2 + 5·P2 − r1 lies in (0, 6·P2) — no reduction of r1 mod P2 needed.
      // One conditional subtract of 3·P2 brings it below 3·P2 < 2^31, the
      // multiplicand bound NeonShoupMul needs.
      uint32x4_t d2 = vsubq_u32(vaddq_u32(r2, k.fiveP2), r1);
      d2 = vminq_u32(d2, vsubq_u32(d2, k.threeP2));
      u2 = NeonShoupMul<F2>(k.inv12, k.inv12s, d2);

      uint32x4_t r1m3 = vminq_u32(r1, vsubq_u32(r1, k.p3)); // r1 mod P3
      uint32x4_t t3 = NeonShoupMul<F3>(k.inv13, k.inv13s, NeonSub<F3>(r3, r1m3));
      u3 = NeonShoupMul<F3>(k.inv23, k.inv23s, NeonSub<F3>(t3, u2));
    }
#endif

    // Limbs [ls, le) of sum_i x_i · 2^(32 i) restricted to coefficients
    // 2·ls .. 2·le − 1, computed with carry-in 0. Coefficients at index
    // ≥ coeffCount are zero (only the last limb of an odd-length product).
    // Returns the carry out of limb le − 1; it is < 2^59.
    inline ULong GarnerCarryRange64(const UInt *f1, const UInt *f2, const UInt *f3,
                                    SizeT coeffCount, DataT *out, SizeT ls, SizeT le,
                                    const InvTable &inv)
    {
      ULong carry = 0;
      SizeT i = ls;
#if BIGMATH_NEON
      const GarnerNeonConsts k(inv);
      const ULong p1 = P1;
      const ULong p1p2 = (ULong)inv.p1p2;
      for (; i + 2 <= le && (ULong)2 * i + 4 <= coeffCount; i += 2)
      {
        const SizeT c = 2 * i;
        uint32x4_t u2v, u3v;
        NeonGarnerDigits(vld1q_u32(f1 + c), vld1q_u32(f2 + c), vld1q_u32(f3 + c), k, u2v, u3v);
        UInt u2[4], u3[4];
        vst1q_u32(u2, u2v);
        vst1q_u32(u3, u3v);
        for (SizeT h = 0; h < 2; ++h)
        {
          ULong128 x0 = (ULong128)p1p2 * u3[2 * h] + ((ULong)f1[c + 2 * h] + p1 * u2[2 * h]);
          ULong128 x1 = (ULong128)p1p2 * u3[2 * h + 1] + ((ULong)f1[c + 2 * h + 1] + p1 * u2[2 * h + 1]);
          ULong128 t = x0 + (x1 << 32) + carry; // < 2^123
          out[i + h] = (DataT)t;
          carry = (ULong)(t >> 64);
        }
      }
#endif
      for (; i < le; ++i)
      {
        const ULong c = (ULong)2 * i;
        ULong128 x0 = c < coeffCount ? Garner(f1[c], f2[c], f3[c], inv) : (ULong128)0;
        ULong128 x1 = c + 1 < coeffCount ? Garner(f1[c + 1], f2[c + 1], f3[c + 1], inv) : (ULong128)0;
        ULong128 t = x0 + (x1 << 32) + carry;
        out[i] = (DataT)t;
        carry = (ULong)(t >> 64);
      }
      return carry;
    }

    // A carry chain over out[0 .. nLimbs) cut into per-thread limb ranges.
    // range(ls, le) writes limbs [ls, le) with carry-in 0 and returns the
    // carry out of limb le − 1. Returns the carry out of the top limb. Splits
    // across the pool from BIGMATH_NTT_FINALIZE_PARALLEL_MIN limbs.
    template <typename Range>
    inline ULong CarryChunks64(DataT *out, SizeT nLimbs, const Range &range)
    {
      constexpr SizeT MaxChunks = 64;
      SizeT chunks = 1;
      if (nLimbs >= (SizeT)BIGMATH_NTT_FINALIZE_PARALLEL_MIN)
        chunks = std::min<SizeT>({(SizeT)ParallelNumThreads(), MaxChunks,
                                  nLimbs / ((SizeT)BIGMATH_NTT_FINALIZE_PARALLEL_MIN / 4)});
      if (chunks <= 1)
        return range((SizeT)0, nLimbs);

      // Even range length keeps every range on a 4-coefficient NEON boundary.
      const SizeT per = (((nLimbs + chunks - 1) / chunks) + 1) & ~(SizeT)1;
      chunks = (nLimbs + per - 1) / per;
      ULong carries[MaxChunks];
      ULong *carriesPtr = carries;
      ParallelDo((Int)chunks, [&range, carriesPtr, per, nLimbs](Int s, Int e) {
        for (Int kIdx = s; kIdx < e; ++kIdx)
        {
          SizeT ls = (SizeT)kIdx * per;
          SizeT le = std::min(nLimbs, ls + per);
          carriesPtr[kIdx] = range(ls, le);
        }
      });

      // Stitch: add each range's carry-out at the next range's first limb.
      ULong top = carries[chunks - 1];
      for (SizeT kIdx = 0; kIdx + 1 < chunks; ++kIdx)
      {
        ULong c = carries[kIdx];
        for (SizeT i = (kIdx + 1) * per; c != 0; ++i)
        {
          if (i == nLimbs)
          {
            top += c; // a ripple that ran off the top (c is 1 here)
            break;
          }
          ULong sum = (ULong)out[i] + c;
          c = sum < c ? 1 : 0;
          out[i] = (DataT)sum;
        }
      }
      return top;
    }

    // out[0 .. nLimbs) = the low nLimbs limbs of sum_i x_i · 2^(32 i); returns
    // the carry out of the top limb.
    inline ULong GarnerCarry64(const UInt *f1, const UInt *f2, const UInt *f3,
                               SizeT coeffCount, DataT *out, SizeT nLimbs)
    {
      const InvTable &inv = GarnerInverses();
      return CarryChunks64(out, nLimbs, [=, &inv](SizeT ls, SizeT le) {
        return GarnerCarryRange64(f1, f2, f3, coeffCount, out, ls, le, inv);
      });
    }

    inline std::vector<DataT> FinalizeProduct(const std::vector<UInt> &fa1,
                                              const std::vector<UInt> &fa2,
                                              const std::vector<UInt> &fa3,
                                              ULong coeffCount,
                                              BaseT base,
                                              SizeT reserveLimbs)
    {
      if (base == Base2_64)
      {
        // Two 32-bit coefficients per 64-bit limb. The product needs
        // ceil(coeffCount / 2) limbs plus whatever carry leaves the top.
        const SizeT nLimbs = (SizeT)((coeffCount + 1) / 2);
        std::vector<DataT> result;
        result.reserve(std::max<SizeT>(reserveLimbs, nLimbs + 1));
        result.resize(nLimbs);
        ULong carry = GarnerCarry64(fa1.data(), fa2.data(), fa3.data(), (SizeT)coeffCount,
                                    result.data(), nLimbs);
        if (carry != 0)
          result.push_back((DataT)carry);
        TrimZeros(result);
        return result;
      }

      const InvTable &inv = GarnerInverses();
      std::vector<DataT> result;
      result.reserve(reserveLimbs);
      {
        // Base2_32: 1 coefficient per output limb (already 32-bit).
        ULong128 carry = 0;
        for (SizeT i = 0; i < (SizeT)coeffCount; ++i)
        {
          ULong128 total = Garner(fa1[i], fa2[i], fa3[i], inv) + carry;
          result.push_back((DataT)(total & 0xFFFFFFFFULL));
          carry = total >> 32;
        }
        while (carry > 0)
        {
          result.push_back((DataT)(carry & 0xFFFFFFFFULL));
          carry >>= 32;
        }
      }

      TrimZeros(result);
      return result;
    }


    // ─── Shared whole-transform batch orchestration ─────────────────────────
    // The plain (non-MFA) forward/pointwise/inverse ParallelDo blocks were
    // duplicated across Multiply, MultiplyMod2km1 and the prepared-operand
    // path; any future edit must happen here once.

    // Forward transforms: always the three a-side buffers; the three b-side
    // buffers too when fb1 != nullptr (one 6-unit batch). Squaring and the
    // prepared paths pass fb1 == nullptr for a 3-unit batch.
    inline void CrtForwardBatch(std::vector<UInt> &fa1, std::vector<UInt> &fa2, std::vector<UInt> &fa3,
                                std::vector<UInt> *fb1, std::vector<UInt> *fb2, std::vector<UInt> *fb3,
                                const Plan<F1> &plan1, const Plan<F2> &plan2, const Plan<F3> &plan3)
    {
#if BIGMATH_USE_THREADS
      if ((Int)fa1.size() < (Int)BIGMATH_NTT_PARALLEL_MIN)
      {
        Forward<F1>(fa1, plan1);
        Forward<F2>(fa2, plan2);
        Forward<F3>(fa3, plan3);
        if (fb1 != nullptr)
        {
          Forward<F1>(*fb1, plan1);
          Forward<F2>(*fb2, plan2);
          Forward<F3>(*fb3, plan3);
        }
        return;
      }
      std::vector<UInt> *bufs[6] = {&fa1, fb1, &fa2, fb2, &fa3, fb3};
      const Plan<F1> *p1 = &plan1;
      const Plan<F2> *p2 = &plan2;
      const Plan<F3> *p3 = &plan3;
      if (fb1 == nullptr)
      {
        auto body = [bufs, p1, p2, p3](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx)
          {
            switch (idx)
            {
              case 0: Forward<F1>(*bufs[0], *p1); break;
              case 1: Forward<F2>(*bufs[2], *p2); break;
              case 2: Forward<F3>(*bufs[4], *p3); break;
            }
          }
        };
        ParallelDo(3, body);
      }
      else
      {
        auto body = [bufs, p1, p2, p3](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx)
          {
            switch (idx)
            {
              case 0: Forward<F1>(*bufs[0], *p1); break;
              case 1: Forward<F1>(*bufs[1], *p1); break;
              case 2: Forward<F2>(*bufs[2], *p2); break;
              case 3: Forward<F2>(*bufs[3], *p2); break;
              case 4: Forward<F3>(*bufs[4], *p3); break;
              case 5: Forward<F3>(*bufs[5], *p3); break;
            }
          }
        };
        ParallelDo(6, body);
      }
#else
      Forward<F1>(fa1, plan1);
      Forward<F2>(fa2, plan2);
      Forward<F3>(fa3, plan3);
      if (fb1 != nullptr)
      {
        Forward<F1>(*fb1, plan1);
        Forward<F2>(*fb2, plan2);
        Forward<F3>(*fb3, plan3);
      }
#endif
    }

    // Pack + forward in one batch: each (operand, prime) unit packs its own
    // residue plane (and zero-pads it) before transforming it, so the pack and
    // the padding run on the pool instead of serially on the caller. `b` is
    // null when squaring (three units instead of six). The buffers are sized
    // here; every element of each plane is written by its unit.
    inline void CrtPackForwardBatch(const std::vector<DataT> &a, const std::vector<DataT> *b,
                                    BaseT base, Int n,
                                    std::vector<UInt> &fa1, std::vector<UInt> &fa2, std::vector<UInt> &fa3,
                                    std::vector<UInt> &fb1, std::vector<UInt> &fb2, std::vector<UInt> &fb3,
                                    const Plan<F1> &plan1, const Plan<F2> &plan2, const Plan<F3> &plan3)
    {
      fa1.resize(n);
      fa2.resize(n);
      fa3.resize(n);
      if (b != nullptr)
      {
        fb1.resize(n);
        fb2.resize(n);
        fb3.resize(n);
      }
      UInt *const pa1 = fa1.data(), *const pa2 = fa2.data(), *const pa3 = fa3.data();
      UInt *const pb1 = b ? fb1.data() : nullptr, *const pb2 = b ? fb2.data() : nullptr, *const pb3 = b ? fb3.data() : nullptr;
      const std::vector<DataT> *const srcA = &a;
      const Plan<F1> *const p1 = &plan1;
      const Plan<F2> *const p2 = &plan2;
      const Plan<F3> *const p3 = &plan3;
      auto unit = [=](Int idx) {
        switch (idx)
        {
          case 0: PackPlane<F1>(*srcA, base, pa1, n); ForwardPtr<F1>(pa1, n, *p1); break;
          case 1: PackPlane<F2>(*srcA, base, pa2, n); ForwardPtr<F2>(pa2, n, *p2); break;
          case 2: PackPlane<F3>(*srcA, base, pa3, n); ForwardPtr<F3>(pa3, n, *p3); break;
          case 3: PackPlane<F1>(*b, base, pb1, n); ForwardPtr<F1>(pb1, n, *p1); break;
          case 4: PackPlane<F2>(*b, base, pb2, n); ForwardPtr<F2>(pb2, n, *p2); break;
          case 5: PackPlane<F3>(*b, base, pb3, n); ForwardPtr<F3>(pb3, n, *p3); break;
        }
      };
      const Int units = b != nullptr ? 6 : 3;
#if BIGMATH_USE_THREADS
      if (n >= (Int)BIGMATH_NTT_PARALLEL_MIN)
      {
        ParallelDo(units, [&unit](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx) unit(idx);
        });
        return;
      }
#endif
      for (Int idx = 0; idx < units; ++idx) unit(idx);
    }

    // Pointwise product into the a-side buffers. Pass the a-side pointers as
    // the b-side to square in place.
    inline void CrtPointwiseMulInto(UInt *p1a, UInt *p2a, UInt *p3a,
                                    const UInt *p1b, const UInt *p2b, const UInt *p3b,
                                    Int n)
    {
      auto body = [p1a, p1b, p2a, p2b, p3a, p3b](Int s, Int e) {
        for (Int i = s; i < e; ++i)
        {
          p1a[i] = F1::Mul(p1a[i], p1b[i]);
          p2a[i] = F2::Mul(p2a[i], p2b[i]);
          p3a[i] = F3::Mul(p3a[i], p3b[i]);
        }
      };
      if ((SizeT)n >= ParallelMinSize()) ParallelFor(n, body);
      else body(0, n);
    }

    // Inverse transforms of the three (already pointwise-multiplied) buffers.
    inline void CrtInverseBatch(std::vector<UInt> &f1, std::vector<UInt> &f2, std::vector<UInt> &f3,
                                const Plan<F1> &plan1, const Plan<F2> &plan2, const Plan<F3> &plan3)
    {
#if BIGMATH_USE_THREADS
      const Int n = (Int)f1.size();
      if (n < (Int)BIGMATH_NTT_PARALLEL_MIN)
      {
        Inverse<F1>(f1, plan1);
        Inverse<F2>(f2, plan2);
        Inverse<F3>(f3, plan3);
        return;
      }
      const Plan<F1> *p1 = &plan1;
      const Plan<F2> *p2 = &plan2;
      const Plan<F3> *p3 = &plan3;
      if (n >= (Int)BIGMATH_NTT_SPLIT_INVERSE_MIN && ParallelNumThreads() >= 4)
      {
        // Six units per phase instead of three whole transforms — see
        // InverseSplitPtr. The half plans are fetched here, on the
        // dispatching thread (plans are acquired before dispatch).
        const Int h = n >> 1;
        const Plan<F1> *h1 = &GetPlan<F1, G1>(h);
        const Plan<F2> *h2 = &GetPlan<F2, G2>(h);
        const Plan<F3> *h3 = &GetPlan<F3, G3>(h);
        UInt *const a1 = f1.data(), *const a2 = f2.data(), *const a3 = f3.data();
        ParallelDo(6, [=](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx)
          {
            const Int off = (idx & 1) ? h : 0;
            switch (idx >> 1)
            {
              case 0: InversePtr<F1>(a1 + off, h, *h1, /*scale=*/ false); break;
              case 1: InversePtr<F2>(a2 + off, h, *h2, /*scale=*/ false); break;
              case 2: InversePtr<F3>(a3 + off, h, *h3, /*scale=*/ false); break;
            }
          }
        });
        // The barrier between the phases is load-bearing: the top layer
        // reads both finished halves of its prime.
        const Int mid = (h >> 1) & ~(Int)3; // keep the NEON 4-lane groups aligned
        ParallelDo(6, [=](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx)
          {
            const Int j0 = (idx & 1) ? mid : 0;
            const Int j1 = (idx & 1) ? h : mid;
            switch (idx >> 1)
            {
              case 0: InverseTopLayerScaled<F1>(a1, n, *p1, j0, j1); break;
              case 1: InverseTopLayerScaled<F2>(a2, n, *p2, j0, j1); break;
              case 2: InverseTopLayerScaled<F3>(a3, n, *p3, j0, j1); break;
            }
          }
        });
        return;
      }
      std::vector<UInt> *bufs[3] = {&f1, &f2, &f3};
      auto body = [bufs, p1, p2, p3](Int s, Int e) {
        for (Int idx = s; idx < e; ++idx)
        {
          switch (idx)
          {
            case 0: Inverse<F1>(*bufs[0], *p1); break;
            case 1: Inverse<F2>(*bufs[1], *p2); break;
            case 2: Inverse<F3>(*bufs[2], *p3); break;
          }
        }
      };
      ParallelDo(3, body);
#else
      Inverse<F1>(f1, plan1);
      Inverse<F2>(f2, plan2);
      Inverse<F3>(f3, plan3);
#endif
    }

    // ─── Truncated length: product mod (x^H + 1) and mod (x^L − 1) ─────────
    //
    // Write the product C (c coefficients, H < c ≤ H + L) as
    // C = C_lo + x^H · C_hi with deg C_lo < H and deg C_hi < L. Then
    //
    //   C1 = C mod (x^H + 1):  C1[j] = C_lo[j] − C_hi[j]   (C_hi[j] = 0 for j ≥ L)
    //   C2 = C mod (x^L − 1):  C2[i] = Σ_k C_lo[i + kL] + C_hi[i]     (L | H)
    //
    // so with F[i] = Σ_k C1[i + kL] = Σ_k C_lo[i + kL] − C_hi[i]:
    //
    //   C_hi[i] = (C2[i] − F[i]) / 2,      C_lo[j] = C1[j] + C_hi[j].
    //
    // C1 is a negacyclic product: weight both inputs by ψ^j (ψ a primitive
    // 2H-th root, ψ^H = −1), take the ordinary length-H cyclic product, and
    // unweight by ψ^−j. The weights ψ^±j are exactly the forward/inverse
    // twiddle tables of the length-2H plan the full transform would have
    // used, so no new tables are needed. C2 is a plain length-L cyclic
    // product of the inputs folded modulo x^L − 1.
    //
    // Each prime's plane is laid out [ lower : H | cyc : L ], which after
    // reconstruction is [ C_lo | C_hi ] — the product's coefficients in
    // order, ready for the Garner/carry pass.

    // L for a product of `coeffCount` coefficients against H, or 0 when the
    // count is outside the band H < c ≤ H + (H >> bandShift). No limits on H
    // itself — the two callers below apply their own.
    inline Int TruncBandCyclicLength(ULong coeffCount, Int H, int bandShift)
    {
      if (coeffCount <= (ULong)H)
        return 0;
      const ULong over = coeffCount - (ULong)H;
      if (over > (ULong)(H >> bandShift))
        return 0;
      return std::max<Int>((Int)std::bit_ceil(over), 64);
    }

    // Cyclic-part length L for the H + L scheme, or 0 to use the full
    // power-of-two transform. H = bit_ceil(coeffCount) / 2.
    inline Int TruncCyclicLength(ULong coeffCount, Int H)
    {
#if BIGMATH_NTT_TRUNC
      if (H < (Int)BIGMATH_NTT_TRUNC_MIN || H > (Int)BIGMATH_NTT_TRUNC_MAX)
        return 0;
      return TruncBandCyclicLength(coeffCount, H, BIGMATH_NTT_TRUNC_BAND_SHIFT);
#else
      (void)coeffCount;
      (void)H;
      return 0;
#endif
    }

    // Cyclic-part length L for the single-thread cache-blocked variant
    // (CrtTruncSerial), or 0 when it does not apply: more than one thread, or
    // H outside (BIGMATH_NTT_TRUNC_MAX, BIGMATH_NTT_TRUNC_SERIAL_MAX].
    inline Int TruncSerialCyclicLength(ULong coeffCount, Int H)
    {
#if BIGMATH_NTT_TRUNC && BIGMATH_NTT_MFA
      if (H <= (Int)BIGMATH_NTT_TRUNC_MAX || H > (Int)BIGMATH_NTT_TRUNC_SERIAL_MAX)
        return 0;
      if (ParallelNumThreads() > 1)
        return 0;
      return TruncBandCyclicLength(coeffCount, H, BIGMATH_NTT_TRUNC_SERIAL_BAND_SHIFT);
#else
      (void)coeffCount;
      (void)H;
      return 0;
#endif
    }

    // dst[0..n) = (dst[0..n) + src[0..n)) mod P.
    template <typename F>
    inline void AddModInto(UInt *dst, const UInt *src, SizeT n)
    {
      SizeT i = 0;
#if BIGMATH_NEON
      for (; i + 4 <= n; i += 4)
        vst1q_u32(dst + i, NeonAdd<F>(vld1q_u32(dst + i), vld1q_u32(src + i)));
#endif
      for (; i < n; ++i)
        dst[i] = F::Add(dst[i], src[i]);
    }

    // The large serial H + L transform needs only the negacyclic weights
    // from its length-2H plan. For j = 256*b + t, store ψ^(256*b) and ψ^t.
    template <typename F>
    struct FactoredWeights
    {
      static constexpr SizeT Block = 256;
#if BIGMATH_NEON
      static constexpr SizeT Stride = 2;  // interleaved Shoup companion
#else
      static constexpr SizeT Stride = 1;
#endif
      std::vector<UInt> offsets, blocks;

      static std::vector<UInt> Powers(UInt step, SizeT count)
      {
        std::vector<UInt> table(count * Stride);
        UInt value = 1;
        for (SizeT i = 0; i < count; ++i)
        {
          table[Stride * i] = value;
#if BIGMATH_NEON
          table[Stride * i + 1] = ShoupCompanion<F>(value);
#endif
          value = F::Mul(value, step);
        }
        return table;
      }

      FactoredWeights(UInt root, SizeT h)
          : offsets(Powers(root, std::min(Block, h))),
            blocks(Powers(F::Power(root, Block), (h + Block - 1) / Block)) {}

      UInt Apply(SizeT j, UInt x) const
      {
        return F::Mul(F::Mul(x, offsets[Stride * (j % Block)]), blocks[Stride * (j / Block)]);
      }
#if BIGMATH_NEON
      // Callers advance four aligned coefficients at a time within a block.
      uint32x4_t Apply(SizeT j, uint32x4_t x) const
      {
        const auto offset = vld2q_u32(offsets.data() + Stride * (j % Block));
        const UInt *block = blocks.data() + Stride * (j / Block);
        return NeonShoupMul<F>(vdupq_n_u32(block[0]), vdupq_n_u32(block[1]),
                               NeonShoupMul<F>(offset.val[0], offset.val[1], x));
      }
#endif
    };

    template <typename F>
    struct FullWeights
    {
      const UInt *roots;
#if BIGMATH_NEON
      const UInt *interleaved;
      uint32x4_t Apply(SizeT j, uint32x4_t x) const
      {
        const auto tw = vld2q_u32(interleaved + 2 * j);
        return NeonShoupMul<F>(tw.val[0], tw.val[1], x);
      }
#endif
      UInt Apply(SizeT j, UInt x) const { return F::Mul(x, roots[j]); }
    };

    template <bool Inverse, typename F>
    inline FullWeights<F> WeightDirection(const Plan<F> &plan)
    {
      return {Inverse ? plan.inverseRoots.data() : plan.forwardRoots.data()
#if BIGMATH_NEON
              , Inverse ? plan.inverseIl.data() : plan.forwardIl.data()
#endif
      };
    }

    template <typename F>
    struct NegacyclicWeights
    {
      FactoredWeights<F> forward, inverse;
      NegacyclicWeights(UInt root, SizeT h) : forward(root, h), inverse(F::Inv(root), h) {}
    };

    template <bool Inverse, typename F>
    inline const FactoredWeights<F> &WeightDirection(const NegacyclicWeights<F> &weights)
    {
      return Inverse ? weights.inverse : weights.forward;
    }

    template <typename F, UInt G>
    inline const NegacyclicWeights<F> &GetNegacyclicWeights(Int h)
    {
      static thread_local std::unordered_map<Int, NegacyclicWeights<F>> cache;
      auto it = cache.find(h);
      if (it != cache.end()) return it->second;
      const UInt root = F::Power(G, (F::Prime - 1) / (UInt)(2 * h));
      return cache.try_emplace(h, root, (SizeT)h).first->second;
    }

    // Pack one operand for one prime into [ lower | cyc ]:
    //   lower[j] = (a_j − a_{j+H}) · ψ^j,   cyc[i] = Σ_k a_{i+kL}   (mod P).
    template <typename F, typename Weights>
    inline void TruncPackWeightFold(const std::vector<DataT> &v, BaseT base, UInt *plane,
                                    Int H, Int L, const Weights &plan2H)
    {
      UInt *lower = plane;
      UInt *cyc = plane + H;
      const SizeT cpl = CoeffsPerLimb(base);
      const SizeT coeffs = (SizeT)v.size() * cpl;
      const SizeT lowCount = std::min<SizeT>(coeffs, (SizeT)H);

      // Residues of the first min(coeffs, H) coefficients; zero padding to H.
      PackPlaneLimbs<F>(v.data(), lowCount / cpl, base, lower, H);

      // Fold those residues modulo x^L − 1.
      std::memset(cyc, 0, (SizeT)L * sizeof(UInt));
      for (SizeT k = 0; k < lowCount; k += (SizeT)L)
        AddModInto<F>(cyc, lower + k, std::min<SizeT>((SizeT)L, lowCount - k));

      // Coefficients at index ≥ H (only an operand longer than H, i.e. a
      // skewed product): x^H ≡ −1 in the lower part, and they fold into cyc
      // like any other coefficient.
      for (SizeT j = (SizeT)H; j < coeffs; ++j)
      {
        UInt raw;
        if (base == Base2_64)
          raw = (j & 1) ? (UInt)(v[j >> 1] >> 32) : (UInt)(v[j >> 1] & 0xFFFFFFFFULL);
        else
          raw = (UInt)v[j];
        const UInt r = raw % F::Prime;
        const SizeT lo = j - (SizeT)H;
        lower[lo] = F::Sub(lower[lo], r);
        const SizeT ci = j & ((SizeT)L - 1);
        cyc[ci] = F::Add(cyc[ci], r);
      }

      // Negacyclic weights ψ^j on the non-zero prefix.
      const auto &weights = WeightDirection<false>(plan2H);
      SizeT j = 0;
#if BIGMATH_NEON
      for (; j + 4 <= lowCount; j += 4)
      {
        vst1q_u32(lower + j, weights.Apply(j, vld1q_u32(lower + j)));
      }
#endif
      for (; j < lowCount; ++j)
        lower[j] = weights.Apply(j, lower[j]);
    }

    // After both inverse transforms (each scaled): unweight the lower part,
    // then recombine into [ C_lo | C_hi ] as derived above.
    template <typename F, typename Weights>
    inline void TruncUnweightReconstruct(UInt *plane, Int H, Int L, const Weights &plan2H)
    {
      UInt *lower = plane;
      UInt *cyc = plane + H;
      const auto &weights = WeightDirection<true>(plan2H);
      for (Int k = 0; k < H; k += L)
      {
        Int i = 0;
#if BIGMATH_NEON
        for (; i + 4 <= L; i += 4)
        {
          uint32x4_t x = weights.Apply(k + i, vld1q_u32(lower + k + i));
          vst1q_u32(lower + k + i, x);                              // C1[k + i]
          vst1q_u32(cyc + i, NeonSub<F>(vld1q_u32(cyc + i), x));    // C2[i] − F[i], accumulating
        }
#endif
        for (; i < L; ++i)
        {
          const UInt x = weights.Apply(k + i, lower[k + i]);
          lower[k + i] = x;
          cyc[i] = F::Sub(cyc[i], x);
        }
      }
      for (Int i = 0; i < L; ++i)
      {
        const UInt d = cyc[i];
        const UInt half = (d & 1) ? (UInt)((d + F::Prime) >> 1) : (d >> 1); // d / 2 mod P; d + P < 2^32
        cyc[i] = half;                       // C_hi[i]
        lower[i] = F::Add(lower[i], half);   // C_lo[i] = C1[i] + C_hi[i]
      }
    }

    // Pack + weight + fold + both forward transforms, one unit per
    // (operand, prime). `b` is null when squaring. plan1..3 are the plans for
    // length 2H.
    inline void CrtTruncForwardBatch(const std::vector<DataT> &a, const std::vector<DataT> *b,
                                     BaseT base, Int H, Int L,
                                     std::vector<UInt> &fa1, std::vector<UInt> &fa2, std::vector<UInt> &fa3,
                                     std::vector<UInt> &fb1, std::vector<UInt> &fb2, std::vector<UInt> &fb3,
                                     const Plan<F1> &plan1, const Plan<F2> &plan2, const Plan<F3> &plan3)
    {
      const SizeT len = (SizeT)H + (SizeT)L;
      fa1.resize(len);
      fa2.resize(len);
      fa3.resize(len);
      if (b != nullptr)
      {
        fb1.resize(len);
        fb2.resize(len);
        fb3.resize(len);
      }
      UInt *const pa1 = fa1.data(), *const pa2 = fa2.data(), *const pa3 = fa3.data();
      UInt *const pb1 = b ? fb1.data() : nullptr, *const pb2 = b ? fb2.data() : nullptr, *const pb3 = b ? fb3.data() : nullptr;
      const std::vector<DataT> *const srcA = &a;
      const Plan<F1> *const w1 = &plan1, *const h1 = &GetPlan<F1, G1>(H), *const l1 = &GetPlan<F1, G1>(L);
      const Plan<F2> *const w2 = &plan2, *const h2 = &GetPlan<F2, G2>(H), *const l2 = &GetPlan<F2, G2>(L);
      const Plan<F3> *const w3 = &plan3, *const h3 = &GetPlan<F3, G3>(H), *const l3 = &GetPlan<F3, G3>(L);
      auto unit = [=](Int idx) {
        const std::vector<DataT> &src = idx < 3 ? *srcA : *b;
        switch (idx)
        {
          case 0: case 3:
          {
            UInt *p = idx < 3 ? pa1 : pb1;
            TruncPackWeightFold<F1>(src, base, p, H, L, *w1);
            ForwardPtr<F1>(p, H, *h1);
            ForwardPtr<F1>(p + H, L, *l1);
            break;
          }
          case 1: case 4:
          {
            UInt *p = idx < 3 ? pa2 : pb2;
            TruncPackWeightFold<F2>(src, base, p, H, L, *w2);
            ForwardPtr<F2>(p, H, *h2);
            ForwardPtr<F2>(p + H, L, *l2);
            break;
          }
          case 2: case 5:
          {
            UInt *p = idx < 3 ? pa3 : pb3;
            TruncPackWeightFold<F3>(src, base, p, H, L, *w3);
            ForwardPtr<F3>(p, H, *h3);
            ForwardPtr<F3>(p + H, L, *l3);
            break;
          }
        }
      };
      const Int units = b != nullptr ? 6 : 3;
#if BIGMATH_USE_THREADS
      if (H >= (Int)BIGMATH_NTT_PARALLEL_MIN)
      {
        ParallelDo(units, [&unit](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx) unit(idx);
        });
        return;
      }
#endif
      for (Int idx = 0; idx < units; ++idx) unit(idx);
    }

    // Both inverse transforms plus the unweight/recombine pass for the three
    // (pointwise-multiplied) planes. plan1..3 are the plans for length 2H.
    inline void CrtTruncInverseBatch(std::vector<UInt> &f1, std::vector<UInt> &f2, std::vector<UInt> &f3,
                                     Int H, Int L,
                                     const Plan<F1> &plan1, const Plan<F2> &plan2, const Plan<F3> &plan3)
    {
      UInt *const a1 = f1.data(), *const a2 = f2.data(), *const a3 = f3.data();
      const Plan<F1> *const w1 = &plan1, *const h1 = &GetPlan<F1, G1>(H), *const l1 = &GetPlan<F1, G1>(L);
      const Plan<F2> *const w2 = &plan2, *const h2 = &GetPlan<F2, G2>(H), *const l2 = &GetPlan<F2, G2>(L);
      const Plan<F3> *const w3 = &plan3, *const h3 = &GetPlan<F3, G3>(H), *const l3 = &GetPlan<F3, G3>(L);
      auto whole = [=](Int prime) {
        switch (prime)
        {
          case 0:
            InversePtr<F1>(a1, H, *h1, /*scale=*/ true);
            InversePtr<F1>(a1 + H, L, *l1, /*scale=*/ true);
            TruncUnweightReconstruct<F1>(a1, H, L, *w1);
            break;
          case 1:
            InversePtr<F2>(a2, H, *h2, /*scale=*/ true);
            InversePtr<F2>(a2 + H, L, *l2, /*scale=*/ true);
            TruncUnweightReconstruct<F2>(a2, H, L, *w2);
            break;
          case 2:
            InversePtr<F3>(a3, H, *h3, /*scale=*/ true);
            InversePtr<F3>(a3 + H, L, *l3, /*scale=*/ true);
            TruncUnweightReconstruct<F3>(a3, H, L, *w3);
            break;
        }
      };
#if BIGMATH_USE_THREADS
      if (H >= (Int)BIGMATH_NTT_SPLIT_INVERSE_MIN && ParallelNumThreads() >= 4)
      {
        // Same six-unit split as CrtInverseBatch for the length-H part; the
        // small cyclic inverse rides with the second half's unit.
        const Int h = H >> 1;
        const Plan<F1> *const q1 = &GetPlan<F1, G1>(h);
        const Plan<F2> *const q2 = &GetPlan<F2, G2>(h);
        const Plan<F3> *const q3 = &GetPlan<F3, G3>(h);
        ParallelDo(6, [=](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx)
          {
            const bool second = (idx & 1) != 0;
            const Int off = second ? h : 0;
            switch (idx >> 1)
            {
              case 0:
                InversePtr<F1>(a1 + off, h, *q1, /*scale=*/ false);
                if (second) InversePtr<F1>(a1 + H, L, *l1, /*scale=*/ true);
                break;
              case 1:
                InversePtr<F2>(a2 + off, h, *q2, /*scale=*/ false);
                if (second) InversePtr<F2>(a2 + H, L, *l2, /*scale=*/ true);
                break;
              case 2:
                InversePtr<F3>(a3 + off, h, *q3, /*scale=*/ false);
                if (second) InversePtr<F3>(a3 + H, L, *l3, /*scale=*/ true);
                break;
            }
          }
        });
        const Int mid = (h >> 1) & ~(Int)3;
        ParallelDo(6, [=](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx)
          {
            const Int j0 = (idx & 1) ? mid : 0;
            const Int j1 = (idx & 1) ? h : mid;
            switch (idx >> 1)
            {
              case 0: InverseTopLayerScaled<F1>(a1, H, *h1, j0, j1); break;
              case 1: InverseTopLayerScaled<F2>(a2, H, *h2, j0, j1); break;
              case 2: InverseTopLayerScaled<F3>(a3, H, *h3, j0, j1); break;
            }
          }
        });
        ParallelDo(3, [=](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx)
          {
            switch (idx)
            {
              case 0: TruncUnweightReconstruct<F1>(a1, H, L, *w1); break;
              case 1: TruncUnweightReconstruct<F2>(a2, H, L, *w2); break;
              case 2: TruncUnweightReconstruct<F3>(a3, H, L, *w3); break;
            }
          }
        });
        return;
      }
      if (H >= (Int)BIGMATH_NTT_PARALLEL_MIN)
      {
        ParallelDo(3, [&whole](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx) whole(idx);
        });
        return;
      }
#endif
      for (Int idx = 0; idx < 3; ++idx) whole(idx);
    }

    struct PreparedOperand
    {
      BaseT base = Base2_32;
      SizeT operandLimbs = 0;
      SizeT maxOtherLimbs = 0;
      SizeT coeffsPerLimb = 1;
      ULong operandCoeffSize = 0;
      Int n = 0;
      std::vector<UInt> f1;
      std::vector<UInt> f2;
      std::vector<UInt> f3;

      bool Empty() const { return n == 0 || f1.empty(); }
    };

    inline PreparedOperand PrepareOperand(const std::vector<DataT> &operand,
                                          SizeT maxOtherLimbs,
                                          BaseT base)
    {
      PreparedOperand prepared;
      prepared.base = base;
      prepared.operandLimbs = (SizeT)operand.size();
      prepared.maxOtherLimbs = maxOtherLimbs;
      prepared.coeffsPerLimb = CoeffsPerLimb(base);

      if (IsZero(operand))
        return prepared;
      if (maxOtherLimbs == 0)
        throw std::invalid_argument("prepared NTT maxOtherLimbs must be non-zero");

      prepared.operandCoeffSize = (ULong)operand.size() * prepared.coeffsPerLimb;
      ULong maxOtherCoeffSize = (ULong)maxOtherLimbs * prepared.coeffsPerLimb;
      ULong maxCoeffCount = prepared.operandCoeffSize + maxOtherCoeffSize - 1;
      // Check the ceiling BEFORE narrowing to Int: bit_ceil > 2^30 wraps
      // negative in Int and would sail past the guard.
      ULong nWide = std::max<ULong>(2, std::bit_ceil(maxCoeffCount));
      if (nWide > (1u << 26))
        throw std::invalid_argument(
            "NttCrt: operands exceed the CRT NTT length ceiling (2^26 coefficients, "
            "~640M decimal digits) — P2/P3 have 2-adic order 2^26, so longer "
            "transforms would silently compute with invalid roots");
      prepared.n = (Int)nWide;

      prepared.f1.assign(prepared.n, 0);
      prepared.f2.assign(prepared.n, 0);
      prepared.f3.assign(prepared.n, 0);
      PackOperand(operand, base, prepared.f1, prepared.f2, prepared.f3);

      const auto &plan1 = GetPlan<F1, G1>(prepared.n);
      const auto &plan2 = GetPlan<F2, G2>(prepared.n);
      const auto &plan3 = GetPlan<F3, G3>(prepared.n);

      CrtForwardBatch(prepared.f1, prepared.f2, prepared.f3,
                      nullptr, nullptr, nullptr, plan1, plan2, plan3);

      return prepared;
    }

    inline std::vector<DataT> Multiply(const PreparedOperand &prepared,
                                       const std::vector<DataT> &other)
    {
      if (prepared.Empty() || IsZero(other))
        return std::vector<DataT>();
      if (other.size() > prepared.maxOtherLimbs)
        throw std::invalid_argument("prepared NTT operand exceeds maxOtherLimbs");

      ULong otherCoeffSize = (ULong)other.size() * prepared.coeffsPerLimb;
      ULong coeffCount = prepared.operandCoeffSize + otherCoeffSize - 1;

      static thread_local std::vector<UInt> fb1, fb2, fb3;
      fb1.assign(prepared.n, 0);
      fb2.assign(prepared.n, 0);
      fb3.assign(prepared.n, 0);
      PackOperand(other, prepared.base, fb1, fb2, fb3);

      const auto &plan1 = GetPlan<F1, G1>(prepared.n);
      const auto &plan2 = GetPlan<F2, G2>(prepared.n);
      const auto &plan3 = GetPlan<F3, G3>(prepared.n);

      CrtForwardBatch(fb1, fb2, fb3, nullptr, nullptr, nullptr, plan1, plan2, plan3);

      CrtPointwiseMulInto(fb1.data(), fb2.data(), fb3.data(),
                          prepared.f1.data(), prepared.f2.data(), prepared.f3.data(),
                          prepared.n);

      CrtInverseBatch(fb1, fb2, fb3, plan1, plan2, plan3);

      return FinalizeProduct(
          fb1, fb2, fb3, coeffCount, prepared.base, prepared.operandLimbs + other.size() + 2);
    }

#if BIGMATH_NTT_MFA
    // ─── Matrix Fourier Algorithm (Bailey 6-step) ───────────────────────────
    //
    // Recursive 2D decomposition for large NTTs that exceed L2 cache. For
    // length N = N1 * N2:
    //   1. Transpose viewing buffer as N2×N1 row-major → N1×N2 row-major
    //   2. N1 forward sub-FFTs of length N2 along (stride-1) rows
    //   3. Cross-twiddle: a[i*N2 + k2] *= ω_N^{i·k2}, i∈[0,N1), k2∈[0,N2)
    //   4. Transpose N1×N2 → N2×N1
    //   5. N2 forward sub-FFTs of length N1 along rows
    //
    // Sub-FFTs recurse via MFA when sub-size > LEAF, else hit the radix-8
    // leaf via ForwardPtr. Output is in an MFA-permuted ordering, identical
    // across both operands, so pointwise multiply is index-aligned. The
    // inverse mirrors the six steps in reverse with ω_N^{-1} twiddles and
    // returns natural-order output.
    //
    // Scaling: each leaf inverse applies 1/leaf_n. With balanced splits the
    // product of leaf scalings telescopes to 1/N — no extra top-level fixup.

    // Bit-reversal permutation table for sub-FFT length m. The MFA leaf
    // forward (DIF) writes output at bit-reversed positions, and the leaf
    // inverse (DIT) consumes bit-reversed input. The cross-twiddle must be
    // applied at the *logical* k2 index, but addressed by bit-reversed
    // storage position. Iterating k2 naturally and storing at br_table[k2]
    // keeps the exponent monotone (idx += i, never wraps within < n).
    inline const std::vector<Int> &GetBitReverseTable(Int m)
    {
      static thread_local std::unordered_map<Int, std::vector<Int>> cache;
      auto it = cache.find(m);
      if (it != cache.end()) return it->second;
      Int logm = __builtin_ctz((unsigned)m);
      std::vector<Int> tbl((SizeT)m);
      for (Int x = 0; x < m; ++x)
      {
        Int y = 0;
        for (Int b = 0; b < logm; ++b)
          y |= ((x >> b) & 1) << (logm - 1 - b);
        tbl[(SizeT)x] = y;
      }
      return cache.emplace(m, std::move(tbl)).first->second;
    }

    // Tiled out-of-place transpose, rows×cols → cols×rows.
    inline void Transpose(const UInt *src, UInt *dst, Int rows, Int cols)
    {
      constexpr Int TILE = 32;
      for (Int rb = 0; rb < rows; rb += TILE)
      {
        Int rEnd = std::min<Int>(rb + TILE, rows);
        for (Int cb = 0; cb < cols; cb += TILE)
        {
          Int cEnd = std::min<Int>(cb + TILE, cols);
          for (Int r = rb; r < rEnd; ++r)
          {
            const UInt *srcRow = src + (SizeT)r * cols;
            for (Int c = cb; c < cEnd; ++c)
              dst[(SizeT)c * rows + r] = srcRow[c];
          }
        }
      }
    }

    // Forward decl so BuildMfaPlanTree can reference it. Defined later.
    inline void MFAFactor(Int n, Int &n1, Int &n2);

    // A long MFA node only applies cross-twiddles between its two axes.
    // Only the leaves need a butterfly plan.
    struct MfaRoots
    {
      std::vector<UInt> forward;
      std::vector<UInt> inverse;
    };

#if BIGMATH_NEON
    constexpr Int MfaTwiddleBlock = 32;
    inline bool UseFactoredMfaRoots(Int n2) { return n2 >= 256; }
#endif

    template <typename F, UInt G>
    inline std::vector<UInt> BuildMfaRoots(Int n, [[maybe_unused]] Int n2, bool invert)
    {
#if BIGMATH_NEON
      if (UseFactoredMfaRoots(n2))
      {
        // For storage index b*Block + j, reversing the bits gives
        // reverse(b) + groups*reverse(j). Cache the two factors of
        // omega^(row*frequency), each with its Shoup companion. A row
        // needs just Block + groups multipliers for its n2 coefficients.
        constexpr Int Block = MfaTwiddleBlock;
        const Int groups = n2 / Block;
        const Int stride = 2 * (Block + groups);
        const Int n1 = n / n2;
        const auto &laneOrder = GetBitReverseTable(Block);
        const auto &blockOrder = GetBitReverseTable(groups);
        UInt root = F::Power(G, (F::Prime - 1) / (UInt)n);
        if (invert) root = F::Inv(root);
        std::vector<UInt> factors((SizeT)n1 * stride);
        UInt step = 1;
        for (Int r = 0; r < n1; ++r)
        {
          UInt *row = factors.data() + (SizeT)r * stride;
          const UInt laneStep = F::Power(step, (UInt)groups);
          UInt w = 1;
          for (Int j : laneOrder)
          {
            row[2 * j] = w;
            row[2 * j + 1] = ShoupCompanion<F>(w);
            w = F::Mul(w, laneStep);
          }
          w = 1;
          for (Int b : blockOrder)
          {
            row[2 * (Block + b)] = w;
            row[2 * (Block + b) + 1] = ShoupCompanion<F>(w);
            w = F::Mul(w, step);
          }
          step = F::Mul(step, root);
        }
        return factors;
      }
#endif
      return BuildRoots<F, G>(n, invert);
    }

    // Build on the calling thread. Workers read these immutable tables while
    // the caller owns their thread-local caches.
    template <typename F>
    struct MfaPlanTree
    {
      std::unordered_map<Int, const Plan<F> *> leaves;
      std::unordered_map<Int, const MfaRoots *> cross;
      const Plan<F> &Leaf(Int sz) const { return *leaves.at(sz); }
      const MfaRoots &Roots(Int sz) const { return *cross.at(sz); }
    };

    template <typename F, UInt G>
    inline void BuildMfaPlanTree(Int n, MfaPlanTree<F> &tree)
    {
      if (n <= BIGMATH_NTT_MFA_LEAF)
      {
        tree.leaves[n] = &GetPlan<F, G>(n);
        return;
      }
      if (tree.cross.count(n)) return;
      Int n1, n2;
      MFAFactor(n, n1, n2);
      tree.cross[n] = &CachedTransformTable<MfaRoots>(n, [n, n2] {
        return MfaRoots{BuildMfaRoots<F, G>(n, n2, false), BuildMfaRoots<F, G>(n, n2, true)};
      });
      BuildMfaPlanTree<F, G>(n1, tree);
      if (n2 != n1) BuildMfaPlanTree<F, G>(n2, tree);
    }

    template <typename F>
    inline void ForwardMFA(UInt *a, Int n, UInt *scratch, bool parallel,
                           const MfaPlanTree<F> &tree);

    template <typename F>
    inline void InverseMFA(UInt *a, Int n, UInt *scratch, bool parallel,
                           const MfaPlanTree<F> &tree);

    template <typename F>
    inline void ForwardMFAMul(UInt *b, Int n, UInt *scratch, bool parallel,
                              const MfaPlanTree<F> &tree, UInt *acc);

    // Wrappers for callers that have not already assembled the plan tree.
    template <typename F, UInt G>
    inline void ForwardMFA(UInt *a, Int n, UInt *scratch, bool parallel = true)
    {
      MfaPlanTree<F> tree;
      BuildMfaPlanTree<F, G>(n, tree);
      ForwardMFA<F>(a, n, scratch, parallel, tree);
    }

    template <typename F, UInt G>
    inline void ForwardMFAMul(UInt *b, Int n, UInt *scratch, UInt *acc, bool parallel = true)
    {
      MfaPlanTree<F> tree;
      BuildMfaPlanTree<F, G>(n, tree);
      ForwardMFAMul<F>(b, n, scratch, parallel, tree, acc);
    }

    template <typename F, UInt G>
    inline void InverseMFA(UInt *a, Int n, UInt *scratch, bool parallel = true)
    {
      MfaPlanTree<F> tree;
      BuildMfaPlanTree<F, G>(n, tree);
      InverseMFA<F>(a, n, scratch, parallel, tree);
    }

    // Choose N1 (slower, outer) and N2 (faster, inner) so both fit cache.
    // Convention: N1 ≥ N2, equal for log-even, N1 = 2·N2 for log-odd.
    inline void MFAFactor(Int n, Int &n1, Int &n2)
    {
      Int logN = __builtin_ctz((unsigned)n);
#ifdef BIGMATH_NTT_MFA_FORCE_LOG_N1
      Int logN1 = BIGMATH_NTT_MFA_FORCE_LOG_N1;
      if (logN1 >= logN) logN1 = logN - 1;
      Int logN2 = logN - logN1;
#else
      Int logN2 = logN >> 1;
      Int logN1 = logN - logN2;
#endif
      n1 = (Int)1 << logN1;
      n2 = (Int)1 << logN2;
    }

    inline bool UseMfaForShape(SizeT lhsLimbs, SizeT rhsLimbs, Int n)
    {
      SizeT minLimbs = std::min(lhsLimbs, rhsLimbs);
      SizeT maxLimbs = std::max(lhsLimbs, rhsLimbs);

      // Very skewed products were the shape that lost the earlier early-MFA
      // win when the global threshold was raised to protect the balanced
      // 10M-20M band. Give those shapes a lower cutover without touching the
      // balanced path.
      bool highSkew = (minLimbs > 0) && (maxLimbs / minLimbs >= 8);
      if (highSkew)
        return n >= (BIGMATH_NTT_MFA_THRESHOLD >> 1);
      return n >= BIGMATH_NTT_MFA_THRESHOLD;
    }

    // Cross-twiddle for one MFA row. After a leaf DIF sub-FFT, storage
    // position k2_br holds the value for logical k2 = br(k2_br). The correct
    // factor at storage [rowIndex, k2_br] is therefore ω_N^{±rowIndex·br(k2_br)}.
    //
    // The caller applies this immediately after the row FFT so the twiddle
    // pass is fused into step 2 / step 5 instead of costing a separate plane
    // sweep.
    template <typename F>
    inline void MfaTwiddleApplyRow(UInt *row, Int rowIndex, Int n2, Int nFull, const UInt *roots, const Int *br)
    {
      if (rowIndex == 0) return;

#if BIGMATH_NEON
      if (UseFactoredMfaRoots(n2))
      {
        constexpr Int Block = MfaTwiddleBlock;
        const Int groups = n2 / Block;
        const UInt *lanes = roots + (SizeT)rowIndex * 2 * (Block + groups);
        const UInt *blocks = lanes + 2 * Block;
        for (Int b = 0; b < groups; ++b)
        {
          const auto w = vdupq_n_u32(blocks[2 * b]);
          const auto ws = vdupq_n_u32(blocks[2 * b + 1]);
          UInt *values = row + b * Block;
          for (Int j = 0; j < Block; j += 4)
          {
            const auto lane = vld2q_u32(lanes + 2 * j);
            const auto x = NeonShoupMul<F>(lane.val[0], lane.val[1], vld1q_u32(values + j));
            vst1q_u32(values + j, NeonShoupMul<F>(w, ws, x));
          }
        }
        return;
      }
#endif
      const Int half = nFull >> 1;
      const UInt P = F::Prime;
      Int posCount = (half + rowIndex - 1) / rowIndex;
      if (posCount > n2) posCount = n2;

      Int idx = 0;
      for (Int k2 = 0; k2 < posCount; ++k2)
      {
        Int pos = br[k2];
        row[pos] = F::Mul(row[pos], roots[idx]);
        idx += rowIndex;
      }

      if (posCount < n2)
      {
        Int negIdx = idx - half;
        for (Int k2 = posCount; k2 < n2; ++k2)
        {
          Int pos = br[k2];
          row[pos] = F::Mul(row[pos], (UInt)(P - roots[negIdx]));
          negIdx += rowIndex;
        }
      }
    }

#if BIGMATH_NTT_MFA_FUSE
    // ─── Fused MFA (single level, n1,n2 <= LEAF) ─────────────────────────────
    // Collapses each (full transpose + full row-FFT) pair into one streaming
    // pass. A block of TILE output rows is gathered into an L2-resident tile,
    // FFT'd there, then written once — eliminating the scratch round-trip that
    // dominates the bandwidth-bound regime. Per-row math is identical to the
    // unfused step2/step5, so output is bit-exact.

    // Distance in UInt entries between consecutive rows of a tile holding
    // sub-FFT rows of `rowLen` entries (see BIGMATH_NTT_MFA_FUSE_TILE_PAD).
    inline SizeT MfaTilePitch(Int rowLen)
    {
      return (SizeT)rowLen + (SizeT)BIGMATH_NTT_MFA_FUSE_TILE_PAD;
    }

    // Forward stage A (step1+step2): for output rows r in [r0s,r0e), scratch row
    // r = FFT_n2(column r of a viewed n2×n1) then cross-twiddle(r). Reads `a`,
    // writes `scratch`.
    template <typename F>
    inline void FusedForwardA(const UInt *a, UInt *scratch, Int n, Int n1, Int n2,
                              const Plan<F> &planN2, const UInt *fwdRoots,
                              const Int *br, Int r0s, Int r0e)
    {
      constexpr Int TILE = BIGMATH_NTT_MFA_FUSE_TILE;
      static thread_local std::vector<UInt> tilebuf;
      const SizeT pitch = MfaTilePitch(n2);
      tilebuf.resize((SizeT)TILE * pitch);
      UInt *tile = tilebuf.data();
      for (Int r0 = r0s; r0 < r0e; r0 += TILE)
      {
        Int tr = std::min<Int>(TILE, r0e - r0);
        // Gather: tile[t*n2 + c] = a[c*n1 + (r0+t)], c in [0,n2).
        for (Int c = 0; c < n2; ++c)
        {
          const UInt *s = a + (SizeT)c * n1 + r0;
          UInt *d = tile + c;
          for (Int t = 0; t < tr; ++t) d[(SizeT)t * pitch] = s[t];
        }
        for (Int t = 0; t < tr; ++t)
        {
          UInt *row = tile + (SizeT)t * pitch;
          ForwardPtr<F>(row, n2, planN2);
          MfaTwiddleApplyRow<F>(row, r0 + t, n2, n, fwdRoots, br);
          std::copy(row, row + n2, scratch + (SizeT)(r0 + t) * n2);
        }
      }
    }

    // Forward stage B (step4+step5): for output rows rr in [r0s,r0e), a row rr =
    // FFT_n1(column rr of scratch viewed n1×n2). Reads `scratch`, writes `a`.
    template <typename F>
    inline void FusedForwardB(UInt *a, const UInt *scratch, Int n1, Int n2,
                              const Plan<F> &planN1, Int r0s, Int r0e)
    {
      constexpr Int TILE = BIGMATH_NTT_MFA_FUSE_TILE;
      static thread_local std::vector<UInt> tilebuf;
      const SizeT pitch = MfaTilePitch(n1);
      tilebuf.resize((SizeT)TILE * pitch);
      UInt *tile = tilebuf.data();
      for (Int r0 = r0s; r0 < r0e; r0 += TILE)
      {
        Int tr = std::min<Int>(TILE, r0e - r0);
        // Gather: tile[t*n1 + j] = scratch[j*n2 + (r0+t)], j in [0,n1).
        for (Int j = 0; j < n1; ++j)
        {
          const UInt *s = scratch + (SizeT)j * n2 + r0;
          UInt *d = tile + j;
          for (Int t = 0; t < tr; ++t) d[(SizeT)t * pitch] = s[t];
        }
        for (Int t = 0; t < tr; ++t)
        {
          UInt *row = tile + (SizeT)t * pitch;
          ForwardPtr<F>(row, n1, planN1);
          std::copy(row, row + n1, a + (SizeT)(r0 + t) * n1);
        }
      }
    }

    // Zero-copy view of an operand's limbs as CRT coefficients (pass fusion
    // F3, mem_pass_fusion.md). PackOperand materialized each prime's
    // coefficient plane in a separate serial sweep (plus a full-plane zero
    // fill for the padding); the fused stage A can instead pack each
    // gathered element on the fly — the limb array IS the input plane.
    // Eliminates, per operand per prime, one plane zero-fill and one packed
    // plane write+read, and moves the % P work into the parallel gather.
    struct PackedOperandView
    {
      const DataT *v = nullptr;
      SizeT limbs = 0;
      ULong coeffCount = 0; // limbs * coeffsPerLimb
      bool base64 = false;

      PackedOperandView(const std::vector<DataT> &op, BaseT base)
          : v(op.data()),
            limbs((SizeT)op.size()),
            coeffCount((ULong)op.size() * CoeffsPerLimb(base)),
            base64(base == Base2_64)
      {
      }

      // The 32-bit coefficient at idx < coeffCount as the limbs store it.
      inline UInt Raw(SizeT idx) const
      {
        if (base64)
          return (UInt)(v[idx >> 1] >> ((idx & 1) * 32));
        return (UInt)v[idx];
      }

      template <typename F>
      inline UInt Coeff(SizeT idx) const
      {
        if ((ULong)idx >= coeffCount) return 0;
        return Raw(idx) % F::Prime;
      }

      // The coefficient at idx < coeffCount before reduction: congruent to
      // Coeff(idx) and small enough to multiply by one residue or to sum
      // 2^30 of them in 64 bits.
      template <typename F>
      inline ULong Value(SizeT idx) const { return Raw(idx); }
    };

    // The same limbs with every coefficient shifted down by 2^31 into
    // [−2^31, 2^31). Products of such terms are a quarter the size of the
    // unsigned ones, which doubles the operand length the CRT modulus resolves
    // (CrtMaxBiasedTerms); FinalizeBiasedProduct adds the shift back.
    // Indices past the operand stay zero, like the unsigned view.
    struct BiasedOperandView : PackedOperandView
    {
      static constexpr ULong Shift = (ULong)1 << 31;

      explicit BiasedOperandView(const PackedOperandView &view) : PackedOperandView(view) {}

      template <typename F>
      static constexpr UInt Offset() { return F::Prime - (UInt)(Shift % F::Prime); }

      template <typename F>
      inline UInt Coeff(SizeT idx) const
      {
        if ((ULong)idx >= coeffCount) return 0;
        return F::Add(Raw(idx) % F::Prime, Offset<F>());
      }

      template <typename F>
      inline ULong Value(SizeT idx) const { return (ULong)Raw(idx) + Offset<F>(); }
    };
    static_assert(((ULong128)0xFFFFFFFFu + BiasedOperandView::Offset<F1>()) * (P1 - 1) < ((ULong128)1 << 64) &&
                  ((ULong128)0xFFFFFFFFu + BiasedOperandView::Offset<F2>()) * (P2 - 1) < ((ULong128)1 << 64) &&
                  ((ULong128)0xFFFFFFFFu + BiasedOperandView::Offset<F3>()) * (P3 - 1) < ((ULong128)1 << 64),
                  "an unreduced coefficient times a residue must fit 64 bits");

    // All coefficient readers for a prime share this worker's packing tile.
    // A gather finishes with the tile before the next gather begins.
    template <typename F>
    inline UInt *MfaPackTile(SizeT size)
    {
      static thread_local std::vector<UInt> storage;
      storage.resize(size);
      return storage.data();
    }

    // Forward stage A reading straight from the operand's limbs via a
    // PackedOperandView instead of a pre-packed coefficient plane. Identical
    // math: source element for tile slot (t, c) is coefficient c*n1+(r0+t),
    // exactly what the plane-based gather read.
    template <typename F, typename Operand>
    inline void FusedForwardAPack(const Operand &src, UInt *scratch,
                                  Int n, Int n1, Int n2,
                                  const Plan<F> &planN2, const UInt *fwdRoots,
                                  const Int *br, Int r0s, Int r0e)
    {
      constexpr Int TILE = BIGMATH_NTT_MFA_FUSE_TILE;
      const SizeT pitch = MfaTilePitch(n2);
      UInt *tile = MfaPackTile<F>((SizeT)TILE * pitch);
      for (Int r0 = r0s; r0 < r0e; r0 += TILE)
      {
        Int tr = std::min<Int>(TILE, r0e - r0);
        // Gather: tile[t*n2 + c] = coeff(c*n1 + (r0+t)), c in [0,n2).
        for (Int c = 0; c < n2; ++c)
        {
          SizeT base = (SizeT)c * n1 + r0;
          UInt *d = tile + c;
          for (Int t = 0; t < tr; ++t) d[(SizeT)t * pitch] = src.template Coeff<F>(base + t);
        }
        for (Int t = 0; t < tr; ++t)
        {
          UInt *row = tile + (SizeT)t * pitch;
          ForwardPtr<F>(row, n2, planN2);
          MfaTwiddleApplyRow<F>(row, r0 + t, n2, n, fwdRoots, br);
          std::copy(row, row + n2, scratch + (SizeT)(r0 + t) * n2);
        }
      }
    }

    // Forward stage B with the pointwise product fused into the scatter
    // (mem_pass_fusion.md F1). Row math is identical to FusedForwardB, but
    // instead of writing operand B's spectrum to its own plane, each finished
    // row is multiplied into the already-transformed operand-A plane:
    // acc[i] = acc[i] * Bhat[i]. B's spectrum never reaches DRAM and the
    // separate full-plane pointwise sweep disappears.
    // The caller must finish A's forward for these same rows first.
    template <typename F>
    inline void FusedForwardBMul(UInt *acc, const UInt *scratch, Int n1, Int n2,
                                 const Plan<F> &planN1, Int r0s, Int r0e)
    {
      constexpr Int TILE = BIGMATH_NTT_MFA_FUSE_TILE;
      static thread_local std::vector<UInt> tilebuf;
      const SizeT pitch = MfaTilePitch(n1);
      tilebuf.resize((SizeT)TILE * pitch);
      UInt *tile = tilebuf.data();
      for (Int r0 = r0s; r0 < r0e; r0 += TILE)
      {
        Int tr = std::min<Int>(TILE, r0e - r0);
        // Gather: tile[t*n1 + j] = scratch[j*n2 + (r0+t)], j in [0,n1).
        for (Int j = 0; j < n1; ++j)
        {
          const UInt *s = scratch + (SizeT)j * n2 + r0;
          UInt *d = tile + j;
          for (Int t = 0; t < tr; ++t) d[(SizeT)t * pitch] = s[t];
        }
        for (Int t = 0; t < tr; ++t)
        {
          UInt *row = tile + (SizeT)t * pitch;
          ForwardPtr<F>(row, n1, planN1);
          UInt *d = acc + (SizeT)(r0 + t) * n1;
          for (Int j = 0; j < n1; ++j) d[j] = F::Mul(d[j], row[j]);
        }
      }
    }

    // Inverse stage A (invStep5+rev4): for input rows r in [r0s,r0e), inv-FFT
    // a row r (length n1), scatter to scratch[c*n2 + r]. Reads `a`, writes
    // `scratch`.
    template <typename F>
    inline void FusedInverseA(const UInt *a, UInt *scratch, Int n1, Int n2,
                              const Plan<F> &planN1, Int r0s, Int r0e)
    {
      constexpr Int TILE = BIGMATH_NTT_MFA_FUSE_TILE;
      static thread_local std::vector<UInt> tilebuf;
      const SizeT pitch = MfaTilePitch(n1);
      tilebuf.resize((SizeT)TILE * pitch);
      UInt *tile = tilebuf.data();
      for (Int r0 = r0s; r0 < r0e; r0 += TILE)
      {
        Int tr = std::min<Int>(TILE, r0e - r0);
        for (Int t = 0; t < tr; ++t)
        {
          UInt *row = tile + (SizeT)t * pitch;
          std::copy(a + (SizeT)(r0 + t) * n1, a + (SizeT)(r0 + t) * n1 + n1, row);
          InversePtr<F>(row, n1, planN1, /*scale=*/ true);
        }
        // Scatter: scratch[c*n2 + (r0+t)] = tile[t*n1 + c], c in [0,n1).
        for (Int c = 0; c < n1; ++c)
        {
          const UInt *s = tile + c;
          UInt *d = scratch + (SizeT)c * n2 + r0;
          for (Int t = 0; t < tr; ++t) d[t] = s[(SizeT)t * pitch];
        }
      }
    }

    // Inverse stage B (invStep2+rev1): for input rows r in [r0s,r0e), inv-FFT
    // scratch row r (length n2) + cross-twiddle(r), scatter to a[c*n1 + r].
    // Reads `scratch`, writes `a`.
    template <typename F>
    inline void FusedInverseB(UInt *a, const UInt *scratch, Int n, Int n1, Int n2,
                              const Plan<F> &planN2, const UInt *invRoots,
                              const Int *br, Int r0s, Int r0e)
    {
      constexpr Int TILE = BIGMATH_NTT_MFA_FUSE_TILE;
      static thread_local std::vector<UInt> tilebuf;
      const SizeT pitch = MfaTilePitch(n2);
      tilebuf.resize((SizeT)TILE * pitch);
      UInt *tile = tilebuf.data();
      for (Int r0 = r0s; r0 < r0e; r0 += TILE)
      {
        Int tr = std::min<Int>(TILE, r0e - r0);
        for (Int t = 0; t < tr; ++t)
        {
          UInt *row = tile + (SizeT)t * pitch;
          std::copy(scratch + (SizeT)(r0 + t) * n2, scratch + (SizeT)(r0 + t) * n2 + n2, row);
          // Inverse cross-twiddle must precede the inverse row FFT (adjoint
          // of the forward fusion order; the bit-reversed addressing refers
          // to the pre-DIT layout). Same daaf3f4 ordering bug as the
          // non-tiled inverse step.
          MfaTwiddleApplyRow<F>(row, r0 + t, n2, n, invRoots, br);
          InversePtr<F>(row, n2, planN2, /*scale=*/ true);
        }
        // Scatter: a[c*n1 + (r0+t)] = tile[t*n2 + c], c in [0,n2).
        for (Int c = 0; c < n2; ++c)
        {
          const UInt *s = tile + c;
          UInt *d = a + (SizeT)c * n1 + r0;
          for (Int t = 0; t < tr; ++t) d[t] = s[(SizeT)t * pitch];
        }
      }
    }
#endif // BIGMATH_NTT_MFA_FUSE

    template <typename F>
    inline void ForwardMFA(UInt *a, Int n, UInt *scratch, bool parallel,
                           const MfaPlanTree<F> &tree)
    {
      if (n <= BIGMATH_NTT_MFA_LEAF)
      {
        ForwardPtr<F>(a, n, tree.Leaf(n));
        return;
      }

      Int n1, n2;
      MFAFactor(n, n1, n2);

#if BIGMATH_NTT_MFA_FUSE
      if (n1 <= BIGMATH_NTT_MFA_LEAF && n2 <= BIGMATH_NTT_MFA_LEAF)
      {
        const MfaRoots &roots = tree.Roots(n);
        const Plan<F> &planN2_f = tree.Leaf(n2);
        const Plan<F> &planN1_f = tree.Leaf(n1);
        const UInt *fwdRoots_f = roots.forward.data();
        const Int *br_f = GetBitReverseTable(n2).data();
        auto fa = [a, scratch, n, n1, n2, &planN2_f, fwdRoots_f, br_f](Int s, Int e) {
          FusedForwardA<F>(a, scratch, n, n1, n2, planN2_f, fwdRoots_f, br_f, s, e);
        };
        auto fb = [a, scratch, n1, n2, &planN1_f](Int s, Int e) {
          FusedForwardB<F>(a, scratch, n1, n2, planN1_f, s, e);
        };
#if BIGMATH_USE_THREADS
        if (parallel && (SizeT)n1 >= ParallelMinSize()) ParallelFor(n1, fa); else fa(0, n1);
        if (parallel && (SizeT)n2 >= ParallelMinSize()) ParallelFor(n2, fb); else fb(0, n2);
#else
        fa(0, n1);
        fb(0, n2);
#endif
        return;
      }
#endif

      // Step 1: a (n2×n1) → scratch (n1×n2)
      Transpose(a, scratch, n2, n1);

      // Step 2: n1 forward sub-FFTs of length n2 on rows of scratch.
      // Each sub-call uses the matching slice of `a` as its own scratch.
      const MfaRoots &roots = tree.Roots(n);
      const Plan<F> *planN2 = n2 <= BIGMATH_NTT_MFA_LEAF ? &tree.Leaf(n2) : nullptr;
      const std::vector<Int> &brTable = GetBitReverseTable(n2);
      const Int *br = brTable.data();
      auto step2 = [scratch, a, n2, parallel, &roots, planN2, &tree, br, n](Int rStart, Int rEnd) {
        for (Int r = rStart; r < rEnd; ++r)
        {
          UInt *row = scratch + (SizeT)r * n2;
          if (n2 <= BIGMATH_NTT_MFA_LEAF)
            ForwardPtr<F>(row, n2, *planN2);
          else
            ForwardMFA<F>(row, n2, a + (SizeT)r * n2, parallel, tree);
          MfaTwiddleApplyRow<F>(row, r, n2, n, roots.forward.data(), br);
        }
      };
#if BIGMATH_USE_THREADS
      if (parallel && (SizeT)n1 >= ParallelMinSize()) ParallelFor(n1, step2);
      else step2(0, n1);
#else
      step2(0, n1);
#endif

      // Step 4: scratch (n1×n2) → a (n2×n1)
      Transpose(scratch, a, n1, n2);

      // Step 5: n2 forward sub-FFTs of length n1 on rows of a.
      const Plan<F> *planN1 = n1 <= BIGMATH_NTT_MFA_LEAF ? &tree.Leaf(n1) : nullptr;
      auto step5 = [a, scratch, n1, parallel, planN1, &tree](Int rStart, Int rEnd) {
        for (Int r = rStart; r < rEnd; ++r)
        {
          UInt *row = a + (SizeT)r * n1;
          if (n1 <= BIGMATH_NTT_MFA_LEAF)
            ForwardPtr<F>(row, n1, *planN1);
          else
            ForwardMFA<F>(row, n1, scratch + (SizeT)r * n1, parallel, tree);
        }
      };
#if BIGMATH_USE_THREADS
      if (parallel && (SizeT)n2 >= ParallelMinSize()) ParallelFor(n2, step5);
      else step5(0, n2);
#else
      step5(0, n2);
#endif
    }

    // Forward-transform `b` and multiply its spectrum into `acc` (operand A's
    // already-transformed plane) in one pass — mem_pass_fusion.md F1. When the
    // single-level fused path doesn't apply (n at/below the leaf, a recursive
    // split, or FUSE disabled), falls back to the plain forward followed by a
    // pointwise loop; the MFA gate keeps every in-practice length (2^24..2^26)
    // on the fused branch.
    template <typename F>
    inline void ForwardMFAMul(UInt *b, Int n, UInt *scratch, bool parallel,
                              const MfaPlanTree<F> &tree, UInt *acc)
    {
#if BIGMATH_NTT_MFA_FUSE
      if (n > BIGMATH_NTT_MFA_LEAF)
      {
        Int n1, n2;
        MFAFactor(n, n1, n2);
        if (n1 <= BIGMATH_NTT_MFA_LEAF && n2 <= BIGMATH_NTT_MFA_LEAF)
        {
          const MfaRoots &roots = tree.Roots(n);
          const Plan<F> &planN2_f = tree.Leaf(n2);
          const Plan<F> &planN1_f = tree.Leaf(n1);
          const UInt *fwdRoots_f = roots.forward.data();
          const Int *br_f = GetBitReverseTable(n2).data();
          auto fa = [b, scratch, n, n1, n2, &planN2_f, fwdRoots_f, br_f](Int s, Int e) {
            FusedForwardA<F>(b, scratch, n, n1, n2, planN2_f, fwdRoots_f, br_f, s, e);
          };
          auto fb = [acc, scratch, n1, n2, &planN1_f](Int s, Int e) {
            FusedForwardBMul<F>(acc, scratch, n1, n2, planN1_f, s, e);
          };
#if BIGMATH_USE_THREADS
          if (parallel && (SizeT)n1 >= ParallelMinSize()) ParallelFor(n1, fa); else fa(0, n1);
          if (parallel && (SizeT)n2 >= ParallelMinSize()) ParallelFor(n2, fb); else fb(0, n2);
#else
          fa(0, n1);
          fb(0, n2);
#endif
          return;
        }
      }
#endif
      ForwardMFA<F>(b, n, scratch, parallel, tree);
      for (Int i = 0; i < n; ++i) acc[i] = F::Mul(acc[i], b[i]);
    }

    template <typename F>
    inline void InverseMFA(UInt *a, Int n, UInt *scratch, bool parallel,
                           const MfaPlanTree<F> &tree)
    {
      if (n <= BIGMATH_NTT_MFA_LEAF)
      {
        InversePtr<F>(a, n, tree.Leaf(n), /*scale=*/ true);
        return;
      }

      Int n1, n2;
      MFAFactor(n, n1, n2);

#if BIGMATH_NTT_MFA_FUSE
      if (n1 <= BIGMATH_NTT_MFA_LEAF && n2 <= BIGMATH_NTT_MFA_LEAF)
      {
        const MfaRoots &roots = tree.Roots(n);
        const Plan<F> &planN2_f = tree.Leaf(n2);
        const Plan<F> &planN1_f = tree.Leaf(n1);
        const UInt *invRoots_f = roots.inverse.data();
        const Int *br_f = GetBitReverseTable(n2).data();
        auto ia = [a, scratch, n1, n2, &planN1_f](Int s, Int e) {
          FusedInverseA<F>(a, scratch, n1, n2, planN1_f, s, e);
        };
        auto ib = [a, scratch, n, n1, n2, &planN2_f, invRoots_f, br_f](Int s, Int e) {
          FusedInverseB<F>(a, scratch, n, n1, n2, planN2_f, invRoots_f, br_f, s, e);
        };
#if BIGMATH_USE_THREADS
        if (parallel && (SizeT)n2 >= ParallelMinSize()) ParallelFor(n2, ia); else ia(0, n2);
        if (parallel && (SizeT)n1 >= ParallelMinSize()) ParallelFor(n1, ib); else ib(0, n1);
#else
        ia(0, n2);
        ib(0, n1);
#endif
        return;
      }
#endif

      // Reverse step 5: n2 inverse sub-FFTs of length n1 on rows of a.
      const Plan<F> *planN1 = n1 <= BIGMATH_NTT_MFA_LEAF ? &tree.Leaf(n1) : nullptr;
      auto invStep5 = [a, scratch, n1, parallel, planN1, &tree](Int rStart, Int rEnd) {
        for (Int r = rStart; r < rEnd; ++r)
        {
          UInt *row = a + (SizeT)r * n1;
          if (n1 <= BIGMATH_NTT_MFA_LEAF)
            InversePtr<F>(row, n1, *planN1, /*scale=*/ true);
          else
            InverseMFA<F>(row, n1, scratch + (SizeT)r * n1, parallel, tree);
        }
      };
#if BIGMATH_USE_THREADS
      if (parallel && (SizeT)n2 >= ParallelMinSize()) ParallelFor(n2, invStep5);
      else invStep5(0, n2);
#else
      invStep5(0, n2);
#endif

      // Reverse step 4: a (n2×n1) → scratch (n1×n2)
      Transpose(a, scratch, n2, n1);

      // Reverse step 2: n1 inverse sub-FFTs of length n2 on rows of scratch.
      const MfaRoots &roots = tree.Roots(n);
      const Plan<F> *planN2 = n2 <= BIGMATH_NTT_MFA_LEAF ? &tree.Leaf(n2) : nullptr;
      const std::vector<Int> &brTable = GetBitReverseTable(n2);
      const Int *br = brTable.data();
      auto invStep2 = [scratch, a, n2, parallel, &roots, planN2, &tree, br, n](Int rStart, Int rEnd) {
        for (Int r = rStart; r < rEnd; ++r)
        {
          UInt *row = scratch + (SizeT)r * n2;
          // The inverse cross-twiddle must PRECEDE the inverse row FFT: it is
          // the adjoint of the forward order (row FFT, then twiddle), and
          // MfaTwiddleApplyRow's bit-reversed addressing refers to the
          // pre-DIT layout. The fusion commit (daaf3f4, PR #72) applied it
          // after the row inverse, breaking forward·inverse identity at every
          // MFA size — latent because the 2^24 gate sits above all routine
          // correctness tests.
          MfaTwiddleApplyRow<F>(row, r, n2, n, roots.inverse.data(), br);
          if (n2 <= BIGMATH_NTT_MFA_LEAF)
            InversePtr<F>(row, n2, *planN2, /*scale=*/ true);
          else
            InverseMFA<F>(row, n2, a + (SizeT)r * n2, parallel, tree);
        }
      };
#if BIGMATH_USE_THREADS
      if (parallel && (SizeT)n1 >= ParallelMinSize()) ParallelFor(n1, invStep2);
      else invStep2(0, n1);
#else
      invStep2(0, n1);
#endif

      // Reverse step 1: scratch (n1×n2) → a (n2×n1)
      Transpose(scratch, a, n1, n2);
    }
#endif // BIGMATH_NTT_MFA

#if BIGMATH_NTT_MFA
    // True when n is in the single-level fused MFA window (n1, n2 <= LEAF —
    // with LEAF = 2^13 this covers every admissible CRT length up to 2^26).
    inline bool MfaFusedWindow(Int n)
    {
#if BIGMATH_NTT_MFA_FUSE
      if (n <= BIGMATH_NTT_MFA_LEAF) return false;
      Int n1, n2;
      MFAFactor(n, n1, n2);
      return n1 <= BIGMATH_NTT_MFA_LEAF && n2 <= BIGMATH_NTT_MFA_LEAF;
#else
      (void)n;
      return false;
#endif
    }

#if BIGMATH_NTT_MFA_FUSE
    // Split the tiles of all planes across the pool. At each barrier a worker
    // can own part of a plane, so the number of primes or operands does
    // not cap parallelism. Adjacent ranges never share a transpose tile.
    template <typename Body>
    inline void MfaParallelRows(Int planes, Int rows, Body body)
    {
      constexpr Int Tile = BIGMATH_NTT_MFA_FUSE_TILE;
      const Int tiles = (rows + Tile - 1) / Tile;
      ParallelDo(planes * tiles, [&](Int start, Int end) {
        while (start < end)
        {
          const Int plane = start / tiles;
          const Int first = plane * tiles;
          const Int stop = std::min(end, first + tiles);
          body(plane, (start - first) * Tile, std::min(rows, (stop - first) * Tile));
          start = stop;
        }
      });
    }

    // At large lengths a single prime supplies enough row tiles to occupy
    // the pool. Reuse the two transpose planes across all three primes.
    template <typename F, typename Operand>
    inline void MfaFusedProductPrime(const Operand &a, const Operand &b,
                                     UInt *out, UInt *scratchA, UInt *scratchB, Int n,
                                     const MfaPlanTree<F> &tree)
    {
      Int n1, n2;
      MFAFactor(n, n1, n2);
      const auto &roots = tree.Roots(n);
      const auto &plan1 = tree.Leaf(n1);
      const auto &plan2 = tree.Leaf(n2);
      const Int *br = GetBitReverseTable(n2).data();

      MfaParallelRows(2, n1, [&](Int operand, Int start, Int end) {
        FusedForwardAPack<F>(operand == 0 ? a : b, operand == 0 ? scratchA : scratchB,
                             n, n1, n2, plan2, roots.forward.data(), br, start, end);
      });
      MfaParallelRows(1, n2, [&](Int, Int start, Int end) {
        FusedForwardB<F>(out, scratchA, n1, n2, plan1, start, end);
        FusedForwardBMul<F>(out, scratchB, n1, n2, plan1, start, end);
      });
      MfaParallelRows(1, n2, [&](Int, Int start, Int end) {
        FusedInverseA<F>(out, scratchA, n1, n2, plan1, start, end);
      });
      MfaParallelRows(1, n1, [&](Int, Int start, Int end) {
        FusedInverseB<F>(out, scratchA, n, n1, n2, plan2, roots.inverse.data(), br, start, end);
      });
    }

    // Readers for length-H products modulo x^H - w^H and length-L products
    // modulo x^L - 1. Both are ordinary cyclic convolutions after packing.
    // An operand longer than `length` wraps: with weights w^i the wrapped
    // terms pick up w^H per turn, which is the reduction modulo x^H - w^H.
    template <size_t Period, typename Operand>
    struct TwistedOperandView
    {
      const Operand &src;
      const std::array<UInt, Period> &weights;
      SizeT length;

      template <typename F>
      UInt Coeff(SizeT i) const
      {
        if ((ULong)i >= src.coeffCount) return 0;
        // The unreduced coefficient times a weight fits 64 bits: one
        // reduction per term.
        UInt sum = (UInt)(src.template Value<F>(i) * weights[i % Period] % F::Prime);
        for (ULong u = (ULong)i + length; u < src.coeffCount; u += length)
          sum = F::Add(sum, (UInt)(src.template Value<F>((SizeT)u) * weights[u % Period] % F::Prime));
        return sum;
      }
    };

    // An operand longer than `length` reduced modulo x^length - 1 as it is
    // read. Each turn is another stream through the limbs, so this suits the
    // few turns of length = H; shorter cyclic lengths use FoldOperand.
    template <typename Operand>
    struct FoldedOperandView
    {
      const Operand &src;
      SizeT length;

      template <typename F>
      UInt Coeff(SizeT i) const
      {
        if ((ULong)i >= src.coeffCount) return 0;
        ULong sum = src.template Value<F>(i);
        for (ULong u = (ULong)i + length; u < src.coeffCount; u += length)
          sum += src.template Value<F>((SizeT)u);
        return (UInt)(sum % F::Prime);
      }
    };

    // For L < H, the upper halves of the two scratch planes are free during
    // the cyclic convolution. Fold in contiguous passes before the transpose.
    struct ResidueView
    {
      const UInt *data;
      template <typename F> UInt Coeff(SizeT i) const { return data[i]; }
    };

    // folded[0 .. L) = the operand modulo x^L − 1. Every turn of L
    // coefficients is read front to back; the unreduced coefficients of one
    // index sum in 64 bits and are reduced once.
    template <typename F, typename Operand>
    inline void FoldOperand(const Operand &src, Int L, UInt *folded)
    {
      ParallelFor(L, [&](Int first, Int end) {
        constexpr Int Tile = 1024;
        ULong sums[Tile];
        for (Int i0 = first; i0 < end; i0 += Tile)
        {
          const Int n = std::min<Int>(Tile, end - i0);
          std::fill_n(sums, n, (ULong)0);
          for (ULong u0 = (ULong)i0; u0 < src.coeffCount; u0 += (ULong)L)
          {
            const Int count = (Int)std::min<ULong>((ULong)n, src.coeffCount - u0);
            for (Int x = 0; x < count; ++x) sums[x] += src.template Value<F>((SizeT)(u0 + x));
          }
          for (Int x = 0; x < n; ++x) folded[i0 + x] = (UInt)(sums[x] % F::Prime);
        }
      });
    }

    template <typename F, UInt G, typename Operand>
    inline void CyclicProductViews(const Operand &a, const Operand &b, Int n,
                                   UInt *out, UInt *sa, UInt *sb)
    {
      if (MfaFusedWindow(n))
      {
        MfaPlanTree<F> tree;
        BuildMfaPlanTree<F, G>(n, tree);
        MfaFusedProductPrime<F>(a, b, out, sa, sb, n, tree);
        return;
      }
      const auto &plan = GetPlan<F, G>(n);
      for (Int i = 0; i < n; ++i)
      {
        out[i] = a.template Coeff<F>(i);
        sa[i] = b.template Coeff<F>(i);
      }
      ForwardPtr<F>(out, n, plan);
      ForwardPtr<F>(sa, n, plan);
      for (Int i = 0; i < n; ++i) out[i] = F::Mul(out[i], sa[i]);
      InversePtr<F>(out, n, plan, /*scale=*/ true);
    }

    // ─── Products above one transform: several residues of length H ─────────
    //
    // A product of c = k·H + tail coefficients (0 < tail ≤ H) is recovered
    // from k weighted cyclic products of length H and one plain cyclic
    // product of length L ≥ tail, L a power of two that divides H. Write the
    // product in blocks of H coefficients, C = C_0 + x^H·C_1 + … + x^(kH)·C_k,
    // with `tail` coefficients in C_k.
    //
    //   * Weighting both inputs by w_j^i turns a cyclic product into a product
    //     modulo x^H − s_j with s_j = w_j^H. After unweighting,
    //     R_j[i] = Σ_m s_j^m · C_m[i]: the polynomial whose coefficients are
    //     the blocks at index i, evaluated at s_j. Each prime supplies w of
    //     odd order q (TwistPeriod), so the s_j are distinct q-th roots of
    //     unity other than 1 and the weights repeat with period q.
    //   * Modulo x^L − 1 every x^(mH) is 1, so B[i] = Σ_m fold_L(C_m)[i]: the
    //     polynomial of the folded blocks evaluated at 1.
    //
    // The k + 1 nodes s_0 … s_(k−1) and 1 determine the degree-k polynomial
    // of folded blocks. Its leading coefficient is fold_L(C_k), which is C_k
    // itself because the top block fits L. With C_k known, the k values
    // R_j[i] − s_j^k · C_k[i] give C_0[i] … C_(k−1)[i] through the inverse
    // Vandermonde matrix of the s_j. For k = 1 this reads A = C_lo + s·C_hi
    // and B − fold_L(A) = (1 − s)·fold_L(C_hi).

    // Weighted residues in one product. Every prime has at least this many
    // q-th roots of unity besides 1, which belongs to the cyclic tail.
    constexpr Int MaxCeilingResidues = 6;

    template <typename F>
    constexpr UInt TwistPeriod = (F::Prime - 1) % 15 == 0 ? 15 : (F::Prime - 1) % 9 == 0 ? 9 : 7;
    static_assert((P1 - 1) % TwistPeriod<F1> == 0 && TwistPeriod<F1> > (UInt)MaxCeilingResidues);
    static_assert((P2 - 1) % TwistPeriod<F2> == 0 && TwistPeriod<F2> > (UInt)MaxCeilingResidues);
    static_assert((P3 - 1) % TwistPeriod<F3> == 0 && TwistPeriod<F3> > (UInt)MaxCeilingResidues);

    // The reconstruction sums products of reduced values in 64 bits and
    // reduces after every fourth term: four products plus a reduced
    // carry-over stay below 2^64 for every prime.
    static_assert(P1 > P2 && P1 > P3 &&
                  (ULong128)4 * (P1 - 1) * (P1 - 1) + P1 < ((ULong128)1 << 64));

    using CeilingNodes = std::array<UInt, MaxCeilingResidues + 1>;

    // Coefficients, lowest degree first, of the polynomial that is 1 at
    // nodes[j] and 0 at the other nodes[0 .. count).
    template <typename F>
    inline CeilingNodes LagrangeBasis(const CeilingNodes &nodes, Int count, Int j)
    {
      CeilingNodes poly{};
      poly[0] = 1;
      UInt scale = 1;
      Int degree = 0;
      for (Int l = 0; l < count; ++l)
      {
        if (l == j) continue;
        // poly *= x − nodes[l]
        for (Int d = ++degree; d > 0; --d)
          poly[d] = F::Sub(poly[d - 1], F::Mul(poly[d], nodes[l]));
        poly[0] = F::Sub(0, F::Mul(poly[0], nodes[l]));
        scale = F::Mul(scale, F::Sub(nodes[j], nodes[l]));
      }
      scale = F::Inv(scale);
      for (UInt &c : poly) c = F::Mul(c, scale);
      return poly;
    }

    // Weights and reconstruction coefficients for k weighted residues of
    // length H. Every coefficient that meets a weighted residue has the
    // unweighting factor w_j^(−i) multiplied in, so it depends on i mod q.
    template <typename F>
    struct CeilingTables
    {
      // Indices reconstructed per pass over the 64-bit partial sums.
      static constexpr Int Tile = 512;
      using Cycle = std::array<UInt, TwistPeriod<F>>; // by i mod q
      // A cycle laid out along a tile: the coefficients for indices
      // i .. i + Tile are Run[i mod q ..].
      using Run = std::array<UInt, Tile + TwistPeriod<F>>;
      using Row = std::array<UInt, MaxCeilingResidues>;

      std::array<Cycle, MaxCeilingResidues> weights{}; // w_j^i
      // Top block C_k from the cyclic tail and the folded residues.
      UInt topTail = 0;
      std::array<Run, MaxCeilingResidues> top{};
      // C_k at an index whose C_0 is already known (see `extra` below).
      UInt knownLow = 0;
      std::array<Cycle, MaxCeilingResidues> known{};
      // Block m from residue j, block[m·k + j], and from the top block.
      std::vector<Run> block;
      Row blockTop{};
    };

    template <typename F, UInt G>
    inline CeilingTables<F> BuildCeilingTables(Int H, Int k)
    {
      constexpr UInt q = TwistPeriod<F>;
      CeilingTables<F> t;
      t.block.resize((SizeT)k * k);
      const UInt zeta = F::Power(G, (F::Prime - 1) / q);
      std::array<typename CeilingTables<F>::Cycle, MaxCeilingResidues> unweight{};
      CeilingNodes nodes{};
      for (Int j = 0; j < k; ++j)
      {
        const UInt w = F::Power(zeta, (UInt)j + 1);
        t.weights[j][0] = 1;
        for (UInt r = 1; r < q; ++r) t.weights[j][r] = F::Mul(t.weights[j][r - 1], w);
        for (UInt r = 0; r < q; ++r) unweight[j][r] = t.weights[j][(q - r) % q];
        nodes[j] = t.weights[j][(UInt)H % q];
      }
      nodes[k] = 1;

      t.topTail = LagrangeBasis<F>(nodes, k + 1, k)[k];
      for (Int j = 0; j < k; ++j)
      {
        const UInt lead = LagrangeBasis<F>(nodes, k + 1, j)[k];
        // Column j of the inverse Vandermonde matrix of s_0 … s_(k−1).
        const CeilingNodes basis = LagrangeBasis<F>(nodes, k, j);
        // With C_0 known, (R_j − C_0) / s_j has C_k as its leading coefficient.
        const UInt known = F::Mul(basis[k - 1], F::Inv(nodes[j]));
        const UInt wrap = F::Power(nodes[j], (UInt)k);
        t.knownLow = F::Sub(t.knownLow, known);
        for (Int m = 0; m < k; ++m)
          t.blockTop[m] = F::Sub(t.blockTop[m], F::Mul(basis[m], wrap));
        for (UInt r = 0; r < q; ++r) t.known[j][r] = F::Mul(known, unweight[j][r]);
        for (UInt x = 0; x < t.top[j].size(); ++x)
        {
          t.top[j][x] = F::Mul(lead, unweight[j][x % q]);
          for (Int m = 0; m < k; ++m)
            t.block[(SizeT)m * k + j][x] = F::Mul(basis[m], unweight[j][x % q]);
        }
      }
      return t;
    }

    // One prime of the product above: out receives the k·H + tail
    // coefficients; sa and sb are scratch planes of H entries. H and L are
    // powers of two, 1 ≤ k ≤ MaxCeilingResidues and tail ≤ L ≤ H. The
    // operands may be longer than H.
    template <typename F, UInt G, typename Operand>
    inline void CeilingProductPrime(const Operand &a, const Operand &b,
                                    Int H, Int k, Int L, Int tail, UInt *out, UInt *sa, UInt *sb)
    {
      constexpr UInt q = TwistPeriod<F>;
      const CeilingTables<F> t = BuildCeilingTables<F, G>(H, k);
      UInt *const top = out + (SizeT)k * (SizeT)H;

      // When tail just exceeds L/2, recover the excess from the low
      // coefficients of the ordinary product. Multiplying two prefixes of
      // 'extra' coefficients needs just 2*extra-1 convolution coefficients.
      Int extra = 0;
      if (L >= 128 && tail > L / 2 && tail - L / 2 <= L / 8)
      {
        L /= 2;
        extra = tail - L;
        auto lowA = a, lowB = b;
        lowA.coeffCount = std::min<ULong>(a.coeffCount, (ULong)extra);
        lowB.coeffCount = std::min<ULong>(b.coeffCount, (ULong)extra);
        const Int n = std::max<Int>(64, (Int)std::bit_ceil((UInt)(2 * extra - 1)));
        CyclicProductViews<F, G>(lowA, lowB, n, out, sa, sb);
        std::copy_n(out, extra, top + L);
      }

      // The cyclic tail goes first. Only its first min(tail, L) entries are
      // needed; they wait in the top block while the weighted residues fill
      // the k blocks below with both scratch planes free again.
      if (std::max(a.coeffCount, b.coeffCount) <= (ULong)L)
        CyclicProductViews<F, G>(a, b, L, out, sa, sb);
      else if (L < H)
      {
        FoldOperand<F>(a, L, sa + L);
        FoldOperand<F>(b, L, sb + L);
        CyclicProductViews<F, G>(ResidueView{sa + L}, ResidueView{sb + L}, L, out, sa, sb);
      }
      else
        CyclicProductViews<F, G>(FoldedOperandView<Operand>{a, (SizeT)L},
                                 FoldedOperandView<Operand>{b, (SizeT)L}, L, out, sa, sb);
      const Int foldedCount = std::min(tail, L);
      std::copy_n(out, foldedCount, top);

      for (Int j = 0; j < k; ++j)
      {
        const TwistedOperandView<q, Operand> ta{a, t.weights[j], (SizeT)H}, tb{b, t.weights[j], (SizeT)H};
        CyclicProductViews<F, G>(ta, tb, H, out + (SizeT)j * (SizeT)H, sa, sb);
      }

      // The passes below accumulate 64-bit sums of products over a tile of
      // indices and reduce them after every fourth term.
      constexpr Int Tile = CeilingTables<F>::Tile;
      const auto reduce = [](ULong *sums, Int n) {
        for (Int x = 0; x < n; ++x) sums[x] %= F::Prime;
      };

      // Top block. Index i reads its own residue class modulo L in the
      // blocks below and writes only top[i] and top[L + i].
      ParallelFor(foldedCount, [&](Int first, Int end) {
        ULong sums[Tile];
        for (Int i0 = first; i0 < end; i0 += Tile)
        {
          const Int n = std::min<Int>(Tile, end - i0);
          for (Int x = 0; x < n; ++x) sums[x] = (ULong)t.topTail * top[i0 + x];
          Int pending = 1;
          for (Int j = 0; j < k; ++j)
          {
            const UInt *residue = out + (SizeT)j * (SizeT)H;
            for (Int u0 = i0; u0 < H; u0 += L)
            {
              const UInt *src = residue + u0;
              const UInt *w = t.top[j].data() + (UInt)u0 % q;
              for (Int x = 0; x < n; ++x) sums[x] += (ULong)w[x] * src[x];
              if (++pending == 4)
              {
                reduce(sums, n);
                pending = 0;
              }
            }
          }
          for (Int x = 0; x < n; ++x)
          {
            const Int i = i0 + x;
            const UInt folded = (UInt)(sums[x] % F::Prime);
            if (i < extra)
            {
              // folded = C_k[i] + C_k[L + i]. The known low coefficient
              // separates them: the residues at i hold C_0[i] … C_k[i] only.
              UInt high = F::Mul(t.knownLow, top[L + i]);
              for (Int j = 0; j < k; ++j)
                high = F::Add(high, F::Mul(t.known[j][(UInt)i % q], out[(SizeT)j * (SizeT)H + i]));
              top[i] = high;
              top[L + i] = F::Sub(folded, high);
            }
            else
              top[i] = folded;
          }
        }
      });

      // Blocks below the top. A tile reads index range [i0, i0 + n) of every
      // residue before it overwrites the same range with the blocks.
      ParallelFor(H, [&](Int first, Int end) {
        ULong sums[MaxCeilingResidues][Tile];
        for (Int i0 = first; i0 < end; i0 += Tile)
        {
          const Int n = std::min<Int>(Tile, end - i0);
          const Int withTop = std::clamp(tail - i0, 0, n);
          const UInt phase = (UInt)i0 % q;
          for (Int m = 0; m < k; ++m)
          {
            const UInt scale = t.blockTop[m];
            for (Int x = 0; x < withTop; ++x) sums[m][x] = (ULong)scale * top[i0 + x];
            for (Int x = withTop; x < n; ++x) sums[m][x] = 0;
          }
          Int pending = 1;
          for (Int j = 0; j < k; ++j)
          {
            const UInt *src = out + (SizeT)j * (SizeT)H + i0;
            for (Int m = 0; m < k; ++m)
            {
              const UInt *w = t.block[(SizeT)m * k + j].data() + phase;
              ULong *sum = sums[m];
              for (Int x = 0; x < n; ++x) sum[x] += (ULong)w[x] * src[x];
            }
            if (++pending == 4)
            {
              for (Int m = 0; m < k; ++m) reduce(sums[m], n);
              pending = 0;
            }
          }
          for (Int m = 0; m < k; ++m)
          {
            UInt *dst = out + (SizeT)m * (SizeT)H + i0;
            for (Int x = 0; x < n; ++x) dst[x] = (UInt)(sums[m][x] % F::Prime);
          }
        }
      });
    }

    // Forward stage A packs both operands for all three primes. Stage B
    // transforms matching rows of a and b, multiplying b into a on scatter.
    // a/b are zero-copy limb views (pass fusion F3 — stage A packs each
    // gathered element on the fly; no pre-packed input planes, no zero
    // fills). fa = {fa1, fa2, fa3} output planes — written in full by stage
    // B before anything reads them, so the caller only needs capacity, not
    // zeroed contents. scrs = six scratch planes (same: fully written by
    // stage A before stage B reads). Caller must pre-build the plan trees
    // in the dispatching thread and guarantee MfaFusedWindow(n).
    inline void MfaFusedForwardPointwise(const PackedOperandView &a,
                                         const PackedOperandView &b,
                                         UInt *const fa[3], UInt *const scrs[6], Int n,
                                         const MfaPlanTree<F1> &tree1,
                                         const MfaPlanTree<F2> &tree2,
                                         const MfaPlanTree<F3> &tree3)
    {
      Int n1 = 0, n2 = 0;
      MFAFactor(n, n1, n2);
      const UInt *fwd1 = tree1.Roots(n).forward.data();
      const UInt *fwd2 = tree2.Roots(n).forward.data();
      const UInt *fwd3 = tree3.Roots(n).forward.data();
      const Plan<F1> *pn2_1 = &tree1.Leaf(n2), *pn1_1 = &tree1.Leaf(n1);
      const Plan<F2> *pn2_2 = &tree2.Leaf(n2), *pn1_2 = &tree2.Leaf(n1);
      const Plan<F3> *pn2_3 = &tree3.Leaf(n2), *pn1_3 = &tree3.Leaf(n1);
      const Int *br = GetBitReverseTable(n2).data();

      MfaParallelRows(6, n1, [&](Int plane, Int start, Int end) {
        switch (plane)
        {
          case 0: FusedForwardAPack<F1>(a, scrs[0], n, n1, n2, *pn2_1, fwd1, br, start, end); break;
          case 1: FusedForwardAPack<F1>(b, scrs[1], n, n1, n2, *pn2_1, fwd1, br, start, end); break;
          case 2: FusedForwardAPack<F2>(a, scrs[2], n, n1, n2, *pn2_2, fwd2, br, start, end); break;
          case 3: FusedForwardAPack<F2>(b, scrs[3], n, n1, n2, *pn2_2, fwd2, br, start, end); break;
          case 4: FusedForwardAPack<F3>(a, scrs[4], n, n1, n2, *pn2_3, fwd3, br, start, end); break;
          case 5: FusedForwardAPack<F3>(b, scrs[5], n, n1, n2, *pn2_3, fwd3, br, start, end); break;
        }
      });

      MfaParallelRows(3, n2, [&](Int prime, Int start, Int end) {
        switch (prime)
        {
          case 0:
            FusedForwardB<F1>(fa[0], scrs[0], n1, n2, *pn1_1, start, end);
            FusedForwardBMul<F1>(fa[0], scrs[1], n1, n2, *pn1_1, start, end);
            break;
          case 1:
            FusedForwardB<F2>(fa[1], scrs[2], n1, n2, *pn1_2, start, end);
            FusedForwardBMul<F2>(fa[1], scrs[3], n1, n2, *pn1_2, start, end);
            break;
          case 2:
            FusedForwardB<F3>(fa[2], scrs[4], n1, n2, *pn1_3, start, end);
            FusedForwardBMul<F3>(fa[2], scrs[5], n1, n2, *pn1_3, start, end);
            break;
        }
      });
    }

    // Stage B reads scattered results from every stage A range. Returning
    // from the first parallel dispatch provides the required barrier.
    inline void MfaFusedInverse(UInt *const bufs[3], UInt *const scrs[3], Int n,
                                const MfaPlanTree<F1> &tree1,
                                const MfaPlanTree<F2> &tree2,
                                const MfaPlanTree<F3> &tree3)
    {
      Int n1 = 0, n2 = 0;
      MFAFactor(n, n1, n2);
      const UInt *inv1 = tree1.Roots(n).inverse.data();
      const UInt *inv2 = tree2.Roots(n).inverse.data();
      const UInt *inv3 = tree3.Roots(n).inverse.data();
      const Plan<F1> *pn2_1 = &tree1.Leaf(n2), *pn1_1 = &tree1.Leaf(n1);
      const Plan<F2> *pn2_2 = &tree2.Leaf(n2), *pn1_2 = &tree2.Leaf(n1);
      const Plan<F3> *pn2_3 = &tree3.Leaf(n2), *pn1_3 = &tree3.Leaf(n1);
      const Int *br = GetBitReverseTable(n2).data();

      MfaParallelRows(3, n2, [&](Int prime, Int start, Int end) {
        switch (prime)
        {
          case 0: FusedInverseA<F1>(bufs[0], scrs[0], n1, n2, *pn1_1, start, end); break;
          case 1: FusedInverseA<F2>(bufs[1], scrs[1], n1, n2, *pn1_2, start, end); break;
          case 2: FusedInverseA<F3>(bufs[2], scrs[2], n1, n2, *pn1_3, start, end); break;
        }
      });
      MfaParallelRows(3, n1, [&](Int prime, Int start, Int end) {
        switch (prime)
        {
          case 0: FusedInverseB<F1>(bufs[0], scrs[0], n, n1, n2, *pn2_1, inv1, br, start, end); break;
          case 1: FusedInverseB<F2>(bufs[1], scrs[1], n, n1, n2, *pn2_2, inv2, br, start, end); break;
          case 2: FusedInverseB<F3>(bufs[2], scrs[2], n, n1, n2, *pn2_3, inv3, br, start, end); break;
        }
      });
    }
#endif // BIGMATH_NTT_MFA_FUSE
#endif // BIGMATH_NTT_MFA

#if BIGMATH_NTT_MFA
    // ─── Truncated product on one thread, cache-blocked ─────────────────────
    //
    // Same arithmetic as CrtTruncForwardBatch / CrtTruncInverseBatch (product
    // mod x^H + 1 and mod x^L − 1, recombined), for H above
    // BIGMATH_NTT_TRUNC_MAX. Two differences:
    //
    //   * Each transform of length ≥ mfaMin goes through the fused MFA stages
    //     instead of a whole-array radix-8 chain, whose strided layers leave
    //     the cache at these lengths. Shorter ones (a small L) keep the chain.
    //   * One prime runs start to finish — pack, forward, pointwise, inverse,
    //     recombine — before the next begins. Only one b-side plane and one
    //     MFA scratch plane are live, and b's spectrum is multiplied into a's
    //     as its last stage writes it (ForwardMFAMul): 4 planes of H + L
    //     entries and one of H, against 9 planes of 2H for the full fused
    //     transform.
    //
    // Everything runs on the calling thread; see BIGMATH_NTT_TRUNC_SERIAL_MAX
    // for why the threaded build does not use it. `b` is null when squaring.
    // pa receives [ C_lo : H | C_hi : L ] for this prime; pb needs H + L
    // entries, scratch H. plan2H supplies the weights ψ^±j.
    template <typename F, UInt G>
    inline void TruncSerialPrime(const std::vector<DataT> &a, const std::vector<DataT> *b,
                                 BaseT base, Int H, Int L, Int mfaMin,
                                 UInt *pa, UInt *pb, UInt *scratch, const NegacyclicWeights<F> &plan2H)
    {
      const bool mfaH = H >= mfaMin && H > (Int)BIGMATH_NTT_MFA_LEAF;
      const bool mfaL = L >= mfaMin && L > (Int)BIGMATH_NTT_MFA_LEAF;
      MfaPlanTree<F> tree;
      if (mfaH) BuildMfaPlanTree<F, G>(H, tree);
      if (mfaL) BuildMfaPlanTree<F, G>(L, tree);
      const Plan<F> *planH = mfaH ? nullptr : &GetPlan<F, G>(H);
      const Plan<F> *planL = mfaL ? nullptr : &GetPlan<F, G>(L);

      TruncPackWeightFold<F>(a, base, pa, H, L, plan2H);
      if (mfaH) ForwardMFA<F>(pa, H, scratch, /*parallel=*/ false, tree);
      else ForwardPtr<F>(pa, H, *planH);
      if (mfaL) ForwardMFA<F>(pa + H, L, scratch, /*parallel=*/ false, tree);
      else ForwardPtr<F>(pa + H, L, *planL);

      if (b != nullptr)
      {
        TruncPackWeightFold<F>(*b, base, pb, H, L, plan2H);
        if (mfaH)
          ForwardMFAMul<F>(pb, H, scratch, /*parallel=*/ false, tree, pa);
        else
        {
          ForwardPtr<F>(pb, H, *planH);
          for (Int i = 0; i < H; ++i) pa[i] = F::Mul(pa[i], pb[i]);
        }
        if (mfaL)
          ForwardMFAMul<F>(pb + H, L, scratch, /*parallel=*/ false, tree, pa + H);
        else
        {
          ForwardPtr<F>(pb + H, L, *planL);
          for (Int i = 0; i < L; ++i) pa[H + i] = F::Mul(pa[H + i], pb[H + i]);
        }
      }
      else
      {
        const SizeT len = (SizeT)H + (SizeT)L;
        for (SizeT i = 0; i < len; ++i) pa[i] = F::Mul(pa[i], pa[i]);
      }

      if (mfaH) InverseMFA<F>(pa, H, scratch, /*parallel=*/ false, tree);
      else InversePtr<F>(pa, H, *planH, /*scale=*/ true);
      if (mfaL) InverseMFA<F>(pa + H, L, scratch, /*parallel=*/ false, tree);
      else InversePtr<F>(pa + H, L, *planL, /*scale=*/ true);
      TruncUnweightReconstruct<F>(pa, H, L, plan2H);
    }

    // The three primes in turn. f1..f3 receive the residue planes for the
    // Garner/carry pass; fb and scratch are workspace.
    inline void CrtTruncSerial(const std::vector<DataT> &a, const std::vector<DataT> *b,
                               BaseT base, Int H, Int L, Int mfaMin,
                               std::vector<UInt> &f1, std::vector<UInt> &f2, std::vector<UInt> &f3,
                               std::vector<UInt> &fb, std::vector<UInt> &scratch)
    {
      const SizeT len = (SizeT)H + (SizeT)L;
      f1.resize(len);
      f2.resize(len);
      f3.resize(len);
      if (b != nullptr) fb.resize(len);
      scratch.resize((SizeT)H);
      TruncSerialPrime<F1, G1>(a, b, base, H, L, mfaMin, f1.data(), fb.data(), scratch.data(), GetNegacyclicWeights<F1, G1>(H));
      TruncSerialPrime<F2, G2>(a, b, base, H, L, mfaMin, f2.data(), fb.data(), scratch.data(), GetNegacyclicWeights<F2, G2>(H));
      TruncSerialPrime<F3, G3>(a, b, base, H, L, mfaMin, f3.data(), fb.data(), scratch.data(), GetNegacyclicWeights<F3, G3>(H));
    }

    // The product through CrtTruncSerial with explicit lengths. NttCrt::Multiply
    // calls CrtTruncSerial itself when TruncSerialCyclicLength says so; this
    // entry lets tests drive the same pipeline in any build and at lengths
    // small enough for a unit test (mfaMin lowers the MFA cutover).
    // Requires H a power of two ≥ 64, L a power of two with 64 ≤ L ≤ H, and a
    // product of more than H and at most H + L coefficients. Pass the same
    // vector twice to square.
    inline std::vector<DataT> MultiplyTruncSerial(const std::vector<DataT> &a,
                                                  const std::vector<DataT> &b,
                                                  BaseT base, Int H, Int L,
                                                  Int mfaMin = BIGMATH_NTT_MFA_THRESHOLD)
    {
      const ULong coeffCount = (ULong)(a.size() + b.size()) * CoeffsPerLimb(base) - 1;
      if (H < 64 || L < 64 || L > H || !std::has_single_bit((ULong)H) || !std::has_single_bit((ULong)L) ||
          H > (1 << 25) || coeffCount <= (ULong)H || coeffCount > (ULong)H + (ULong)L)
        throw std::invalid_argument("NttCrt::MultiplyTruncSerial: lengths do not fit the product");
      static thread_local std::vector<UInt> f1, f2, f3, fb, scratch;
      CrtTruncSerial(a, &a == &b ? nullptr : &b, base, H, L, mfaMin, f1, f2, f3, fb, scratch);
      return FinalizeProduct(f1, f2, f3, coeffCount, base, a.size() + b.size() + 2);
    }
#endif // BIGMATH_NTT_MFA

    // Linear and cyclic products run successively inside Newton division.
    // Both finish reconstructing into an owned limb vector before returning.
    struct ProductWorkspace
    {
      std::vector<UInt> fa1, fb1, fa2, fb2, fa3, fb3;
      std::array<std::vector<UInt>, 6> transpose;
    };

    inline ProductWorkspace &GetProductWorkspace()
    {
      static thread_local ProductWorkspace buffers;
      return buffers;
    }

#if BIGMATH_NTT_TRUNC && BIGMATH_NTT_MFA && BIGMATH_NTT_MFA_FUSE
    // ─── Finalize for biased coefficients ───────────────────────────────────
    //
    // With h = 2^31, the sum of (a_i − h)·(b_j − h) over the n index pairs
    // that meet at one product coefficient c is c − h·(Σ a_i + Σ b_j) + h^2·n.
    // The residue planes hold that sum, so c is its centered CRT value plus h
    // times the two window sums minus h^2·n.

    // Sum of the raw coefficients of one operand that meet product
    // coefficient k: indices max(0, k − other + 1) … min(k, count − 1).
    struct CoefficientWindow
    {
      const PackedOperandView &v;
      ULong other; // coefficient count of the other operand
      ULong sum = 0;

      CoefficientWindow(const PackedOperandView &view, ULong otherCount, ULong k)
          : v(view), other(otherCount)
      {
        ULong first = k + 1 > other ? k + 1 - other : 0;
        ULong end = std::min<ULong>(k + 1, v.coeffCount);
        if (v.base64 && first < end)
        {
          // Whole limbs in the middle: both halves of each, one pass.
          if (first & 1) sum += v.Raw((SizeT)first++);
          if (first < end && (end & 1)) sum += v.Raw((SizeT)--end);
          for (ULong limb = first / 2; limb < end / 2; ++limb)
            sum += (v.v[limb] & 0xFFFFFFFFULL) + (v.v[limb] >> 32);
          return;
        }
        for (ULong i = first; i < end; ++i) sum += v.Raw((SizeT)i);
      }

      // From product coefficient k to k + 1.
      inline void Advance(ULong k)
      {
        if (k + 1 < v.coeffCount) sum += v.Raw((SizeT)(k + 1));
        if (k + 1 >= other) sum -= v.Raw((SizeT)(k + 1 - other));
      }
    };

    struct StandardCrt
    {
      using Inverses = InvTable;
      static constexpr ULong128 Modulus = CrtModulus;
      static const Inverses &GetInverses() { return GarnerInverses(); }
      static ULong128 Value(UInt r1, UInt r2, UInt r3, const Inverses &inv)
      {
        return Garner(r1, r2, r3, inv);
      }
    };

    // Limbs [ls, le) of the product from the residues of the biased
    // convolution, computed with carry-in 0. Returns the carry out of limb
    // le − 1. For up to 2^30 terms, each true coefficient is below 2^94.
    // In radix B = 2^32 the carry is at most terms·(B − 1) < 2^62:
    // (terms·(B − 1)^2 + terms·(B − 1)) / B = terms·(B − 1).
    template <typename Crt = StandardCrt>
    inline ULong BiasedCarryRange(const UInt *f1, const UInt *f2, const UInt *f3,
                                  const PackedOperandView &a, const PackedOperandView &b,
                                  DataT *out, SizeT ls, SizeT le, SizeT perLimb, const typename Crt::Inverses &inv)
    {
      constexpr ULong128 half = (Crt::Modulus - 1) / 2;
      const ULong coeffCount = a.coeffCount + b.coeffCount - 1;
      ULong k = (ULong)ls * perLimb;
      CoefficientWindow wa(a, b.coeffCount, k), wb(b, a.coeffCount, k);
      // Recover product coefficient k from the centered CRT residue, then
      // move the windows to k + 1.
      auto coefficient = [&](ULong128 x) {
        const ULong first = k + 1 > b.coeffCount ? k + 1 - b.coeffCount : 0;
        const ULong terms = std::min<ULong>(k + 1, a.coeffCount) - first;
        if (x > half) x -= Crt::Modulus;
        x += ((ULong128)(wa.sum + wb.sum) << 31) - ((ULong128)terms << 62);
        wa.Advance(k);
        wb.Advance(k);
        ++k;
        return x;
      };

      ULong carry = 0;
      SizeT i = ls;
#if BIGMATH_NEON
      if constexpr (std::is_same_v<Crt, StandardCrt>)
      {
        if (perLimb == 2)
        {
          const GarnerNeonConsts consts(inv);
          const ULong p1 = P1;
          const ULong p1p2 = (ULong)inv.p1p2;
          for (; i + 2 <= le && k + 4 <= coeffCount; i += 2)
          {
            const UInt *r1 = f1 + k;
            uint32x4_t u2v, u3v;
            NeonGarnerDigits(vld1q_u32(r1), vld1q_u32(f2 + k), vld1q_u32(f3 + k), consts, u2v, u3v);
            UInt u2[4], u3[4];
            vst1q_u32(u2, u2v);
            vst1q_u32(u3, u3v);
            for (SizeT h = 0; h < 2; ++h)
            {
              const ULong128 x0 = coefficient((ULong128)p1p2 * u3[2 * h] + ((ULong)r1[2 * h] + p1 * u2[2 * h]));
              const ULong128 x1 =
                  coefficient((ULong128)p1p2 * u3[2 * h + 1] + ((ULong)r1[2 * h + 1] + p1 * u2[2 * h + 1]));
              const ULong128 t = x0 + (x1 << 32) + carry; // < 2^125
              out[i + h] = (DataT)t;
              carry = (ULong)(t >> 64);
            }
          }
        }
      }
#endif
      for (; i < le; ++i)
      {
        // Even with two coefficients, x0 + (x1 << 32) + carry < 2^127.
        ULong128 t = carry;
        for (SizeT h = 0; h < perLimb && k < coeffCount; ++h)
          t += coefficient(Crt::Value(f1[k], f2[k], f3[k], inv)) << (32 * h);
        if (perLimb == 2)
        {
          out[i] = (DataT)t;
          carry = (ULong)(t >> 64);
        }
        else
        {
          out[i] = (DataT)(t & 0xFFFFFFFFULL);
          carry = (ULong)(t >> 32);
        }
      }
      return carry;
    }

    // FinalizeProduct for residue planes computed from BiasedOperandView
    // inputs; a and b are the unbiased views of the same operands.
    template <typename Crt = StandardCrt>
    inline std::vector<DataT> FinalizeBiasedProduct(const std::vector<UInt> &fa1,
                                                    const std::vector<UInt> &fa2,
                                                    const std::vector<UInt> &fa3,
                                                    const PackedOperandView &a,
                                                    const PackedOperandView &b,
                                                    BaseT base, SizeT reserveLimbs)
    {
      const SizeT perLimb = CoeffsPerLimb(base);
      const ULong coeffCount = a.coeffCount + b.coeffCount - 1;
      const SizeT nLimbs = (SizeT)((coeffCount + perLimb - 1) / perLimb);
      const auto &inv = Crt::GetInverses();
      std::vector<DataT> result;
      result.reserve(std::max<SizeT>(reserveLimbs, nLimbs + 2));
      result.resize(nLimbs);
      DataT *out = result.data();
      auto range = [&](SizeT ls, SizeT le) {
        return BiasedCarryRange<Crt>(fa1.data(), fa2.data(), fa3.data(), a, b, out, ls, le, perLimb, inv);
      };
      // The chunk stitch ripples whole 64-bit limbs.
      ULong carry = perLimb == 2 ? CarryChunks64(out, nLimbs, range) : range((SizeT)0, nLimbs);
      const unsigned limbBits = 32 * perLimb;
      while (carry != 0)
      {
        result.push_back(limbBits == 64 ? (DataT)carry : (DataT)(carry & 0xFFFFFFFFULL));
        carry = limbBits == 64 ? 0 : carry >> 32;
      }
      TrimZeros(result);
      return result;
    }

    // The product of a and b through CeilingProductPrime with transforms no
    // longer than H: more than H and at most (MaxCeilingResidues + 1)·H
    // coefficients. `biased` shifts the coefficients (BiasedOperandView); the
    // shorter operand may then have up to CrtMaxBiasedTerms coefficients
    // instead of CrtMaxTerms.
    inline std::vector<DataT> MultiplyAtCrtCeiling(const std::vector<DataT> &a,
                                                  const std::vector<DataT> &b, BaseT base, Int H,
                                                  bool biased)
    {
      const PackedOperandView va(a, base), vb(b, base);
      const ULong coeffCount = va.coeffCount + vb.coeffCount - 1;
      const ULong blocks = (coeffCount - 1) / (ULong)H;
      if (H < 64 || !std::has_single_bit((UInt)H) || blocks < 1 || blocks > (ULong)MaxCeilingResidues)
        throw std::invalid_argument("NttCrt::MultiplyAtCrtCeiling: product does not fit the residues");
      const Int k = (Int)blocks;
      const Int tail = (Int)(coeffCount - blocks * (ULong)H);
      const Int L = std::max<Int>(64, (Int)std::bit_ceil((UInt)tail));
      auto &[fa1, fb1, fa2, fb2, fa3, fb3, mfaScratch] = GetProductWorkspace();
      fa1.resize(coeffCount); fa2.resize(coeffCount); fa3.resize(coeffCount);
      // Reuse the same two work planes as the preceding transform band.
      const bool single = ParallelNumThreads() == 1;
      auto &scratchA = single ? fb1 : mfaScratch[0];
      auto &scratchB = single ? fb2 : mfaScratch[1];
      scratchA.resize(H); scratchB.resize(H);
      UInt *sa = scratchA.data(), *sb = scratchB.data();
      UInt *const planes[3] = {fa1.data(), fa2.data(), fa3.data()};
      const auto primes = [&](const auto &xa, const auto &xb) {
        CeilingProductPrime<F1, G1>(xa, xb, H, k, L, tail, planes[0], sa, sb);
        CeilingProductPrime<F2, G2>(xa, xb, H, k, L, tail, planes[1], sa, sb);
        CeilingProductPrime<F3, G3>(xa, xb, H, k, L, tail, planes[2], sa, sb);
      };
      const SizeT reserveLimbs = (SizeT)(a.size() + b.size() + 2);
      if (biased)
      {
        primes(BiasedOperandView(va), BiasedOperandView(vb));
        return FinalizeBiasedProduct(fa1, fa2, fa3, va, vb, base, reserveLimbs);
      }
      primes(va, vb);
      return FinalizeProduct(fa1, fa2, fa3, coeffCount, base, reserveLimbs);
    }

    // Unsigned coefficients while the shorter operand allows them: their
    // finalize pass is the cheaper one.
    inline std::vector<DataT> MultiplyAtCrtCeiling(const std::vector<DataT> &a,
                                                  const std::vector<DataT> &b, BaseT base, Int H)
    {
      const ULong shorter = (ULong)std::min(a.size(), b.size()) * CoeffsPerLimb(base);
      return MultiplyAtCrtCeiling(a, b, base, H, shorter > CrtMaxTerms);
    }
#endif

    inline std::vector<DataT> Multiply(const std::vector<DataT> &a,
                                       const std::vector<DataT> &b, BaseT base);

    // Spread the tail across all blocks so neighboring block pairs can reuse
    // the same transform plans and scratch capacities.
    inline size_t ProductBlockLength(size_t limbs, size_t room)
    {
      const size_t count = 1 + (limbs - 1) / room;
      return 1 + (limbs - 1) / count;
    }

    inline const std::vector<DataT> &ProductBlock(const std::vector<DataT> &v, size_t offset,
                                                  size_t length, std::vector<DataT> &scratch)
    {
      if (length == v.size()) return v;
      const size_t end = std::min(v.size(), offset + length);
      scratch.assign(v.begin() + offset, v.begin() + end);
      return scratch;
    }

    // Multiply by overlap-add in the limb radix. The shorter operand is cut
    // into blocks of at most shorterLimbs limbs and the longer one into
    // blocks that keep the limbs of a pair within pairLimbs (at least twice
    // shorterLimbs); product(left, right) multiplies one pair. Each partial
    // result is fully carry-normalized before adding it to the preallocated
    // destination.
    template <typename Product>
    inline std::vector<DataT> MultiplyBlocked(const std::vector<DataT> &a,
                                              const std::vector<DataT> &b, BaseT base,
                                              size_t pairLimbs, size_t shorterLimbs,
                                              const Product &product)
    {
      if (a.empty() || b.empty()) return {};
      if (a.size() + b.size() > (size_t)std::numeric_limits<Int>::max())
        throw std::length_error("CRT product exceeds the supported limb index range");
      // The square below pairs equal blocks, so it cuts both sides alike.
      const bool square = &a == &b;
      const size_t shortBlock = ProductBlockLength(std::min(a.size(), b.size()), shorterLimbs);
      const size_t longBlock = square ? shortBlock
          : ProductBlockLength(std::max(a.size(), b.size()), pairLimbs - shortBlock);
      const size_t aBlock = a.size() < b.size() ? shortBlock : longBlock;
      const size_t bBlock = a.size() < b.size() ? longBlock : shortBlock;
      std::vector<DataT> result(a.size() + b.size(), 0);
      std::vector<DataT> aStorage, bStorage;
      for (size_t i = 0; i < a.size(); i += aBlock)
      {
        const auto &left = ProductBlock(a, i, aBlock, aStorage);
        for (size_t j = square ? i : 0; j < b.size(); j += bBlock)
        {
          const auto &right = square && i == j ? left : ProductBlock(b, j, bBlock, bStorage);
          const std::vector<DataT> partial = product(left, right);
          if (partial.empty()) continue;
          const SizeT offset = (SizeT)(i + j), last = (SizeT)partial.size() - 1;
          AddTo(result, offset, offset + last, partial, 0, last, base);
          if (square && i != j)
            AddTo(result, offset, offset + last, partial, 0, last, base);
        }
      }
      TrimZeros(result);
      return result;
    }

    // Block pairs that each fit one transform of transformLimbs limbs.
    inline std::vector<DataT> MultiplyBlocked(const std::vector<DataT> &a,
                                              const std::vector<DataT> &b, BaseT base,
                                              size_t transformLimbs)
    {
      return MultiplyBlocked(a, b, base, transformLimbs, transformLimbs / 2,
                             [base](const std::vector<DataT> &left, const std::vector<DataT> &right) {
                               return Multiply(left, right, base);
                             });
    }

#include "NTTWideCrt.h"

    // A product of more than H coefficients, H being the longest transform.
    //
    // While the shorter operand is longer than H/2 — a single transform could
    // not hold it against a block of equal length — the product is one set of
    // residues (MultiplyAtCrtCeiling): unsigned coefficients up to plainTerms
    // in the shorter operand, biased ones up to biasedTerms. Past that limit
    // or the standard node limit, try the wide prime set. If neither set fits,
    // cut the operands into blocks whose pairs fit the available limits.
    // A shorter operand of at most H/2 coefficients is instead paired with
    // blocks that keep each product within one transform.
    inline std::vector<DataT> MultiplyAboveTransform(const std::vector<DataT> &a,
                                                     const std::vector<DataT> &b, BaseT base, Int H,
                                                     ULong plainTerms = CrtMaxTerms,
                                                     ULong biasedTerms = CrtMaxBiasedTerms,
                                                     ULong wideTerms = BIGMATH_NTT_WIDE_MAX_TERMS,
                                                     ULong wideCoeffs = BIGMATH_NTT_WIDE_MAX_COEFFS)
    {
      const SizeT perLimb = CoeffsPerLimb(base);
      ULong pairCoeffs = (ULong)H, shorterCoeffs = (ULong)H / 2;
#if BIGMATH_NTT_TRUNC && BIGMATH_NTT_MFA && BIGMATH_NTT_MFA_FUSE
      const ULong aCoeffs = (ULong)a.size() * perLimb, bCoeffs = (ULong)b.size() * perLimb;
      const ULong shorter = std::min(aCoeffs, bCoeffs);
      const ULong residueCoeffs = (ULong)(MaxCeilingResidues + 1) * (ULong)H;
      // Transforms of length H must be the cache-blocked ones or small.
      if (shorter > (ULong)H / 2 && (H <= (Int)BIGMATH_NTT_MFA_LEAF || MfaFusedWindow(H)))
      {
        if (shorter <= biasedTerms && aCoeffs + bCoeffs - 1 <= residueCoeffs)
          return MultiplyAtCrtCeiling(a, b, base, H, shorter > plainTerms);
        pairCoeffs = residueCoeffs;
        shorterCoeffs = std::min(biasedTerms, residueCoeffs / 2);
#if BIGMATH_NTT_WIDE
        const Int wideH = std::min(H, Wide::WideResidueLength);
        const ULong wideCapacity = std::min(wideCoeffs, (ULong)(Wide::WideMaxResidues + 1) * wideH);
        const ULong wideShorter = std::min({wideTerms, Wide::WideMaxBiasedTerms, wideCapacity / 2});
        const ULong coeffCount = aCoeffs + bCoeffs - 1;
        if (shorter <= std::min(wideTerms, Wide::WideMaxBiasedTerms) && coeffCount <= wideCapacity)
          return Wide::MultiplyWide(a, b, base, wideH);
        if (wideShorter >= shorterCoeffs && wideCapacity >= pairCoeffs)
        {
          pairCoeffs = wideCapacity;
          shorterCoeffs = wideShorter;
        }
#endif
      }
#else
      (void)plainTerms;
      (void)biasedTerms;
#endif
#if !BIGMATH_NTT_WIDE || !(BIGMATH_NTT_TRUNC && BIGMATH_NTT_MFA && BIGMATH_NTT_MFA_FUSE)
      (void)wideTerms;
      (void)wideCoeffs;
#endif
      return MultiplyBlocked(a, b, base, (size_t)(pairCoeffs / perLimb), (size_t)(shorterCoeffs / perLimb),
                             [&](const std::vector<DataT> &left, const std::vector<DataT> &right) {
                               const ULong coeffs = (ULong)(left.size() + right.size()) * perLimb - 1;
                               return coeffs > (ULong)H
                                          ? MultiplyAboveTransform(left, right, base, H, plainTerms, biasedTerms,
                                                                   wideTerms, wideCoeffs)
                                          : Multiply(left, right, base);
                             });
    }

    // ─── Public Multiply ─────────────────────────────────────────────────────

    inline std::vector<DataT> Multiply(const std::vector<DataT> &a,
                                       const std::vector<DataT> &b,
                                       BaseT base)
    {
      if (IsZero(a) || IsZero(b)) return std::vector<DataT>();
      if (a.size() == 1) return ClassicMultiplication::Multiply(b, a[0], base);
      if (b.size() == 1) return ClassicMultiplication::Multiply(a, b[0], base);

      // Squaring: when both operands are the same vector, the b-side pack and
      // forward transforms are byte-identical to the a-side — skip them and
      // square fa pointwise. (The fused MFA window below packs straight from
      // the limb arrays and fuses the pointwise into its stages, so it keeps
      // the two-operand flow.)
      const bool sameOperand = (&a == &b);

      // Split into 32-bit coefficients. Works for Base2_32 (1 coeff per limb)
      // and Base2_64 (2 coeffs per limb). Pre-CRT bounds: each coefficient is
      // < 2^32, convolution sum at length N is < N · 2^64.
      SizeT coeffsPerLimb = CoeffsPerLimb(base);

      ULong aCoeffSize = (ULong)a.size() * coeffsPerLimb;
      ULong bCoeffSize = (ULong)b.size() * coeffsPerLimb;
      ULong coeffCount = aCoeffSize + bCoeffSize - 1;
      // Every transform must fit the roots supported by all three primes.
      constexpr ULong MaxTransform = (ULong)1 << 26;
      if (coeffCount > MaxTransform)
        return MultiplyAboveTransform(a, b, base, (Int)MaxTransform);
      const Int n = (Int)std::bit_ceil(coeffCount);

      auto &[fa1, fb1, fa2, fb2, fa3, fb3, mfaScratch] = GetProductWorkspace();
      auto multiplySpectra = [&](Int length) {
        CrtPointwiseMulInto(fa1.data(), fa2.data(), fa3.data(),
                            sameOperand ? fa1.data() : fb1.data(),
                            sameOperand ? fa2.data() : fb2.data(),
                            sameOperand ? fa3.data() : fb3.data(), length);
      };

      // Just above a power of two: H + L instead of 2H (see TruncCyclicLength).
      Int truncL = TruncCyclicLength(coeffCount, n >> 1);
#if BIGMATH_NTT_MFA && BIGMATH_NTT_MFA_FUSE
      // At this transform size, skewed operands use less memory and time
      // when the threaded matrix transform packs directly into its tiles.
      if (truncL && n >= (1 << 21) && ParallelNumThreads() > 1 &&
          std::max(a.size(), b.size()) / std::min(a.size(), b.size()) >= 8 &&
          UseMfaForShape(a.size(), b.size(), n) && MfaFusedWindow(n))
        truncL = 0;
#endif
      if (truncL)
      {
        const auto &plan1 = GetPlan<F1, G1>(n);
        const auto &plan2 = GetPlan<F2, G2>(n);
        const auto &plan3 = GetPlan<F3, G3>(n);
        const Int H = n >> 1;
        CrtTruncForwardBatch(a, sameOperand ? nullptr : &b, base, H, truncL,
                             fa1, fa2, fa3, fb1, fb2, fb3, plan1, plan2, plan3);
        multiplySpectra(H + truncL);
        CrtTruncInverseBatch(fa1, fa2, fa3, H, truncL, plan1, plan2, plan3);
        return FinalizeProduct(fa1, fa2, fa3, coeffCount, base, a.size() + b.size() + 2);
      }

#if BIGMATH_NTT_MFA
      // The same scheme above BIGMATH_NTT_TRUNC_MAX, single-threaded only —
      // see CrtTruncSerial. fb1 is the one b-side plane, fb2 the MFA scratch.
      if (const Int serialL = TruncSerialCyclicLength(coeffCount, n >> 1); serialL != 0)
      {
        CrtTruncSerial(a, sameOperand ? nullptr : &b, base, n >> 1, serialL,
                       (Int)BIGMATH_NTT_MFA_THRESHOLD, fa1, fa2, fa3, fb1, fb2);
        return FinalizeProduct(fa1, fa2, fa3, coeffCount, base, a.size() + b.size() + 2);
      }

      if (UseMfaForShape(a.size(), b.size(), n))
      {
        MfaPlanTree<F1> tree1;
        MfaPlanTree<F2> tree2;
        MfaPlanTree<F3> tree3;
        // Acquire plans before dispatch so every row tile can use them directly.
        BuildMfaPlanTree<F1, G1>(n, tree1);
        BuildMfaPlanTree<F2, G2>(n, tree2);
        BuildMfaPlanTree<F3, G3>(n, tree3);

#if BIGMATH_NTT_MFA_FUSE
        if (MfaFusedWindow(n))
        {
          PackedOperandView va(a, base), vb(b, base);
          // Conversion also retains scratch for its worker subtrees and divider chain.
          const bool reuseScratch = n >= (1 << 24) ||
              (n >= (1 << 20) && detail::DecimalConversionScope::Active());
          if (reuseScratch)
          {
            // Earlier small products on this thread may have used six planes.
            for (Int i = 2; i < 6; ++i) std::vector<UInt>().swap(mfaScratch[i]);
          }
          for (Int i = 0; i < (reuseScratch ? 2 : 6); ++i) mfaScratch[i].resize(n);
          fa1.resize(n); fa2.resize(n); fa3.resize(n);
          if (reuseScratch)
          {
            UInt *sa = mfaScratch[0].data(), *sb = mfaScratch[1].data();
            MfaFusedProductPrime<F1>(va, vb, fa1.data(), sa, sb, n, tree1);
            MfaFusedProductPrime<F2>(va, vb, fa2.data(), sa, sb, n, tree2);
            MfaFusedProductPrime<F3>(va, vb, fa3.data(), sa, sb, n, tree3);
          }
          else
          {
            UInt *fap[3] = {fa1.data(), fa2.data(), fa3.data()};
            UInt *scrs[6] = {mfaScratch[0].data(), mfaScratch[1].data(), mfaScratch[2].data(),
                             mfaScratch[3].data(), mfaScratch[4].data(), mfaScratch[5].data()};
            MfaFusedForwardPointwise(va, vb, fap, scrs, n, tree1, tree2, tree3);
            MfaFusedInverse(fap, scrs, n, tree1, tree2, tree3);
          }
        }
        else
#endif
        {
          // F4 audit note (mem_pass_fusion.md): with LEAF = 2^13 the fused
          // window covers every n ≤ 2^26 — the entire admissible CRT range
          // under the ceiling guard above. This multi-level branch is
          // therefore reachable ONLY in -DBIGMATH_NTT_MFA_FUSE=0 builds,
          // where it (and ForwardMFA/InverseMFA's standalone Transpose
          // sweeps) is the intended fallback. Kept, not dead code.
          for (int i = 0; i < 6; ++i) mfaScratch[i].assign(n, 0);
          fa1.assign(n, 0);
          fa2.assign(n, 0);
          fa3.assign(n, 0);
          PackOperand(a, base, fa1, fa2, fa3);
          if (!sameOperand)
          {
            fb1.assign(n, 0); fb2.assign(n, 0); fb3.assign(n, 0);
            PackOperand(b, base, fb1, fb2, fb3);
          }
          UInt *bufs[6]   = {fa1.data(), fb1.data(), fa2.data(), fb2.data(), fa3.data(), fb3.data()};
          UInt *scrs[6]   = {mfaScratch[0].data(), mfaScratch[1].data(), mfaScratch[2].data(),
                             mfaScratch[3].data(), mfaScratch[4].data(), mfaScratch[5].data()};
          if (sameOperand)
          {
            // Squaring: only the three a-side forwards.
            auto fwdBody = [bufs, scrs, n, &tree1, &tree2, &tree3](Int s, Int e) {
              for (Int idx = s; idx < e; ++idx)
              {
                switch (idx)
                {
                  case 0: ForwardMFA<F1>(bufs[0], n, scrs[0], /*parallel=*/false, tree1); break;
                  case 1: ForwardMFA<F2>(bufs[2], n, scrs[1], /*parallel=*/false, tree2); break;
                  case 2: ForwardMFA<F3>(bufs[4], n, scrs[2], /*parallel=*/false, tree3); break;
                }
              }
            };
            ParallelDo(3, fwdBody);
          }
          else
          {
            // Multi-level MFA: original 6-unit forward batch; the pointwise
            // sweep below handles the product.
            auto fwdBody = [bufs, scrs, n, &tree1, &tree2, &tree3](Int s, Int e) {
              for (Int idx = s; idx < e; ++idx)
              {
                switch (idx)
                {
                  case 0: ForwardMFA<F1>(bufs[0], n, scrs[0], /*parallel=*/false, tree1); break;
                  case 1: ForwardMFA<F1>(bufs[1], n, scrs[1], /*parallel=*/false, tree1); break;
                  case 2: ForwardMFA<F2>(bufs[2], n, scrs[2], /*parallel=*/false, tree2); break;
                  case 3: ForwardMFA<F2>(bufs[3], n, scrs[3], /*parallel=*/false, tree2); break;
                  case 4: ForwardMFA<F3>(bufs[4], n, scrs[4], /*parallel=*/false, tree3); break;
                  case 5: ForwardMFA<F3>(bufs[5], n, scrs[5], /*parallel=*/false, tree3); break;
                }
              }
            };
            ParallelDo(6, fwdBody);
          }
          multiplySpectra(n);
          auto invBody = [bufs, scrs, n, &tree1, &tree2, &tree3](Int s, Int e) {
            for (Int idx = s; idx < e; ++idx)
            {
              switch (idx)
              {
                case 0: InverseMFA<F1>(bufs[0], n, scrs[0], /*parallel=*/false, tree1); break;
                case 1: InverseMFA<F2>(bufs[2], n, scrs[1], /*parallel=*/false, tree2); break;
                case 2: InverseMFA<F3>(bufs[4], n, scrs[2], /*parallel=*/false, tree3); break;
              }
            }
          };
          ParallelDo(3, invBody);
        }
        return FinalizeProduct(fa1, fa2, fa3, coeffCount, base, a.size() + b.size() + 2);
      }
#endif
      const auto &plan1 = GetPlan<F1, G1>(n);
      const auto &plan2 = GetPlan<F2, G2>(n);
      const auto &plan3 = GetPlan<F3, G3>(n);
      CrtPackForwardBatch(a, sameOperand ? nullptr : &b, base, n,
                          fa1, fa2, fa3, fb1, fb2, fb3, plan1, plan2, plan3);
      multiplySpectra(n);
      CrtInverseBatch(fa1, fa2, fa3, plan1, plan2, plan3);
      return FinalizeProduct(fa1, fa2, fa3, coeffCount, base, a.size() + b.size() + 2);
    }

    // ─── Cyclic (wrap-around) product: (a · b) mod (B^L − 1) ────────────────
    //
    // An NTT of length N computes coefficient products mod (x^N − 1); with
    // N = L · coeffsPerLimb the wrapped coefficient sums carry-propagate to
    // exactly (a · b) mod (B^L − 1), because B^L ≡ 1 there. The transform is
    // HALF the length the full product would need — and transform length is
    // the serial critical path of a multiply (the three primes already run in
    // parallel), so callers that can reconstruct what they need from the
    // residue get the product at about half the wall-clock cost.
    //
    // Caller contract: L · coeffsPerLimb is a power of two (so the cyclic
    // length is an admissible NTT size), a.size() ≤ L, b.size() ≤ L, and
    // N ≤ 2^26 (the CRT length ceiling: coefficient-sum headroom N·2^64 <
    // p1·p2·p3 ≈ 2^90.5 and the primes' 2-adic order). The cap was 2^22
    // while this path predated MFA support; the MFA permutation cancels
    // between forward and pointwise (both operands share it) and the
    // inverse returns natural order, so the fused MFA pipeline is valid
    // here and engages at n ≥ BIGMATH_NTT_MFA_THRESHOLD (2026-06-12).
    // Returns the canonical residue in [0, B^L − 1), trimmed.
    inline std::vector<DataT> MultiplyMod2km1(const std::vector<DataT> &a,
                                              const std::vector<DataT> &b,
                                              SizeT L,
                                              BaseT base)
    {
      if (IsZero(a) || IsZero(b))
        return std::vector<DataT>{0};

      SizeT coeffsPerLimb = CoeffsPerLimb(base);
      Int n = (Int)((ULong)L * coeffsPerLimb);
      if (n <= 1 || (n & (n - 1)) != 0 || n > (1 << 26))
        throw std::invalid_argument("MultiplyMod2km1: L*coeffsPerLimb must be a power of two <= 2^26");
      if (a.size() > L || b.size() > L)
        throw std::invalid_argument("MultiplyMod2km1: operand exceeds L limbs");

      auto &[fa1, fb1, fa2, fb2, fa3, fb3, mfaScratch] = GetProductWorkspace();

#if BIGMATH_NTT_MFA && BIGMATH_NTT_MFA_FUSE
      // Same fused MFA pipeline as the linear path (forward + fused
      // pointwise via on-the-fly packing, row-chunked inverse). The MFA
      // permutation cancels in the pointwise product and the inverse
      // returns natural order, so the Garner walk below is unchanged. This
      // is what makes the 2^26 cap usable: large cyclic transforms
      // previously ran whole-transform ParallelDo(6/3) units that idle
      // cores (the pre-#109 division profile: ~88% of compute in plain
      // butterfly layers).
      if (n >= BIGMATH_NTT_MFA_THRESHOLD && MfaFusedWindow(n))
      {
        for (int i = 0; i < 6; ++i) mfaScratch[i].resize(n);
        MfaPlanTree<F1> tree1;
        MfaPlanTree<F2> tree2;
        MfaPlanTree<F3> tree3;
        BuildMfaPlanTree<F1, G1>(n, tree1);
        BuildMfaPlanTree<F2, G2>(n, tree2);
        BuildMfaPlanTree<F3, G3>(n, tree3);
        fa1.resize(n); fa2.resize(n); fa3.resize(n);
        UInt *fap[3] = {fa1.data(), fa2.data(), fa3.data()};
        UInt *scrs[6] = {mfaScratch[0].data(), mfaScratch[1].data(), mfaScratch[2].data(),
                         mfaScratch[3].data(), mfaScratch[4].data(), mfaScratch[5].data()};
        PackedOperandView va(a, base), vb(b, base);
        MfaFusedForwardPointwise(va, vb, fap, scrs, n, tree1, tree2, tree3);
        MfaFusedInverse(fap, scrs, n, tree1, tree2, tree3);
      }
      else
#endif
      {
        const auto &plan1 = GetPlan<F1, G1>(n);
        const auto &plan2 = GetPlan<F2, G2>(n);
        const auto &plan3 = GetPlan<F3, G3>(n);
        CrtPackForwardBatch(a, &b, base, n, fa1, fa2, fa3, fb1, fb2, fb3, plan1, plan2, plan3);

        CrtPointwiseMulInto(fa1.data(), fa2.data(), fa3.data(),
                            fb1.data(), fb2.data(), fb3.data(), n);

        CrtInverseBatch(fa1, fa2, fa3, plan1, plan2, plan3);
      }

      // Garner-recombine and carry across exactly N coefficients (L limbs);
      // the carry that spills past limb L wraps back to limb 0 (B^L ≡ 1).
      const InvTable &inv = GarnerInverses();
      std::vector<DataT> r(L, 0);
      ULong128 carry = 0;
      if (base == Base2_64)
      {
        carry = GarnerCarry64(fa1.data(), fa2.data(), fa3.data(), (SizeT)(2 * (ULong)L), r.data(), L);
      }
      else
      {
        for (SizeT i = 0; i < L; ++i)
        {
          ULong128 total = Garner(fa1[i], fa2[i], fa3[i], inv) + carry;
          r[i] = (DataT)(total & 0xFFFFFFFFULL);
          carry = total >> 32;
        }
      }

      const ULong limbMask = (base == Base2_64) ? ~0ULL : 0xFFFFFFFFULL;
      const int limbBits = (base == Base2_64) ? 64 : 32;
      while (carry > 0)
      {
        ULong128 cur = carry;
        carry = 0;
        for (SizeT i = 0; cur > 0; ++i)
        {
          if (i == L)
          {
            carry += cur; // spilled past the top again — wrap once more
            break;
          }
          ULong128 s = (ULong128)(ULong)r[i] + (ULong)(cur & limbMask);
          r[i] = (DataT)((ULong)s & limbMask);
          cur = (cur >> limbBits) + (s >> limbBits);
        }
      }

      // r is now in [0, B^L − 1]; the single non-canonical value is B^L − 1
      // itself (all limbs max), which ≡ 0.
      bool allMax = true;
      for (SizeT i = 0; i < L; ++i)
        if ((ULong)r[i] != limbMask)
        {
          allMax = false;
          break;
        }
      if (allMax)
        return std::vector<DataT>{0};

      TrimZeros(r);
      if (r.empty())
        r.push_back(0);
      return r;
    }
  } // namespace NttCrt
} // namespace BigMath

#endif
