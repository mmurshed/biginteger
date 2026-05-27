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
 * NTT length ceiling is 2^26 ≈ 67M coefficients — supports operand pairs
 * up to ~640M decimal digits. Earlier triple (998244353, 985661441,
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
 * Gated on BIGMATH_NTT_CRT=1. Default off; selected at compile time in
 * NTTMultiplication's dispatch when the macro is set.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

// Default: CRT NTT enabled via size-gated dispatch (threshold 5000 limbs
// sum). Wins 5-15% on large mul / skewed div / parse / ToString when the
// combined operand size crosses the gate. Below the gate, Goldilocks NTT
// runs as before — zero-cost fall-through. Opt out via -DBIGMATH_NTT_CRT=0.
#ifndef BIGMATH_NTT_CRT
#define BIGMATH_NTT_CRT 1
#endif

// Default-on Matrix Fourier Algorithm (Bailey 6-step) for very large NTT
// transforms. Recursive 2D decomposition N = N1*N2 with sub-FFTs sized to
// fit L1. Targets the bandwidth-bound regime where N exceeds L2 (~32 MB on
// M1 Max ≈ 2^23 UInt entries per prime). Threshold below = transform
// length, NOT operand limb count. Opt out via -DBIGMATH_NTT_MFA=0.
#ifndef BIGMATH_NTT_MFA
#define BIGMATH_NTT_MFA 1
#endif

// Below this NTT length, MFA leaves run the existing radix-8 chain in-place.
// Picked so sub-FFTs of 2^13 = 8192 UInt entries (32 KB) fit M1 L1 (128 KB)
// with twiddles + workspace.
#ifndef BIGMATH_NTT_MFA_LEAF
#define BIGMATH_NTT_MFA_LEAF (1 << 13)
#endif

// MFA fires when transform length n exceeds this. Default 2^21 = 2,097,152
// coeffs ≈ 8 MB per-prime buffer — already past L2 at 6-buffer working set.
// Roughly corresponds to operand sums ≥ 1M Base2_64 limbs (≈ 20 M digits).
#ifndef BIGMATH_NTT_MFA_THRESHOLD
#define BIGMATH_NTT_MFA_THRESHOLD (1 << 21)
#endif

#ifndef NTT_MULTIPLICATION_CRT
#define NTT_MULTIPLICATION_CRT

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "../../common/Parallel.h"
#include "../../common/Util.h"
#include "ClassicMultiplication.h"

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

    // ─── NTT core templated by ModField ──────────────────────────────────────

    template <typename F>
    struct Plan
    {
      std::vector<UInt> forwardRoots;
      std::vector<UInt> inverseRoots;
      UInt invSize = 1;
    };

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
      p.invSize = F::Inv((UInt)n);
      return p;
    }

    template <typename F, UInt G>
    inline const Plan<F> &GetPlan(Int n)
    {
      static thread_local std::unordered_map<Int, Plan<F>> cache;
      auto it = cache.find(n);
      if (it != cache.end()) return it->second;
      return cache.emplace(n, BuildPlan<F, G>(n)).first->second;
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
      Int len = n;
      while (len >= 8)
      {
        ForwardRadix8Layer<F>(aPtr, n, len, roots);
        len >>= 3;
      }
      if (len == 4)
        ForwardRadix4Layer<F>(aPtr, n, 4, roots);
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
      if (n >= 2)
      {
        Int logn = __builtin_ctz((unsigned)n);
        Int rem = logn % 3;
        Int len;
        if (rem == 1)
        {
          InverseRadix2Layer<F>(aPtr, n, 2, roots);
          len = 16;
        }
        else if (rem == 2)
        {
          InverseRadix4Layer<F>(aPtr, n, 4, roots);
          len = 32;
        }
        else
        {
          len = 8;
        }
        while (len <= n)
        {
          InverseRadix8Layer<F>(aPtr, n, len, roots);
          len <<= 3;
        }
      }

      if (scale)
      {
        UInt invSize = plan.invSize;
        for (Int i = 0; i < n; ++i) aPtr[i] = F::Mul(aPtr[i], invSize);
      }
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

    inline std::vector<DataT> FinalizeProduct(const std::vector<UInt> &fa1,
                                              const std::vector<UInt> &fa2,
                                              const std::vector<UInt> &fa3,
                                              ULong coeffCount,
                                              BaseT base,
                                              SizeT reserveLimbs)
    {
      const InvTable &inv = GarnerInverses();
      std::vector<DataT> result;
      result.reserve(reserveLimbs);

      if (base == Base2_64)
      {
        // 2 coefficients per output limb, each 32 bits → 64-bit limb. Carry
        // propagation needs more than the per-coeff width; use ULong128.
        ULong128 carry = 0;
        SizeT i = 0;
        for (; i + 1 < (SizeT)coeffCount; i += 2)
        {
          ULong128 total = Garner(fa1[i], fa2[i], fa3[i], inv) + carry;
          ULong lo = (ULong)(total & 0xFFFFFFFFULL);
          carry = total >> 32;

          total = Garner(fa1[i + 1], fa2[i + 1], fa3[i + 1], inv) + carry;
          ULong hi = (ULong)(total & 0xFFFFFFFFULL);
          carry = total >> 32;

          result.push_back((DataT)(lo | (hi << 32)));
        }

        ULong limb_acc = 0;
        int slot = 0;
        if (i < (SizeT)coeffCount)
        {
          ULong128 total = Garner(fa1[i], fa2[i], fa3[i], inv) + carry;
          limb_acc = (ULong)(total & 0xFFFFFFFFULL);
          carry = total >> 32;
          slot = 1;
        }
        while (carry > 0)
        {
          ULong digit = (ULong)(carry & 0xFFFFFFFFULL);
          carry >>= 32;
          limb_acc |= digit << (slot * 32);
          ++slot;
          if (slot == 2)
          {
            result.push_back((DataT)limb_acc);
            limb_acc = 0;
            slot = 0;
          }
        }
        if (slot != 0)
          result.push_back((DataT)limb_acc);
      }
      else
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
      prepared.n = (Int)std::max<ULong>(2, std::bit_ceil(maxCoeffCount));

      prepared.f1.assign(prepared.n, 0);
      prepared.f2.assign(prepared.n, 0);
      prepared.f3.assign(prepared.n, 0);
      PackOperand(operand, base, prepared.f1, prepared.f2, prepared.f3);

      const auto &plan1 = GetPlan<F1, G1>(prepared.n);
      const auto &plan2 = GetPlan<F2, G2>(prepared.n);
      const auto &plan3 = GetPlan<F3, G3>(prepared.n);

#if BIGMATH_USE_THREADS
      {
        PreparedOperand *p = &prepared;
        const Plan<F1> *pl1 = &plan1;
        const Plan<F2> *pl2 = &plan2;
        const Plan<F3> *pl3 = &plan3;
        auto body = [p, pl1, pl2, pl3](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx)
          {
            switch (idx)
            {
              case 0: Forward<F1>(p->f1, *pl1); break;
              case 1: Forward<F2>(p->f2, *pl2); break;
              case 2: Forward<F3>(p->f3, *pl3); break;
            }
          }
        };
        ParallelDo(3, body);
      }
#else
      Forward<F1>(prepared.f1, plan1);
      Forward<F2>(prepared.f2, plan2);
      Forward<F3>(prepared.f3, plan3);
#endif

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

#if BIGMATH_USE_THREADS
      {
        std::vector<UInt> *bufs[3] = {&fb1, &fb2, &fb3};
        const Plan<F1> *p1 = &plan1;
        const Plan<F2> *p2 = &plan2;
        const Plan<F3> *p3 = &plan3;
        auto body = [bufs, p1, p2, p3](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx)
          {
            switch (idx)
            {
              case 0: Forward<F1>(*bufs[0], *p1); break;
              case 1: Forward<F2>(*bufs[1], *p2); break;
              case 2: Forward<F3>(*bufs[2], *p3); break;
            }
          }
        };
        ParallelDo(3, body);
      }
#else
      Forward<F1>(fb1, plan1);
      Forward<F2>(fb2, plan2);
      Forward<F3>(fb3, plan3);
#endif

      {
        UInt *p1a = fb1.data(), *p2a = fb2.data(), *p3a = fb3.data();
        const UInt *p1b = prepared.f1.data();
        const UInt *p2b = prepared.f2.data();
        const UInt *p3b = prepared.f3.data();
        auto body = [p1a, p2a, p3a, p1b, p2b, p3b](Int s, Int e) {
          for (Int i = s; i < e; ++i)
          {
            p1a[i] = F1::Mul(p1a[i], p1b[i]);
            p2a[i] = F2::Mul(p2a[i], p2b[i]);
            p3a[i] = F3::Mul(p3a[i], p3b[i]);
          }
        };
        if ((SizeT)prepared.n >= ParallelMinSize()) ParallelFor(prepared.n, body);
        else body(0, prepared.n);
      }

#if BIGMATH_USE_THREADS
      {
        std::vector<UInt> *bufs[3] = {&fb1, &fb2, &fb3};
        const Plan<F1> *p1 = &plan1;
        const Plan<F2> *p2 = &plan2;
        const Plan<F3> *p3 = &plan3;
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
      }
#else
      Inverse<F1>(fb1, plan1);
      Inverse<F2>(fb2, plan2);
      Inverse<F3>(fb3, plan3);
#endif

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

    // Pre-built plan tree for MFA. Must be populated in the main thread before
    // dispatching ParallelDo: the per-prime thread_local Plan cache in
    // GetPlan() is invisible to worker threads, and a cold-cache GetPlan call
    // from a worker would call BuildRoots → ParallelFor and re-enter the
    // (non-reentrant) pool — deadlock.
    template <typename F>
    struct MfaPlanTree
    {
      std::unordered_map<Int, const Plan<F> *> bySize;
      const Plan<F> &Get(Int sz) const { return *bySize.at(sz); }
    };

    template <typename F, UInt G>
    inline void BuildMfaPlanTree(Int n, MfaPlanTree<F> &tree)
    {
      if (tree.bySize.count(n)) return;
      tree.bySize[n] = &GetPlan<F, G>(n);
      if (n > BIGMATH_NTT_MFA_LEAF)
      {
        Int n1, n2;
        MFAFactor(n, n1, n2);
        BuildMfaPlanTree<F, G>(n1, tree);
        BuildMfaPlanTree<F, G>(n2, tree);
      }
    }

    template <typename F>
    inline void ForwardMFA(UInt *a, Int n, UInt *scratch, bool parallel,
                           const MfaPlanTree<F> &tree);

    template <typename F>
    inline void InverseMFA(UInt *a, Int n, UInt *scratch, bool parallel,
                           const MfaPlanTree<F> &tree);

    // Backward-compat wrappers that build a fresh plan tree per call. Safe
    // only when called from a context that is NOT already inside ParallelDo.
    template <typename F, UInt G>
    inline void ForwardMFA(UInt *a, Int n, UInt *scratch, bool parallel = true)
    {
      MfaPlanTree<F> tree;
      BuildMfaPlanTree<F, G>(n, tree);
      ForwardMFA<F>(a, n, scratch, parallel, tree);
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

    // Cross-twiddle for the forward (forwardRoots) or inverse (inverseRoots)
    // pass. After a leaf DIF sub-FFT on each row, storage position k2_br
    // holds the value for *logical* k2 = br(k2_br). The correct factor at
    // storage [i, k2_br] is therefore ω_N^{±i·br(k2_br)}.
    //
    // Iterate logical k2_nat naturally and write to storage br_table[k2_nat]:
    // the exponent i·k2_nat stays monotone and bounded by i·(n2−1) < n, so no
    // wrap is needed. The leaf inverse DIT consumes the same bit-reversed
    // layout, so the inverse twiddle pass uses the identical addressing.
    template <typename F>
    inline void MfaTwiddleApply(UInt *plane, Int n1, Int n2, Int nFull, const UInt *roots, bool parallel)
    {
      Int half = nFull >> 1;
      UInt P = F::Prime;
      const std::vector<Int> &brTable = GetBitReverseTable(n2);
      const Int *br = brTable.data();
      auto body = [plane, n2, half, roots, P, br](Int iStart, Int iEnd) {
        for (Int i = iStart; i < iEnd; ++i)
        {
          UInt *row = plane + (SizeT)i * n2;
          Int idx = 0;
          for (Int k2 = 0; k2 < n2; ++k2)
          {
            UInt w = (idx < half) ? roots[idx] : (UInt)(P - roots[idx - half]);
            Int pos = br[k2];
            row[pos] = F::Mul(row[pos], w);
            idx += i;
          }
        }
      };
#if BIGMATH_USE_THREADS
      if (parallel && (SizeT)n1 >= ParallelMinSize()) ParallelFor(n1, body);
      else body(0, n1);
#else
      body(0, n1);
#endif
    }

    template <typename F>
    inline void ForwardMFA(UInt *a, Int n, UInt *scratch, bool parallel,
                           const MfaPlanTree<F> &tree)
    {
      if (n <= BIGMATH_NTT_MFA_LEAF)
      {
        ForwardPtr<F>(a, n, tree.Get(n));
        return;
      }

      Int n1, n2;
      MFAFactor(n, n1, n2);

      // Step 1: a (n2×n1) → scratch (n1×n2)
      Transpose(a, scratch, n2, n1);

      // Step 2: n1 forward sub-FFTs of length n2 on rows of scratch.
      // Each sub-call uses the matching slice of `a` as its own scratch.
      const Plan<F> &planN2 = tree.Get(n2);
      auto step2 = [scratch, a, n2, parallel, &planN2, &tree](Int rStart, Int rEnd) {
        for (Int r = rStart; r < rEnd; ++r)
        {
          UInt *row = scratch + (SizeT)r * n2;
          if (n2 <= BIGMATH_NTT_MFA_LEAF)
            ForwardPtr<F>(row, n2, planN2);
          else
            ForwardMFA<F>(row, n2, a + (SizeT)r * n2, parallel, tree);
        }
      };
#if BIGMATH_USE_THREADS
      if (parallel && (SizeT)n1 >= ParallelMinSize()) ParallelFor(n1, step2);
      else step2(0, n1);
#else
      step2(0, n1);
#endif

      // Step 3: cross-twiddle in scratch using length-n forward roots.
      const Plan<F> &planN = tree.Get(n);
      MfaTwiddleApply<F>(scratch, n1, n2, n, planN.forwardRoots.data(), parallel);

      // Step 4: scratch (n1×n2) → a (n2×n1)
      Transpose(scratch, a, n1, n2);

      // Step 5: n2 forward sub-FFTs of length n1 on rows of a.
      const Plan<F> &planN1 = tree.Get(n1);
      auto step5 = [a, scratch, n1, parallel, &planN1, &tree](Int rStart, Int rEnd) {
        for (Int r = rStart; r < rEnd; ++r)
        {
          UInt *row = a + (SizeT)r * n1;
          if (n1 <= BIGMATH_NTT_MFA_LEAF)
            ForwardPtr<F>(row, n1, planN1);
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

    template <typename F>
    inline void InverseMFA(UInt *a, Int n, UInt *scratch, bool parallel,
                           const MfaPlanTree<F> &tree)
    {
      if (n <= BIGMATH_NTT_MFA_LEAF)
      {
        InversePtr<F>(a, n, tree.Get(n), /*scale=*/ true);
        return;
      }

      Int n1, n2;
      MFAFactor(n, n1, n2);

      // Reverse step 5: n2 inverse sub-FFTs of length n1 on rows of a.
      const Plan<F> &planN1 = tree.Get(n1);
      auto invStep5 = [a, scratch, n1, parallel, &planN1, &tree](Int rStart, Int rEnd) {
        for (Int r = rStart; r < rEnd; ++r)
        {
          UInt *row = a + (SizeT)r * n1;
          if (n1 <= BIGMATH_NTT_MFA_LEAF)
            InversePtr<F>(row, n1, planN1, /*scale=*/ true);
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

      // Reverse step 3: inverse cross-twiddle in scratch.
      const Plan<F> &planN = tree.Get(n);
      MfaTwiddleApply<F>(scratch, n1, n2, n, planN.inverseRoots.data(), parallel);

      // Reverse step 2: n1 inverse sub-FFTs of length n2 on rows of scratch.
      const Plan<F> &planN2 = tree.Get(n2);
      auto invStep2 = [scratch, a, n2, parallel, &planN2, &tree](Int rStart, Int rEnd) {
        for (Int r = rStart; r < rEnd; ++r)
        {
          UInt *row = scratch + (SizeT)r * n2;
          if (n2 <= BIGMATH_NTT_MFA_LEAF)
            InversePtr<F>(row, n2, planN2, /*scale=*/ true);
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

    // ─── Public Multiply ─────────────────────────────────────────────────────

    inline std::vector<DataT> Multiply(const std::vector<DataT> &a,
                                       const std::vector<DataT> &b,
                                       BaseT base)
    {
      if (IsZero(a) || IsZero(b)) return std::vector<DataT>();
      if (a.size() == 1) return ClassicMultiplication::Multiply(b, a[0], base);
      if (b.size() == 1) return ClassicMultiplication::Multiply(a, b[0], base);

      // Split into 32-bit coefficients. Works for Base2_32 (1 coeff per limb)
      // and Base2_64 (2 coeffs per limb). Pre-CRT bounds: each coefficient is
      // < 2^32, convolution sum at length N is < N · 2^64.
      SizeT coeffsPerLimb = CoeffsPerLimb(base);

      ULong aCoeffSize = (ULong)a.size() * coeffsPerLimb;
      ULong bCoeffSize = (ULong)b.size() * coeffsPerLimb;
      ULong coeffCount = aCoeffSize + bCoeffSize - 1;
      Int n = (Int)std::bit_ceil(coeffCount);

      // Three parallel transforms.
      static thread_local std::vector<UInt> fa1, fb1, fa2, fb2, fa3, fb3;
      fa1.assign(n, 0); fb1.assign(n, 0);
      fa2.assign(n, 0); fb2.assign(n, 0);
      fa3.assign(n, 0); fb3.assign(n, 0);

      PackOperand(a, base, fa1, fa2, fa3);
      PackOperand(b, base, fb1, fb2, fb3);

      const auto &plan1 = GetPlan<F1, G1>(n);
      const auto &plan2 = GetPlan<F2, G2>(n);
      const auto &plan3 = GetPlan<F3, G3>(n);

#if BIGMATH_NTT_MFA
      const bool useMfa = (n >= BIGMATH_NTT_MFA_THRESHOLD);
      // Six per-task scratch buffers, scoped over forward + pointwise + inverse
      // so the first three can be reused for the three inverse transforms.
      // NOT thread_local: ParallelDo dispatches lambdas to worker threads, and
      // thread_local storage in workers is distinct from the main thread that
      // would have to assign() the buffer. Stack-local vectors with captured
      // pointers cross threads safely.
      std::vector<UInt> mfaScratch[6];
      MfaPlanTree<F1> tree1;
      MfaPlanTree<F2> tree2;
      MfaPlanTree<F3> tree3;
      if (useMfa)
      {
        for (int i = 0; i < 6; ++i) mfaScratch[i].assign(n, 0);
        // Pre-warm all plans in main thread: worker threads cannot call
        // GetPlan() safely from inside ParallelDo (BuildRoots reenters pool).
        BuildMfaPlanTree<F1, G1>(n, tree1);
        BuildMfaPlanTree<F2, G2>(n, tree2);
        BuildMfaPlanTree<F3, G3>(n, tree3);
        UInt *bufs[6]   = {fa1.data(), fb1.data(), fa2.data(), fb2.data(), fa3.data(), fb3.data()};
        UInt *scrs[6]   = {mfaScratch[0].data(), mfaScratch[1].data(), mfaScratch[2].data(),
                           mfaScratch[3].data(), mfaScratch[4].data(), mfaScratch[5].data()};
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
      else
#endif
      {
#if BIGMATH_USE_THREADS
        // Cross-prime parallelism: 6 forwards as one batch.
        std::vector<UInt> *bufs[6] = {&fa1, &fb1, &fa2, &fb2, &fa3, &fb3};
        const Plan<F1> *p1 = &plan1;
        const Plan<F2> *p2 = &plan2;
        const Plan<F3> *p3 = &plan3;
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
#else
        Forward<F1>(fa1, plan1); Forward<F1>(fb1, plan1);
        Forward<F2>(fa2, plan2); Forward<F2>(fb2, plan2);
        Forward<F3>(fa3, plan3); Forward<F3>(fb3, plan3);
#endif
      }

      {
        UInt *p1a = fa1.data(), *p1b = fb1.data();
        UInt *p2a = fa2.data(), *p2b = fb2.data();
        UInt *p3a = fa3.data(), *p3b = fb3.data();
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

#if BIGMATH_NTT_MFA
      if (useMfa)
      {
        // Reuse the first three forward scratches; the other three are freed
        // implicitly when the function returns. Each task gets its own.
        UInt *bufs[3] = {fa1.data(), fa2.data(), fa3.data()};
        UInt *scrs[3] = {mfaScratch[0].data(), mfaScratch[1].data(), mfaScratch[2].data()};
        auto invBody = [bufs, scrs, n, &tree1, &tree2, &tree3](Int s, Int e) {
          for (Int idx = s; idx < e; ++idx)
          {
            switch (idx)
            {
              case 0: InverseMFA<F1>(bufs[0], n, scrs[0], /*parallel=*/false, tree1); break;
              case 1: InverseMFA<F2>(bufs[1], n, scrs[1], /*parallel=*/false, tree2); break;
              case 2: InverseMFA<F3>(bufs[2], n, scrs[2], /*parallel=*/false, tree3); break;
            }
          }
        };
        ParallelDo(3, invBody);
      }
      else
#endif
      {
#if BIGMATH_USE_THREADS
        std::vector<UInt> *bufs[3] = {&fa1, &fa2, &fa3};
        const Plan<F1> *p1 = &plan1;
        const Plan<F2> *p2 = &plan2;
        const Plan<F3> *p3 = &plan3;
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
        Inverse<F1>(fa1, plan1);
        Inverse<F2>(fa2, plan2);
        Inverse<F3>(fa3, plan3);
#endif
      }

      return FinalizeProduct(fa1, fa2, fa3, coeffCount, base, a.size() + b.size() + 2);
    }
  } // namespace NttCrt
} // namespace BigMath

#endif
