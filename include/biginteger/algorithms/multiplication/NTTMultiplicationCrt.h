/**
 * BigMath: Multi-prime CRT NTT multiplication.
 *
 * Alternate to Goldilocks NTT (NTTMultiplication.h). Three 30-bit primes
 * (998244353, 985661441, 754974721) with 32-bit coefficient splitting:
 * each 64-bit input limb → 2 coefficients (vs Goldilocks' 4). Halves the
 * NTT length at the cost of 3 parallel transforms + CRT reconstruction.
 *
 * Trade-off vs Goldilocks:
 *   - 3 transforms × NTT(N/2) work each = ~1.5× more transform layers,
 *     BUT each layer runs on 32-bit modular ops (single MUL + simple
 *     reduce) rather than ULong128 Goldilocks mul+reduce.
 *   - Coefficients fit in 32-bit lanes → better cache density.
 *   - Combined modulus = p1·p2·p3 ≈ 2^90, comfortably holds 32×32 product
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

#ifndef NTT_MULTIPLICATION_CRT
#define NTT_MULTIPLICATION_CRT

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
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

    // Three 30-bit NTT-friendly primes used for the CRT scheme.
    // Each P-1 is divisible by a large power of 2, supporting NTT lengths up
    // to that power (p1: 2^23, p2: 2^22, p3: 2^24).
    constexpr UInt P1 = 998244353;
    constexpr UInt P2 = 985661441;
    constexpr UInt P3 = 754974721;
    constexpr UInt G1 = 3;
    constexpr UInt G2 = 3;
    constexpr UInt G3 = 11;

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
      std::vector<UInt> r(n / 2);
      r[0] = 1;
      for (Int k = 1; k < n / 2; ++k) r[k] = F::Mul(r[k - 1], root);
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
    template <typename F>
    inline void Forward(std::vector<UInt> &a, const Plan<F> &plan)
    {
      Int n = (Int)a.size();
      const UInt *roots = plan.forwardRoots.data();
      UInt *aPtr = a.data();
      for (Int len = n; len > 1; len >>= 1)
      {
        Int halflen = len >> 1;
        Int stride = n / len;
        for (Int i = 0; i < n; i += len)
        {
          for (Int j = 0; j < halflen; ++j)
          {
            UInt u = aPtr[i + j];
            UInt v = aPtr[i + j + halflen];
            aPtr[i + j] = F::Add(u, v);
            aPtr[i + j + halflen] = F::Mul(F::Sub(u, v), roots[j * stride]);
          }
        }
      }
    }

    template <typename F>
    inline void Inverse(std::vector<UInt> &a, const Plan<F> &plan)
    {
      Int n = (Int)a.size();
      const UInt *roots = plan.inverseRoots.data();
      UInt *aPtr = a.data();
      for (Int len = 2; len <= n; len <<= 1)
      {
        Int halflen = len >> 1;
        Int stride = n / len;
        for (Int i = 0; i < n; i += len)
        {
          for (Int j = 0; j < halflen; ++j)
          {
            UInt u = aPtr[i + j];
            UInt v = F::Mul(aPtr[i + j + halflen], roots[j * stride]);
            aPtr[i + j] = F::Add(u, v);
            aPtr[i + j + halflen] = F::Sub(u, v);
          }
        }
      }

      UInt invSize = plan.invSize;
      for (Int i = 0; i < n; ++i) aPtr[i] = F::Mul(aPtr[i], invSize);
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
      bool b64 = (base == Base2_64);
      SizeT coeffsPerLimb = b64 ? 2u : 1u;

      ULong aCoeffSize = (ULong)a.size() * coeffsPerLimb;
      ULong bCoeffSize = (ULong)b.size() * coeffsPerLimb;
      ULong coeffCount = aCoeffSize + bCoeffSize - 1;
      Int n = (Int)std::bit_ceil(coeffCount);

      // Three parallel transforms.
      static thread_local std::vector<UInt> fa1, fb1, fa2, fb2, fa3, fb3;
      fa1.assign(n, 0); fb1.assign(n, 0);
      fa2.assign(n, 0); fb2.assign(n, 0);
      fa3.assign(n, 0); fb3.assign(n, 0);

      auto pack = [&](const std::vector<DataT> &v, std::vector<UInt> &dst1,
                      std::vector<UInt> &dst2, std::vector<UInt> &dst3) {
        if (b64)
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
      };

      pack(a, fa1, fa2, fa3);
      pack(b, fb1, fb2, fb3);

      const auto &plan1 = GetPlan<F1, G1>(n);
      const auto &plan2 = GetPlan<F2, G2>(n);
      const auto &plan3 = GetPlan<F3, G3>(n);

#if BIGMATH_USE_THREADS
      // Cross-prime parallelism: 6 forwards as one batch.
      {
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
      }
#else
      Forward<F1>(fa1, plan1); Forward<F1>(fb1, plan1);
      Forward<F2>(fa2, plan2); Forward<F2>(fb2, plan2);
      Forward<F3>(fa3, plan3); Forward<F3>(fb3, plan3);
#endif

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

#if BIGMATH_USE_THREADS
      {
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
      }
#else
      Inverse<F1>(fa1, plan1);
      Inverse<F2>(fa2, plan2);
      Inverse<F3>(fa3, plan3);
#endif

      // CRT reconstruct + propagate carries into output limbs.
      const InvTable &inv = GarnerInverses();
      std::vector<DataT> result;
      result.reserve(a.size() + b.size() + 2);

      if (b64)
      {
        // 2 coefficients per output limb, each 32 bits → 64-bit limb. Carry
        // propagation needs more than the per-coeff width; use ULong128.
        ULong128 carry = 0;
        ULong limb_acc = 0;
        int slot = 0;
        auto flush = [&result, &limb_acc, &slot](bool force) {
          if (slot == 2 || (force && slot != 0))
          {
            result.push_back((DataT)limb_acc);
            limb_acc = 0;
            slot = 0;
          }
        };
        for (SizeT i = 0; i < (SizeT)coeffCount; ++i)
        {
          ULong128 total = Garner(fa1[i], fa2[i], fa3[i], inv) + carry;
          ULong digit = (ULong)(total & 0xFFFFFFFFULL);
          carry = total >> 32;
          limb_acc |= digit << (slot * 32);
          ++slot;
          flush(false);
        }
        while (carry > 0)
        {
          ULong digit = (ULong)(carry & 0xFFFFFFFFULL);
          carry >>= 32;
          limb_acc |= digit << (slot * 32);
          ++slot;
          flush(false);
        }
        flush(true);
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
  } // namespace NttCrt
} // namespace BigMath

#endif
