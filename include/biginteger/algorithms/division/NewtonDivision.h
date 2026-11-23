#ifndef NEWTON_DIVISION
#define NEWTON_DIVISION

#include <algorithm>
#include <bit>
#include <cstring>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../../common/BitShifts.h"
#include "../../common/Comparator.h"
#include "../../common/Util.h"
#include "../Addition.h"
#include "../Multiplication.h"
#include "../Subtraction.h"
#include "ClassicDivision.h"
#include "FastDivision.h"

// Gate for the cyclic (wrap-around) NTT products inside Newton division.
// Separate from NTT_MULTIPLICATION_THRESHOLD because a cyclic product runs at
// HALF the transform length of the full product the threshold was tuned for —
// its crossover vs Karatsuba sits proportionally lower. Swept 2026-06-11 on
// cold ToString (whose divider-chain reciprocal towers run entirely in this
// band at <= 500k digits): 1280 beats 5120 by ~30% at 100k-200k digits and
// is neutral at 500k+; 640 regresses. Tunable for re-sweeps.
#ifndef BIGMATH_CYCLIC_NTT_THRESHOLD
#define BIGMATH_CYCLIC_NTT_THRESHOLD 1280
#endif

// Gate for the top-limbs quotient estimate (the tryApprox path below). It
// used to read NTT_MULTIPLICATION_THRESHOLD directly; that constant moved
// 1280 -> 640 in the 2026-11-20 multiplication retune, and the division
// frontier has not been re-swept against the faster multiply, so the gate
// keeps the value it was tuned at.
#ifndef BIGMATH_NEWTON_APPROX_THRESHOLD
#define BIGMATH_NEWTON_APPROX_THRESHOLD 1280
#endif

// Blockwise division: a leading block with at most this many quotient limbs
// goes to schoolbook instead of a Newton chunk. A chunk costs a reciprocal
// multiply and a back-multiply at the divisor's size whatever the block's
// quotient length — 60-75 ns per divisor limb from the Newton floor up to
// 65k limbs — while schoolbook costs ~1.3 ns per divisor limb per quotient
// limb. The case that matters is one limb: the normalization shift pushes a
// limb out of a dividend of exactly k*n limbs for one random input in three.
#ifndef BIGMATH_NEWTON_SCHOOLBOOK_HEAD
#define BIGMATH_NEWTON_SCHOOLBOOK_HEAD 32
#endif

namespace BigMath
{
  // Newton-Raphson division.
  // Precomputes an n-limb approximate reciprocal R of the normalized divisor D
  // such that R*D ≈ B^(2n). For na ≤ 2n: one Q = (a*R) >> 2n + small fixup.
  // For na > 2n: blockwise — process top first_chunk ∈ [n+1, 2n] limbs, then slide down
  // by n, threading the remainder as the high part of each next chunk. Cost stays O(M(n))
  // per block, so blockwise total = (na/n) · M(n) ≈ O(M(na)) — beats Knuth's O(na·n) and
  // BZ's O(M(na)·log) at sizes where M(n) is NTT-dominated.
  class NewtonDivision
  {
  public:
    // Number of times this thread's Newton divisions gave up on the
    // reciprocal path (a quotient estimate outside its fixup budget, or a
    // block-size invariant violated) and fell back to schoolbook division.
    // The result is still correct when that happens, only quadratic — which
    // is exactly why it needs counting: a silent fallback once cost an 8x
    // ToString regression. Tests assert it stays at zero.
    static long &FallbackCount()
    {
      static thread_local long count = 0;
      return count;
    }

  private:
    static std::pair<std::vector<DataT>, std::vector<DataT>> Fallback(
        std::span<const DataT> a, std::span<const DataT> b, BaseT base, bool computeRemainder)
    {
      ++FallbackCount();
      return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);
    }

    struct ScratchBuffers
    {
      std::vector<DataT> v0;
      std::vector<DataT> v1;
      std::vector<DataT> v2;
      std::vector<DataT> v3;
      std::vector<DataT> v4;
      std::vector<DataT> v5;
    };

    static ScratchBuffers &Scratch()
    {
      static thread_local ScratchBuffers scratch;
      return scratch;
    }

    // Keep already-normalized inputs as views. Shifted inputs own their
    // storage for the duration of the divide or reciprocal setup.
    static std::span<const DataT> Normalize(std::span<const DataT> input, Int shift,
                                            Int limbBits, std::vector<DataT> &storage)
    {
      if (shift == 0)
        return input;
      storage = ShiftLeftBits(input, shift, limbBits);
      return storage;
    }

  public:
    // ApproxReciprocal: given n-limb normalized D (top bit of D[n-1] set),
    // returns R ≈ floor(B^(2n) / D), n+1 limbs.
    //
    // Contract on 64-bit limbs (pinned by tests/unit/test_newton_division.cpp):
    //   0 ≤ floor(B^(2n) / D) − R ≤ 3
    // i.e. R never overestimates and is within 3 ulps — at every n, without
    // the extra iteration. Implementation: exact reciprocal of the top ≤ 32
    // limbs by schoolbook division, then Newton steps along a guarded
    // precision schedule (see the body).
    //
    // `high_precision` runs one more Newton step at full precision. The
    // legacy 32-in-64 layout still needs it for chunks of 2n limbs (its tower
    // doubles exactly and drifts); 64-bit callers no longer pass it.
    static std::vector<DataT> ApproxReciprocal(std::vector<DataT> const &D, bool high_precision = false)
    {
      return ApproxReciprocalLimbs(D, high_precision);
    }

  private:
    static std::vector<DataT> ApproxReciprocalLimbs(std::span<const DataT> D, bool high_precision)
    {
      SizeT n = (SizeT)D.size();
      auto &scratch = Scratch();
      std::vector<DataT> &R_pad = scratch.v0;
      std::vector<DataT> &T = scratch.v2;
      std::vector<DataT> &two_S = scratch.v3;
      std::vector<DataT> &diff = scratch.v4;
      std::vector<DataT> &RD = scratch.v5;

      // n == 1 base case: R ≈ B^2 / D[0].
      if (n == 1)
      {
#if BIGMATH_LIMB_64
        // B^2 = 2^128 doesn't fit in ULong128; use (2^128 - 1) / D[0]. Result
        // is in [2^64, 2^65 - 1] for normalized D[0] ∈ [2^63, 2^64 - 1], so up
        // to 2 limbs. Newton fixup loop corrects the off-by-one when relevant.
        ULong128 num = ~(ULong128)0;
        ULong128 R = num / D[0];
        std::vector<DataT> result;
        result.push_back((DataT)R);
        if ((ULong)(R >> 64))
          result.push_back((DataT)(R >> 64));
        TrimZerosToOne(result);
        return result;
#else
        ULong128 num = (ULong128)1 << 64; // B^2
        ULong128 R = num / D[0];
        std::vector<DataT> result;
        result.push_back((DataT)(R & 0xFFFFFFFFULL));
        result.push_back((DataT)((R >> 32) & 0xFFFFFFFFULL));
        if (R >> 64)
          result.push_back((DataT)(R >> 64));
        TrimZerosToOne(result);
        return result;
#endif
      }

#if BIGMATH_LIMB_64
      // 64-bit limbs: precision schedule with guard limbs.
      //
      // A Newton step squares the relative error, so a seed good to p limbs
      // yields 2p limbs minus a few bits. The old schedule asked for exactly
      // 2p each time, so that shortfall doubled at every step — the tower
      // drifted by ~B^6 at n = 1024 (measured against floor(B^2n / D)) and
      // callers needing a tight reciprocal paid for one more full-size
      // iteration to repair it. Instead: work back from n, each earlier
      // precision being ceil((p + G) / 2), so every step's doubled accuracy
      // covers its target with G limbs to spare, and start from an EXACT
      // reciprocal of the top limbs (schoolbook division) rather than a
      // one-limb estimate. The result is good to a few ulps at every size
      // without the extra iteration.
      constexpr SizeT RecipGuard = 2;
      constexpr SizeT RecipBase = 32;
      SizeT schedule[72];
      int steps = 0;
      for (SizeT pcur = n;; pcur = (pcur + RecipGuard + 1) / 2)
      {
        schedule[steps++] = pcur;
        if (pcur <= RecipBase)
          break;
      }

      // Base: R = floor(B^(2p) / D_top) for the top p limbs of D.
      SizeT cur_n = schedule[--steps];
      std::vector<DataT> R;
      {
        std::vector<DataT> num(2 * cur_n + 1, 0);
        num[2 * cur_n] = 1;
        R = FastDivision::DivideAndRemainder(std::span<const DataT>(num), D.last(cur_n),
                                             CurrentBase, false).first;
      }
#else
      // Bootstrap from top 2 limbs: R_seed ≈ B^4 / D_top.
      ULong128 D_top = ((ULong128)D[n - 1] << 32) | D[n - 2];
      ULong128 R_seed = ((ULong128)-1) / D_top; // floor((2^128 - 1) / D_top)

      std::vector<DataT> R;
      R.push_back((DataT)(R_seed & 0xFFFFFFFFULL));
      R.push_back((DataT)((R_seed >> 32) & 0xFFFFFFFFULL));
      if (R_seed >> 64)
        R.push_back((DataT)(R_seed >> 64));

      SizeT cur_n = 2;
#endif
      SizeT extra_iters_done = 0;
      const SizeT EXTRA_REFINE_ITERS = high_precision ? 1 : 0;
      while (cur_n < n || extra_iters_done < EXTRA_REFINE_ITERS)
      {
        SizeT new_n;
        SizeT extend;
        if (cur_n < n)
        {
#if BIGMATH_LIMB_64
          new_n = schedule[--steps];
#else
          new_n = std::min((SizeT)(2 * cur_n), n);
#endif
          extend = new_n - cur_n;
        }
        else
        {
          new_n = n;
          extend = 0;
          ++extra_iters_done;
        }

#if BIGMATH_LIMB_64
        if (CurrentBase == Base2_64)
        {
          // One Newton step from precision cur_n to new_n = m, on the
          // UNPADDED seed. With e = m − cur and R_pad = R·B^e the exact step
          //   R_new = R_pad + floor(R_pad · (B^(2m) − D_m·R_pad) / B^(2m))
          // is, since B^(2m) − D_m·R_pad = B^e · E' with
          //   E' = B^(m+cur) − D_m · R,
          // exactly
          //   R_new = R·B^e + floor(R · E' / B^(2cur)).
          // So the step needs (1) the m × (cur+1) product D_m·R — half the
          // size of the padded m × (m+1) product — and (2) the top e+2 limbs
          // of R·E', which only the top e+4 limbs of each factor can reach.
          // The padded form multiplied (m+1) × (2m+1) for (2).
          //
          // |E'| < 4·B^m (the seed is good to a few ulps at cur limbs), so
          // E' fits m+1 limbs; that is what the cyclic window below relies on.
          const SizeT m = new_n;
          const SizeT cur = cur_n;
          const SizeT e = m - cur;
          const auto dTop = D.last(m);
          TrimZerosToOne(R);

          // (1) E' with its sign. Either the plain product, or — when a
          // cyclic transform is strictly shorter — its residue mod B^L − 1,
          // which determines E' because |E'| < B^(m+2) ≤ B^(L−4).
          std::vector<DataT> E;
          bool eNeg = false;
          bool haveE = false;
          {
            const ULong linCoeffs = (ULong)2 * ((ULong)m + R.size());
            const ULong linFull = std::bit_ceil(linCoeffs);
            ULong linCost = linFull;
            {
              const ULong H = linFull / 2;
              const Int Lt = NttCrt::TruncCyclicLength(linCoeffs, (Int)H);
              if (Lt != 0)
                linCost = H + (ULong)Lt;
            }
            const ULong cycCoeffs = std::bit_ceil((ULong)2 * ((ULong)m + 6));
            const SizeT L = (SizeT)(cycCoeffs / 2);
            if ((SizeT)m + (SizeT)R.size() >= BIGMATH_CYCLIC_NTT_THRESHOLD &&
                cycCoeffs < linCost && cycCoeffs <= (1u << 26) && R.size() <= L)
            {
              std::vector<DataT> W = NttCrt::MultiplyMod2km1(dTop, R, L, CurrentBase);
              const std::vector<DataT> M(L, ~(DataT)0); // B^L − 1
              const SizeT rexp = (SizeT)(((ULong)m + cur) % L);
              std::vector<DataT> Br(rexp + 1, 0);
              Br[rexp] = 1;
              E = (Compare(Br, W) >= 0) ? Subtract(Br, W, CurrentBase)
                                        : Subtract(Add(Br, M, CurrentBase), W, CurrentBase);
              // E' ≥ 0 lands below B^(m+2); E' < 0 lands within B^(m+2) of M.
              eNeg = (TrimmedSize(E) > m + 2);
              if (eNeg)
                E = Subtract(M, E, CurrentBase);
              haveE = (TrimmedSize(E) <= m + 2);
            }
          }
          if (!haveE)
          {
            T = detail::MultiplyLimbs(dTop, R, CurrentBase);
            two_S.assign((size_t)m + cur + 1, 0);
            two_S[(size_t)m + cur] = 1; // B^(m+cur)
            eNeg = Compare(two_S, T) < 0;
            E = eNeg ? Subtract(T, two_S, CurrentBase) : Subtract(two_S, T, CurrentBase);
          }
          TrimZerosToOne(E);

          // (2) C = floor(R · |E'| / B^(2cur)) from top slices. Dropping E'
          // below limb cur−2 and R below limb cur−e−3 perturbs C by less
          // than B^−2 each, so C is the exact value or one less.
          std::vector<DataT> C{0};
          if (!IsZero(E))
          {
            const SizeT sE = cur >= 2 ? cur - 2 : 0;
            const SizeT sR = cur >= e + 3 ? cur - e - 3 : 0;
            if (E.size() > sE && R.size() > sR)
            {
              const auto et = std::span<const DataT>(E).subspan(sE);
              const auto rt = std::span<const DataT>(R).subspan(sR);
              std::vector<DataT> P = detail::MultiplyLimbs(rt, et, CurrentBase);
              const SizeT back = 2 * cur - sE - sR;
              if (P.size() > back)
                C.assign(P.begin() + back, P.end());
            }
          }

          // R_new = R·B^e ± C, keeping R an underestimate of the reciprocal.
          R_pad.assign(R.size() + e, 0);
          std::memcpy(R_pad.data() + e, R.data(), R.size() * sizeof(DataT));
          if (!eNeg)
          {
            R = Add(R_pad, C, CurrentBase);
          }
          else
          {
            static const std::vector<DataT> two{2};
            std::vector<DataT> adj = Add(C, two, CurrentBase);
            R = (Compare(R_pad, adj) >= 0) ? Subtract(R_pad, adj, CurrentBase)
                                           : std::vector<DataT>{0};
          }
          TrimZerosToOne(R);
          cur_n = new_n;
          continue;
        }
#endif

        // Pad R as seed at new_n precision: shift up by (new_n - cur_n) limbs.
        R_pad.assign(R.size() + extend, 0);
        std::memcpy(R_pad.data() + extend, R.data(), R.size() * sizeof(DataT));

        // D_new = top new_n limbs of D.
        const auto dTop = D.last(new_n);

        // Wrapped iteration (invertappr style). Algebra: with
        // E = B^(2m) − D_new·R_pad, the exact step
        //   R_new = (R_pad · (2·B^(2m) − D_new·R_pad)) >> 2m
        // decomposes EXACTLY into R_new = R_pad + floor(R_pad·E / B^(2m)).
        // The seed is accurate to cur_n limbs, so |E| < ~B^(2m − cur_n + 1):
        //   (1) E is recovered exactly from the residue (D_new·R_pad) mod
        //       (B^L − 1) at L ≈ 2m − cur_n — a cyclic transform SHORTER than
        //       the full 2m+1-limb product (half at the extra-refine iter,
        //       where cur_n = m). Sign read from residue magnitude, same
        //       window argument as WrappedRemainder.
        //   (2) the correction floor(R_pad·E / B^(2m)) has only m − cur_n + 1
        //       significant limbs, so only the top limbs of R_pad and E
        //       contribute — a (m−cur)×(m−cur)-ish product instead of the
        //       (m+1)×(2m+1) full RD. Truncation costs ≤ ~3 ulp, absorbed by
        //       Newton's quadratic self-correction (non-final iters) or the
        //       DivideChunk fixup budget (final iter).
        bool wrappedIter = false;
#if BIGMATH_NTT_CRT
        if (CurrentBase == Base2_32 || CurrentBase == Base2_64)
        {
          SizeT m = new_n;
          SizeT c = (CurrentBase == Base2_64) ? 2 : 1;
          SizeT Lmin = 2 * m - cur_n + 16;
          ULong nCyc = std::bit_ceil((ULong)(Lmin + 1) * c);
          ULong nLin = std::bit_ceil(((ULong)dTop.size() + R_pad.size()) * c);
          SizeT L = (SizeT)(nCyc / c);
          if (dTop.size() + R_pad.size() >= BIGMATH_CYCLIC_NTT_THRESHOLD &&
              nCyc < nLin && nCyc <= (1u << 26) &&
              dTop.size() <= L && R_pad.size() <= L)
          {
            wrappedIter = true;

            std::vector<DataT> W = NttCrt::MultiplyMod2km1(dTop, R_pad, L, CurrentBase);

            const DataT maxLimb = (CurrentBase == Base2_64) ? (DataT)~0ULL : (DataT)0xFFFFFFFFULL;
            const std::vector<DataT> M(L, maxLimb); // B^L − 1

            // E = (B^(2m) − W) mod M; B^(2m) ≡ B^(2m mod L).
            SizeT r_exp = (SizeT)((2 * (ULong)m) % L);
            std::vector<DataT> Br(r_exp + 1, 0);
            Br[r_exp] = 1;
            std::vector<DataT> E = (Compare(Br, W) >= 0)
                                  ? Subtract(Br, W, CurrentBase)
                                  : Subtract(Add(Br, M, CurrentBase), W, CurrentBase);

            // Sign window: E ≥ 0 lands below B^Lmin, E < 0 lands within
            // B^Lmin of M. (L ≥ Lmin + 1 by construction of nCyc.)
            bool eNeg = (TrimmedSize(E) > Lmin);
            if (eNeg)
              E = Subtract(M, E, CurrentBase); // |E|
            if (TrimmedSize(E) > Lmin)
            {
              // Seed drift exceeds the window — can't trust the sign read.
              // Run this step through the exact path instead.
              wrappedIter = false;
            }

            if (wrappedIter)
            {
              // C = floor(R_pad·E / B^(2m)) from top slices only. Slice
              // guards (20 limbs) sized to the B^16 drift window so the
              // dropped tails stay sub-ulp.
              SizeT xl = (cur_n > 20 && cur_n - 20 < (SizeT)R_pad.size()) ? cur_n - 20 : 0;
              SizeT el = (m > 20) ? m - 20 : 0;
              std::vector<DataT> C{0};
              if (!IsZero(E) && TrimmedSize(E) > el)
              {
                const auto rt = std::span<const DataT>(R_pad).subspan(xl);
                const auto et = std::span<const DataT>(E).subspan(el);
                std::vector<DataT> P = detail::MultiplyLimbs(rt, et, CurrentBase);
                SizeT back = 2 * m - xl - el;
                if (P.size() > back)
                  C.assign(P.begin() + back, P.end());
              }

              if (!eNeg)
              {
                R = Add(R_pad, C, CurrentBase);
              }
              else
              {
                // Subtract with a +2 guard so R keeps the underestimate
                // invariant despite the truncated (under-counted) correction.
                static const std::vector<DataT> two{2};
                std::vector<DataT> adj = Add(C, two, CurrentBase);
                R = (Compare(R_pad, adj) >= 0) ? Subtract(R_pad, adj, CurrentBase)
                                               : std::vector<DataT>{0};
              }
              TrimZerosToOne(R);
            }
          }
        }
#endif

        if (!wrappedIter)
        {
          // T = D_new * R_pad
          T = detail::MultiplyLimbs(dTop, R_pad, CurrentBase);

          // diff = 2*B^(2*new_n) - T  (should be ≥ 0 for valid seeds).
          two_S.assign(2 * new_n + 1, 0);
          two_S[2 * new_n] = 2;

          if (Compare(two_S, T) >= 0)
          {
            diff = Subtract(two_S, T, CurrentBase);
          }
          else
          {
            // R was over-estimate; shouldn't normally trigger. Clamp to 0.
            diff.assign(1, 0);
          }

          // R_new = (R_pad * diff) >> (32 * 2 * new_n)
          RD = Multiply(R_pad, diff, CurrentBase);
          if (RD.size() > 2 * new_n)
          {
            R.assign(RD.begin() + 2 * new_n, RD.end());
          }
          else
          {
            R = std::vector<DataT>{0};
          }
          TrimZerosToOne(R);
        }

        cur_n = new_n;
      }

      return R;
    }

  private:
    // Sentinel for "fixup loop blew the cap" — caller should fall back to FastDivision.
    static SizeT TrimmedSize(std::vector<DataT> const &v)
    {
      SizeT s = (SizeT)v.size();
      while (s > 1 && v[s - 1] == 0)
        --s;
      return s;
    }

    // Wrap-around remainder (GMP mu_div-style): instead of the full product
    // QB = Q·b_norm, compute only its residue W = (Q·b_norm) mod (B^L − 1)
    // with a cyclic NTT of HALF the transform length, and reconstruct
    // rem = chunk − Q·b_norm from it.
    //
    // Exactness: Q is the truncated reciprocal estimate, so the true
    // t = chunk − Q·b_norm satisfies |t| < 9·b_norm < B^(n+1) ≤ B^(L−1) — the
    // residue (chunk − W) mod (B^L − 1) identifies t uniquely, and its sign is
    // readable from the magnitude: t ≥ 0 lands in [0, B^(n+1)) (≤ n+1 limbs),
    // t < 0 lands in (M − 9·B^n, M) (≥ n+2 limbs). L ≥ n+2 by construction.
    //
    // Returns false if a fixup cap blows; caller falls back to FastDivision
    // exactly like the plain path. `fixupLimit` must stay ≪ B so the
    // magnitude window (|t| < (fixupLimit+1)·b_norm < B^(n+1)) holds.
    static bool WrappedRemainder(
        std::span<const DataT> chunk,
        std::span<const DataT> b_norm,
        SizeT L,
        std::vector<DataT> &Q,
        std::vector<DataT> &rem,
        int fixupLimit)
    {
      SizeT n = (SizeT)b_norm.size();
      const DataT maxLimb = (CurrentBase == Base2_64) ? (DataT)~0ULL : (DataT)0xFFFFFFFFULL;
      const std::vector<DataT> M(L, maxLimb); // B^L − 1

      std::vector<DataT> W = NttCrt::MultiplyMod2km1(Q, b_norm, L, CurrentBase);

      // chunk mod M: fold the limbs above L back onto the low L (B^L ≡ 1).
      std::vector<DataT> cm(chunk.begin(), chunk.begin() + std::min((SizeT)chunk.size(), L));
      if (chunk.size() > L)
      {
        cm = detail::AddLimbs(cm, chunk.subspan(L), CurrentBase);
        while (Compare(cm, M) >= 0)
          cm = Subtract(cm, M, CurrentBase);
      }

      // rem_m = (cm − W) mod M
      std::vector<DataT> rem_m = (Compare(cm, W) >= 0)
                                ? Subtract(cm, W, CurrentBase)
                                : Subtract(Add(cm, M, CurrentBase), W, CurrentBase);

      static const std::vector<DataT> one{1};

      // t < 0 (Q overestimated): rem_m ≈ M − |t|, which needs > n+1 limbs.
      int iters = 0;
      while (TrimmedSize(rem_m) > n + 1)
      {
        if (++iters > fixupLimit)
          return false;
        Q = Subtract(Q, one, CurrentBase);
        rem_m = detail::AddLimbs(rem_m, b_norm, CurrentBase);
        while (Compare(rem_m, M) >= 0)
          rem_m = Subtract(rem_m, M, CurrentBase);
      }

      // rem_m now equals chunk − Q·b_norm exactly.
      iters = 0;
      while (Compare(rem_m, b_norm) >= 0)
      {
        if (++iters > fixupLimit)
          return false;
        rem_m = detail::SubtractLimbs(rem_m, b_norm, CurrentBase);
        Q = Add(Q, one, CurrentBase);
      }

      rem = std::move(rem_m);
      TrimZerosToOne(rem);
      return true;
    }

    // Remainder + quotient fixups for one chunk: computes rem = chunk − Q·b_norm,
    // adjusting Q until 0 ≤ rem < b_norm. Routes through the half-length cyclic
    // product when the full Q·b_norm would be CRT-NTT-routed and the cyclic
    // transform is strictly shorter; plain full product otherwise. Returns
    // false when the fixup budget blows.
    static bool RemainderAndFixups(
        std::span<const DataT> chunk,
        std::span<const DataT> b_norm,
        std::vector<DataT> &Q,
        std::vector<DataT> &rem,
        int fixupLimit)
    {
      SizeT n = (SizeT)b_norm.size();

#if BIGMATH_NTT_CRT
      if ((CurrentBase == Base2_32 || CurrentBase == Base2_64) && !IsZero(Q))
      {
        SizeT c = (CurrentBase == Base2_64) ? 2 : 1;
        ULong nCyc = std::bit_ceil((ULong)(n + 2) * c);
        ULong nLinear = std::bit_ceil(((ULong)Q.size() + n) * c);
        SizeT L = (SizeT)(nCyc / c);
        if (Q.size() + n >= BIGMATH_CYCLIC_NTT_THRESHOLD &&
            nCyc < nLinear && nCyc <= (1u << 26) &&
            Q.size() <= L)
          return WrappedRemainder(chunk, b_norm, L, Q, rem, fixupLimit);
      }
#endif

      std::vector<DataT> &QB = Scratch().v2;
      QB = detail::MultiplyLimbs(Q, b_norm, CurrentBase);

      static const std::vector<DataT> one{1};

      int iters = 0;
      while (Compare(QB, chunk) > 0)
      {
        if (++iters > fixupLimit)
          return false;
        Q = Subtract(Q, one, CurrentBase);
        QB = detail::SubtractLimbs(QB, b_norm, CurrentBase);
      }
      rem = detail::SubtractLimbs(chunk, QB, CurrentBase);

      iters = 0;
      while (Compare(rem, b_norm) >= 0)
      {
        if (++iters > fixupLimit)
          return false;
        rem = detail::SubtractLimbs(rem, b_norm, CurrentBase);
        Q = Add(Q, one, CurrentBase);
      }

      return true;
    }

    // The quotient and normalized remainder are valid on success.
    static bool DivideChunk(
        std::span<const DataT> chunk, std::span<const DataT> b_norm,
        std::span<const DataT> R, std::vector<DataT> &Q, std::vector<DataT> &rem)
    {
      SizeT n = (SizeT)b_norm.size();
      std::vector<DataT> &CR = Scratch().v0;

      // Approximate quotient estimate from the TOP n+1 limbs of the chunk
      // (GMP mu_divappr style). b_norm is normalized (top bit set, so
      // b ≥ B^n/2), which bounds the dropped low limbs' contribution to Q by
      // c_lo / b < 2·B^(sh−n) ≤ 2/B — under 1 ulp; with the floor truncations
      // the estimate underestimates Q by ≤ ~3 extra steps, absorbed by the
      // fixup loop. The product shrinks from (chunk + R) to (n+1 + R) limbs;
      // gated on the smaller transform actually being shorter. If the relaxed
      // fixup budget ever blows, retry once with the exact full product before
      // falling back to FastDivision.
      bool tryApprox = false;
#if BIGMATH_NTT_CRT
      if ((CurrentBase == Base2_32 || CurrentBase == Base2_64) &&
          chunk.size() > n + 1)
      {
        SizeT c = (CurrentBase == Base2_64) ? 2 : 1;
        ULong nFull = std::bit_ceil(((ULong)chunk.size() + R.size()) * c);
        ULong nTop = std::bit_ceil(((ULong)n + 1 + R.size()) * c);
        tryApprox = (chunk.size() + R.size() >= (SizeT)BIGMATH_NEWTON_APPROX_THRESHOLD) &&
                    nTop < nFull;
      }
#endif

      for (int attempt = tryApprox ? 0 : 1; attempt < 2; ++attempt)
      {
        SizeT drop = 2 * n;
        if (attempt == 0)
        {
          SizeT sh = (SizeT)chunk.size() - (n + 1);
          CR = detail::MultiplyLimbs(chunk.subspan(sh), R, CurrentBase);
          drop = 2 * n - sh;
        }
        else
        {
          // Q ≈ (chunk * R) >> (2n limbs)
          CR = detail::MultiplyLimbs(chunk, R, CurrentBase);
        }

        if (CR.size() > drop)
          Q.assign(CR.begin() + drop, CR.end());
        else
          Q.assign(1, 0);
        TrimZerosToOne(Q);

        if (RemainderAndFixups(chunk, b_norm, Q, rem, attempt == 0 ? 12 : 8))
        {
          TrimZerosToOne(Q);
          TrimZerosToOne(rem);
          return true;
        }
        // Approx attempt blew its budget — retry once with the exact product.
      }

      return false;
    }

    static std::pair<std::vector<DataT>, std::vector<DataT>> DivideNormalizedWithReciprocal(
        std::span<const DataT> a,
        std::span<const DataT> b,
        BaseT base,
        Int shift,
        std::span<const DataT> a_norm,
        std::span<const DataT> b_norm,
        std::span<const DataT> R,
        bool computeRemainder = true)
    {
      SizeT n = (SizeT)b_norm.size();
      SizeT na = (SizeT)a_norm.size();

      // Two paths:
      //   Single-block (na ≤ 2n): one reciprocal-divide on the whole a.
      //   Blockwise (na > 2n): top chunk ∈ [n+1, 2n] limbs, then slide.
      // The boundary is 2n, not 2n+1: with a 2n+1-limb chunk the truncated
      // (chunk·R)>>2n quotient estimate underestimates Q by up to ~B steps
      // (error ∝ chunk/B^2n, which reaches B once chunk exceeds 2n limbs),
      // overflowing the fixup cap and falling back to quadratic FastDivision.
      // The +1 limb routinely appears from the Knuth normalize shift on a≈2n.
      const bool blockwise = na > 2 * n;
      std::vector<DataT> rem;

      if (!blockwise)
      {
        std::vector<DataT> quotient;
        if (!DivideChunk(a_norm, b_norm, R, quotient, rem))
          return Fallback(a, b, base, computeRemainder);

        if (!computeRemainder)
          return {std::move(quotient), std::vector<DataT>{}};
        ShiftRightBitsInPlace(rem, shift, LimbBitsFor(base));
        TrimZerosToOne(rem);
        return {std::move(quotient), std::move(rem)};
      }

      // ----- Blockwise path -----
      // first_chunk_size in [n+1, 2n], chosen so remaining (na - first_chunk_size) is divisible by n.
      // Formula: ((na - 1) mod n) + 1 + n.  Verified for na in {2n+2, 3n, 3n+1, 4n, …}.
      SizeT first_chunk_size = (SizeT)(((na - 1) % n) + 1 + n);
      SizeT pos_low = na - first_chunk_size;
      auto chunk = a_norm.subspan(pos_low);
      std::vector<DataT> chunkStorage;
      std::vector<DataT> quotient(na - n + 1);
      auto &qBlock = Scratch().v1;
      bool is_first = true;

      // Short leading block: schoolbook (see BIGMATH_NEWTON_SCHOOLBOOK_HEAD).
      // b_norm is normalized, so FastDivision divides it as it stands.
      const bool schoolbookHead =
          base == Base2_64 && first_chunk_size - n <= (SizeT)BIGMATH_NEWTON_SCHOOLBOOK_HEAD;

      while (true)
      {
        if (is_first && schoolbookHead)
        {
          auto head = FastDivision::DivideAndRemainder(chunk, b_norm, base, true);
          TrimZerosToOne(head.first);
          TrimZerosToOne(head.second);
          qBlock = std::move(head.first);
          rem = std::move(head.second);
        }
        else if (!DivideChunk(chunk, b_norm, R, qBlock, rem))
          return Fallback(a, b, base, computeRemainder);

        // Every later block consumes exactly n dividend limbs. The
        // zero-filled output already supplies padding for short blocks.
        const SizeT width = is_first ? first_chunk_size - n + 1 : n;
        if (qBlock.size() > width)
          return Fallback(a, b, base, computeRemainder);
        std::copy(qBlock.begin(), qBlock.end(), quotient.begin() + pos_low);
        is_first = false;

        if (pos_low == 0)
          break;

        // Join the next n dividend limbs and the current remainder in the
        // same working buffer on every iteration.
        pos_low -= n;
        chunkStorage.resize(n + rem.size());
        std::copy_n(a_norm.data() + pos_low, n, chunkStorage.data());
        std::copy(rem.begin(), rem.end(), chunkStorage.begin() + n);
        chunk = chunkStorage;
      }

      TrimZerosToOne(quotient);
      if (!computeRemainder)
        return {std::move(quotient), std::vector<DataT>{}};
      ShiftRightBitsInPlace(rem, shift, LimbBitsFor(base));
      TrimZerosToOne(rem);
      return {std::move(quotient), std::move(rem)};
    }

  public:
    class Divider
    {
    private:
      std::vector<DataT> divisor;
      std::vector<DataT> b_norm;
      std::vector<DataT> reciprocal;
      BaseT base;
      Int shift;
      bool can_use_newton;

    public:
      Divider(std::vector<DataT> const &b, BaseT radix) : divisor(b), base(radix), shift(0), can_use_newton(false)
      {
        TrimZeros(divisor);
        if (IsZero(divisor))
          throw std::invalid_argument("Division by zero");

        if ((base != Base2_32 && base != Base2_64) || divisor.size() <= 1)
          return;

        shift = NormalizationShiftBits(divisor.back(), LimbBitsFor(base));

        const auto normalized = Normalize(divisor, shift, LimbBitsFor(base), b_norm);

        // Precompute the reciprocal once. It is valid for both single-block and blockwise
        // division and avoids recomputing the Newton setup for repeated divisions by the
        // same divisor. (No extra refinement iteration: see ApproxReciprocal's schedule.)
        reciprocal = ApproxReciprocalLimbs(normalized, /*high_precision=*/ !BIGMATH_LIMB_64);
        can_use_newton = true;
      }

      std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainder(
          std::vector<DataT> const &a,
          bool computeRemainder = true) const
      {
        if (IsZero(a))
          return {std::vector<DataT>{0}, computeRemainder ? std::vector<DataT>{0} : std::vector<DataT>()};

        Int cmp = Compare(a, divisor);
        if (cmp < 0)
          return {std::vector<DataT>{0}, computeRemainder ? a : std::vector<DataT>()};
        if (cmp == 0)
          return {std::vector<DataT>{1}, computeRemainder ? std::vector<DataT>{0} : std::vector<DataT>()};

        if (!can_use_newton)
          return FastDivision::DivideAndRemainder(a, divisor, base, computeRemainder);

        std::vector<DataT> shifted;
        const auto a_norm = Normalize(std::span<const DataT>(a).first(LimbDetail::TrimmedLimbCount(a)),
                                      shift, LimbBitsFor(base), shifted);
        const std::span<const DataT> normalizedDivisor = shift == 0 ? divisor : b_norm;

        return DivideNormalizedWithReciprocal(
            a,
            divisor,
            base,
            shift,
            a_norm,
            normalizedDivisor,
            reciprocal,
            computeRemainder);
      }

      void DivideAndRemainderInto(
          std::vector<DataT> const &a,
          std::vector<DataT> &q,
          std::vector<DataT> &r,
          bool computeRemainder = true) const
      {
        auto qr = DivideAndRemainder(a, computeRemainder);
        q = std::move(qr.first);
        r = std::move(qr.second);
      }

      std::vector<DataT> Divide(std::vector<DataT> const &a) const
      {
        return DivideAndRemainder(a, false).first;
      }

      std::vector<DataT> const &Divisor() const
      {
        return divisor;
      }
    };

    static std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainder(
        std::vector<DataT> const &a,
        std::vector<DataT> const &b,
        BaseT base,
        bool computeRemainder = true)
    {
      return DivideAndRemainderLimbs(a, b, base, computeRemainder);
    }

    static std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainderLimbs(
        std::span<const DataT> a, std::span<const DataT> b, BaseT base, bool computeRemainder)
    {
      if (IsZero(b))
        throw std::invalid_argument("Division by zero");
      if (IsZero(a))
        return {std::vector<DataT>{0}, computeRemainder ? std::vector<DataT>{0} : std::vector<DataT>()};

      Int cmp = Compare(a, b);
      if (cmp < 0)
        return {std::vector<DataT>{0}, computeRemainder ? std::vector<DataT>(a.begin(), a.end())
                                                      : std::vector<DataT>()};
      if (cmp == 0)
        return {std::vector<DataT>{1}, computeRemainder ? std::vector<DataT>{0} : std::vector<DataT>()};

      // Newton supports Base2_32 and Base2_64 multi-limb divisors. Other bases
      // and single-limb divisors fall back to FastDivision.
      if ((base != Base2_32 && base != Base2_64) || b.size() <= 1)
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      // Normalize: shift so top bit of b's top limb is set.
      a = a.first(LimbDetail::TrimmedLimbCount(a));
      b = b.first(LimbDetail::TrimmedLimbCount(b));
      const Int limbBits = LimbBitsFor(base);
      Int shift = NormalizationShiftBits(b.back(), limbBits);

      std::vector<DataT> shiftedA, shiftedB;
      const auto a_norm = Normalize(a, shift, limbBits, shiftedA);
      const auto b_norm = Normalize(b, shift, limbBits, shiftedB);

      SizeT n = (SizeT)b_norm.size();
      SizeT na = (SizeT)a_norm.size();

      // The reciprocal is good to a few ulps at every size (guarded precision
      // schedule in ApproxReciprocal), so chunks of up to 2n limbs keep the
      // quotient estimate within the fixup budget without an extra full-size
      // refinement iteration. The legacy 32-in-64 layout still has the old
      // exact-doubling tower and keeps the refinement for na >= 2n.
#if BIGMATH_LIMB_64
      bool need_high_precision = false;
      (void)na;
#else
      bool need_high_precision = (na >= 2 * n);
#endif

      std::vector<DataT> R = ApproxReciprocalLimbs(b_norm, need_high_precision);
      return DivideNormalizedWithReciprocal(a, b, base, shift, a_norm, b_norm, R, computeRemainder);
    }

    static std::vector<DataT> Divide(std::vector<DataT> const &a, std::vector<DataT> const &b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }
  };
}

#endif
