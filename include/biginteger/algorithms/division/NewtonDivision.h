#ifndef NEWTON_DIVISION
#define NEWTON_DIVISION

#include <algorithm>
#include <bit>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;

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
  private:
    struct ScratchBuffers
    {
      vector<DataT> v0;
      vector<DataT> v1;
      vector<DataT> v2;
      vector<DataT> v3;
      vector<DataT> v4;
      vector<DataT> v5;
      vector<DataT> v6;
      vector<vector<DataT>> qPieces;
    };

    static ScratchBuffers &Scratch()
    {
      static thread_local ScratchBuffers scratch;
      return scratch;
    }

    // Bit-shift left by `bits` in [0, LimbBits-1].
    static vector<DataT> ShiftLeftBits(vector<DataT> const &v, Int bits)
    {
      if (bits == 0 || IsZero(v))
        return v;
      vector<DataT> out(v.size() + 1, 0);
      for (SizeT i = 0; i < v.size(); ++i)
      {
#if BIGMATH_LIMB_64
        ULong128 cur = (ULong128)v[i] << bits;
        out[i] |= (DataT)(cur & 0xFFFFFFFFFFFFFFFFULL);
        out[i + 1] |= (DataT)(cur >> 64);
#else
        ULong cur = (ULong)v[i] << bits;
        out[i] |= (DataT)(cur & 0xFFFFFFFFULL);
        out[i + 1] |= (DataT)(cur >> 32);
#endif
      }
      TrimZerosToOne(out);
      return out;
    }

    // Bit-shift right by `bits` in [0, LimbBits-1].
    static vector<DataT> ShiftRightBits(vector<DataT> const &v, Int bits)
    {
      if (bits == 0 || IsZero(v))
        return v;
      vector<DataT> out(v.size(), 0);
      for (Int i = 0; i < (Int)v.size(); ++i)
      {
#if BIGMATH_LIMB_64
        ULong128 cur = (ULong128)v[i];
        if (i + 1 < (Int)v.size())
          cur |= ((ULong128)v[i + 1]) << 64;
        out[i] = (DataT)((cur >> bits) & 0xFFFFFFFFFFFFFFFFULL);
#else
        ULong cur = v[i];
        if (i + 1 < (Int)v.size())
          cur |= ((ULong)v[i + 1]) << 32;
        out[i] = (DataT)((cur >> bits) & 0xFFFFFFFFULL);
#endif
      }
      TrimZerosToOne(out);
      return out;
    }

    // ApproxReciprocal: given n-limb normalized D (top bit of D[n-1] set),
    // returns R such that R*D ≈ B^(2n), off by at most a small constant.
    // R has up to n+1 limbs. Implementation: 2-limb hardware seed, then
    // precision-doubling Newton iteration R_new = R * (2S - D*R) / S.
    //
    // `high_precision`: when true, run one extra Newton iter at full precision after main
    // convergence. Doubles the correct-bit count (Newton is quadratic) — needed only when
    // the caller will use R against an `a` of size 2n+1 (the +1-limb bump band from the
    // normalize shift). Skipped otherwise to avoid a ~20% cost on the common case.
    static vector<DataT> ApproxReciprocal(vector<DataT> const &D, bool high_precision = false)
    {
      SizeT n = (SizeT)D.size();
      auto &scratch = Scratch();
      vector<DataT> &R_pad = scratch.v0;
      vector<DataT> &D_new = scratch.v1;
      vector<DataT> &T = scratch.v2;
      vector<DataT> &two_S = scratch.v3;
      vector<DataT> &diff = scratch.v4;
      vector<DataT> &RD = scratch.v5;

      // n == 1 base case: R ≈ B^2 / D[0].
      if (n == 1)
      {
#if BIGMATH_LIMB_64
        // B^2 = 2^128 doesn't fit in ULong128; use (2^128 - 1) / D[0]. Result
        // is in [2^64, 2^65 - 1] for normalized D[0] ∈ [2^63, 2^64 - 1], so up
        // to 2 limbs. Newton fixup loop corrects the off-by-one when relevant.
        ULong128 num = ~(ULong128)0;
        ULong128 R = num / D[0];
        vector<DataT> result;
        result.push_back((DataT)R);
        if ((ULong)(R >> 64))
          result.push_back((DataT)(R >> 64));
        TrimZerosToOne(result);
        return result;
#else
        ULong128 num = (ULong128)1 << 64; // B^2
        ULong128 R = num / D[0];
        vector<DataT> result;
        result.push_back((DataT)(R & 0xFFFFFFFFULL));
        result.push_back((DataT)((R >> 32) & 0xFFFFFFFFULL));
        if (R >> 64)
          result.push_back((DataT)(R >> 64));
        TrimZerosToOne(result);
        return result;
#endif
      }

#if BIGMATH_LIMB_64
      // 64-bit limbs: bootstrap from top single limb. R_seed ≈ B^2 / D[n-1]
      // gives a 1-2-limb reciprocal at cur_n = 1, then Newton's quadratic
      // doubling reaches full n-limb precision in ceil(log2(n)) iterations.
      ULong128 D_top = (ULong128)D[n - 1];
      ULong128 R_seed = (~(ULong128)0) / D_top; // floor((2^128 - 1) / D[n-1])

      vector<DataT> R;
      R.push_back((DataT)R_seed);
      if ((ULong)(R_seed >> 64))
        R.push_back((DataT)(R_seed >> 64));

      SizeT cur_n = 1;
#else
      // Bootstrap from top 2 limbs: R_seed ≈ B^4 / D_top.
      ULong128 D_top = ((ULong128)D[n - 1] << 32) | D[n - 2];
      ULong128 R_seed = ((ULong128)-1) / D_top; // floor((2^128 - 1) / D_top)

      vector<DataT> R;
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
          new_n = std::min((SizeT)(2 * cur_n), n);
          extend = new_n - cur_n;
        }
        else
        {
          new_n = n;
          extend = 0;
          ++extra_iters_done;
        }

        // Pad R as seed at new_n precision: shift up by (new_n - cur_n) limbs.
        R_pad.assign(R.size() + extend, 0);
        std::memcpy(R_pad.data() + extend, R.data(), R.size() * sizeof(DataT));

        // D_new = top new_n limbs of D.
        D_new.assign(D.begin() + (n - new_n), D.end());

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
          ULong nLin = std::bit_ceil(((ULong)D_new.size() + R_pad.size()) * c);
          SizeT L = (SizeT)(nCyc / c);
          if (D_new.size() + R_pad.size() >= BIGMATH_CYCLIC_NTT_THRESHOLD &&
              nCyc < nLin && nCyc <= (1u << 26) &&
              D_new.size() <= L && R_pad.size() <= L)
          {
            wrappedIter = true;

            vector<DataT> W = NTTMultiplication::MultiplyMod2km1(D_new, R_pad, L, CurrentBase);

            const DataT maxLimb = (CurrentBase == Base2_64) ? (DataT)~0ULL : (DataT)0xFFFFFFFFULL;
            const vector<DataT> M(L, maxLimb); // B^L − 1

            // E = (B^(2m) − W) mod M; B^(2m) ≡ B^(2m mod L).
            SizeT r_exp = (SizeT)((2 * (ULong)m) % L);
            vector<DataT> Br(r_exp + 1, 0);
            Br[r_exp] = 1;
            vector<DataT> E = (Compare(Br, W) >= 0)
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
              vector<DataT> C{0};
              if (!IsZero(E) && TrimmedSize(E) > el)
              {
                vector<DataT> Xt(R_pad.begin() + xl, R_pad.end());
                vector<DataT> Et(E.begin() + el, E.end());
                vector<DataT> P = Multiply(Xt, Et, CurrentBase);
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
                static const vector<DataT> two{2};
                vector<DataT> adj = Add(C, two, CurrentBase);
                R = (Compare(R_pad, adj) >= 0) ? Subtract(R_pad, adj, CurrentBase)
                                               : vector<DataT>{0};
              }
              TrimZerosToOne(R);
            }
          }
        }
#endif

        if (!wrappedIter)
        {
          // T = D_new * R_pad
          T = Multiply(D_new, R_pad, CurrentBase);

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
            R = vector<DataT>{0};
          }
          TrimZerosToOne(R);
        }

        cur_n = new_n;
      }

      return R;
    }

    // Sentinel for "fixup loop blew the cap" — caller should fall back to FastDivision.
    struct DivideResult
    {
      vector<DataT> q;
      vector<DataT> rem; // unshifted (still in normalized space)
      bool ok;
    };

    static SizeT TrimmedSize(vector<DataT> const &v)
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
        vector<DataT> const &chunk,
        vector<DataT> const &b_norm,
        SizeT L,
        vector<DataT> &Q,
        vector<DataT> &rem,
        int fixupLimit)
    {
      SizeT n = (SizeT)b_norm.size();
      const DataT maxLimb = (CurrentBase == Base2_64) ? (DataT)~0ULL : (DataT)0xFFFFFFFFULL;
      const vector<DataT> M(L, maxLimb); // B^L − 1

      vector<DataT> W = NTTMultiplication::MultiplyMod2km1(Q, b_norm, L, CurrentBase);

      // chunk mod M: fold the limbs above L back onto the low L (B^L ≡ 1).
      vector<DataT> cm(chunk.begin(), chunk.begin() + std::min((SizeT)chunk.size(), L));
      if (chunk.size() > L)
      {
        vector<DataT> hi(chunk.begin() + L, chunk.end());
        cm = Add(cm, hi, CurrentBase);
        while (Compare(cm, M) >= 0)
          cm = Subtract(cm, M, CurrentBase);
      }

      // rem_m = (cm − W) mod M
      vector<DataT> rem_m = (Compare(cm, W) >= 0)
                                ? Subtract(cm, W, CurrentBase)
                                : Subtract(Add(cm, M, CurrentBase), W, CurrentBase);

      static const vector<DataT> one{1};

      // t < 0 (Q overestimated): rem_m ≈ M − |t|, which needs > n+1 limbs.
      int iters = 0;
      while (TrimmedSize(rem_m) > n + 1)
      {
        if (++iters > fixupLimit)
          return false;
        Q = Subtract(Q, one, CurrentBase);
        rem_m = Add(rem_m, b_norm, CurrentBase);
        while (Compare(rem_m, M) >= 0)
          rem_m = Subtract(rem_m, M, CurrentBase);
      }

      // rem_m now equals chunk − Q·b_norm exactly.
      iters = 0;
      while (Compare(rem_m, b_norm) >= 0)
      {
        if (++iters > fixupLimit)
          return false;
        rem_m = Subtract(rem_m, b_norm, CurrentBase);
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
        vector<DataT> const &chunk,
        vector<DataT> const &b_norm,
        vector<DataT> &Q,
        vector<DataT> &rem,
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

      vector<DataT> &QB = Scratch().v2;
      QB = Multiply(Q, b_norm, CurrentBase);

      static const vector<DataT> one{1};

      int iters = 0;
      while (Compare(QB, chunk) > 0)
      {
        if (++iters > fixupLimit)
          return false;
        Q = Subtract(Q, one, CurrentBase);
        QB = Subtract(QB, b_norm, CurrentBase);
      }
      rem = Subtract(chunk, QB, CurrentBase);

      iters = 0;
      while (Compare(rem, b_norm) >= 0)
      {
        if (++iters > fixupLimit)
          return false;
        rem = Subtract(rem, b_norm, CurrentBase);
        Q = Add(Q, one, CurrentBase);
      }

      return true;
    }

    // Reciprocal-based divide of a single chunk by b_norm. Used by both single-block and
    // blockwise paths. Returns {q, rem, true} on success, {{}, {}, false} on fixup overflow.
    static DivideResult DivideChunk(
        vector<DataT> const &chunk,
        vector<DataT> const &b_norm,
        vector<DataT> const &R)
    {
      SizeT n = (SizeT)b_norm.size();
      auto &scratch = Scratch();
      vector<DataT> &CR = scratch.v0;
      vector<DataT> &Q = scratch.v1;

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
        tryApprox = (chunk.size() + R.size() >= NTT_MULTIPLICATION_THRESHOLD) &&
                    nTop < nFull;
      }
#endif

      for (int attempt = tryApprox ? 0 : 1; attempt < 2; ++attempt)
      {
        SizeT drop = 2 * n;
        if (attempt == 0)
        {
          SizeT sh = (SizeT)chunk.size() - (n + 1);
          vector<DataT> c_top(chunk.begin() + sh, chunk.end());
          CR = Multiply(c_top, R, CurrentBase);
          drop = 2 * n - sh;
        }
        else
        {
          // Q ≈ (chunk * R) >> (2n limbs)
          CR = Multiply(chunk, R, CurrentBase);
        }

        if (CR.size() > drop)
          Q.assign(CR.begin() + drop, CR.end());
        else
          Q.assign(1, 0);
        TrimZerosToOne(Q);

        vector<DataT> rem;
        if (RemainderAndFixups(chunk, b_norm, Q, rem, attempt == 0 ? 12 : 8))
        {
          TrimZerosToOne(Q);
          TrimZerosToOne(rem);
          return {Q, rem, true};
        }
        // Approx attempt blew its budget — retry once with the exact product.
      }

      return {{}, {}, false};
    }

    static pair<vector<DataT>, vector<DataT>> DivideNormalizedWithReciprocal(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base,
        Int shift,
        vector<DataT> const &a_norm,
        vector<DataT> const &b_norm,
        vector<DataT> const &R,
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
      bool blockwise = (na > 2 * n);
      auto &scratch = Scratch();
      vector<vector<DataT>> &q_pieces = scratch.qPieces;

      if (!blockwise)
      {
        DivideResult res = DivideChunk(a_norm, b_norm, R);
        if (!res.ok)
          return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

        vector<DataT> rem_final;
        if (computeRemainder)
        {
          rem_final = (shift > 0) ? ShiftRightBits(res.rem, shift) : res.rem;
          TrimZerosToOne(rem_final);
        }
        return {res.q, rem_final};
      }

      // ----- Blockwise path -----
      // first_chunk_size in [n+1, 2n], chosen so remaining (na - first_chunk_size) is divisible by n.
      // Formula: ((na - 1) mod n) + 1 + n.  Verified for na in {2n+2, 3n, 3n+1, 4n, …}.
      SizeT first_chunk_size = (SizeT)(((na - 1) % n) + 1 + n);
      SizeT pos_low = na - first_chunk_size;
      vector<DataT> chunk(a_norm.begin() + pos_low, a_norm.end());

      // Collect Q pieces top-down; reverse-concat at end.
      q_pieces.clear();
      q_pieces.reserve((na + n - 1) / n);
      vector<DataT> rem;
      bool is_first = true;

      while (true)
      {
        DivideResult res = DivideChunk(chunk, b_norm, R);
        if (!res.ok)
          return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

        vector<DataT> Q_block = std::move(res.q);
        rem = std::move(res.rem);

        // Non-first blocks must contribute exactly n limbs to Q (the n limbs of A consumed).
        // The rem-strict invariant (rem < b_norm) guarantees Q_block.size() ≤ n there.
        if (!is_first)
        {
          while (Q_block.size() < n)
            Q_block.push_back(0);
          if (Q_block.size() > n)
          {
            // Should be unreachable given the invariant — bail defensively.
            return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);
          }
        }
        q_pieces.push_back(std::move(Q_block));
        is_first = false;

        if (pos_low == 0)
          break;

        // Build next chunk: high = rem (≤ n limbs), low = a_norm[pos_low - n .. pos_low - 1].
        SizeT block_n = n; // by construction, remaining is a multiple of n.
        SizeT next_chunk_size = block_n + (SizeT)rem.size();
        vector<DataT> &next_chunk = scratch.v3;
        next_chunk.assign(next_chunk_size, 0);
        std::memcpy(next_chunk.data(), a_norm.data() + pos_low - block_n, block_n * sizeof(DataT));
        std::memcpy(next_chunk.data() + block_n, rem.data(), rem.size() * sizeof(DataT));
        chunk = std::move(next_chunk);
        pos_low -= block_n;
      }

      // Concat q_pieces — first piece (top) goes to high end of final Q.
      SizeT total_q = 0;
      for (auto &p : q_pieces)
        total_q += (SizeT)p.size();
      vector<DataT> Q;
      Q.reserve(total_q);
      for (auto it = q_pieces.rbegin(); it != q_pieces.rend(); ++it)
        Q.insert(Q.end(), it->begin(), it->end());
      TrimZerosToOne(Q);

      vector<DataT> rem_final;
      if (computeRemainder)
      {
        rem_final = (shift > 0) ? ShiftRightBits(rem, shift) : rem;
        TrimZerosToOne(rem_final);
      }

      return {Q, rem_final};
    }

  public:
    class Divider
    {
    private:
      vector<DataT> divisor;
      vector<DataT> b_norm;
      vector<DataT> reciprocal;
      BaseT base;
      Int shift;
      bool can_use_newton;

    public:
      Divider(vector<DataT> const &b, BaseT radix) : divisor(b), base(radix), shift(0), can_use_newton(false)
      {
        TrimZeros(divisor);
        if (IsZero(divisor))
          throw invalid_argument("Division by zero");

        if ((base != Base2_32 && base != Base2_64) || divisor.size() <= 1)
          return;

        DataT b_top = divisor.back();
#if BIGMATH_LIMB_64
        const DataT topBitMask = 0x8000000000000000ULL;
#else
        const DataT topBitMask = 0x80000000U;
#endif
        while ((b_top & topBitMask) == 0)
        {
          b_top <<= 1;
          ++shift;
        }

        b_norm = (shift > 0) ? ShiftLeftBits(divisor, shift) : divisor;
        TrimZeros(b_norm);

        // Precompute the high-precision reciprocal once. It is valid for both single-block
        // and blockwise division and avoids recomputing the expensive Newton setup for
        // repeated divisions by the same divisor.
        reciprocal = ApproxReciprocal(b_norm, true);
        can_use_newton = true;
      }

      pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
          vector<DataT> const &a,
          bool computeRemainder = true) const
      {
        if (IsZero(a))
          return {vector<DataT>{0}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

        Int cmp = Compare(a, divisor);
        if (cmp < 0)
          return {vector<DataT>{0}, computeRemainder ? a : vector<DataT>()};
        if (cmp == 0)
          return {vector<DataT>{1}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

        if (!can_use_newton)
          return FastDivision::DivideAndRemainder(a, divisor, base, computeRemainder);

        vector<DataT> a_norm = (shift > 0) ? ShiftLeftBits(a, shift) : a;
        TrimZeros(a_norm);

        return DivideNormalizedWithReciprocal(
            a,
            divisor,
            base,
            shift,
            a_norm,
            b_norm,
            reciprocal,
            computeRemainder);
      }

      void DivideAndRemainderInto(
          vector<DataT> const &a,
          vector<DataT> &q,
          vector<DataT> &r,
          bool computeRemainder = true) const
      {
        auto qr = DivideAndRemainder(a, computeRemainder);
        q = std::move(qr.first);
        r = std::move(qr.second);
      }

      vector<DataT> Divide(vector<DataT> const &a) const
      {
        return DivideAndRemainder(a, false).first;
      }

      vector<DataT> const &Divisor() const
      {
        return divisor;
      }
    };

    static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base,
        bool computeRemainder = true)
    {
      if (IsZero(b))
        throw invalid_argument("Division by zero");
      if (IsZero(a))
        return {vector<DataT>{0}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

      Int cmp = Compare(a, b);
      if (cmp < 0)
        return {vector<DataT>{0}, computeRemainder ? a : vector<DataT>()};
      if (cmp == 0)
        return {vector<DataT>{1}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

      // Newton supports Base2_32 and Base2_64 multi-limb divisors. Other bases
      // and single-limb divisors fall back to FastDivision.
      if ((base != Base2_32 && base != Base2_64) || b.size() <= 1)
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      // Normalize: shift so top bit of b's top limb is set.
      DataT b_top = b.back();
      Int shift = 0;
#if BIGMATH_LIMB_64
      const DataT topBitMask = 0x8000000000000000ULL;
#else
      const DataT topBitMask = 0x80000000U;
#endif
      while ((b_top & topBitMask) == 0)
      {
        b_top <<= 1;
        ++shift;
      }

      vector<DataT> a_norm = (shift > 0) ? ShiftLeftBits(a, shift) : a;
      vector<DataT> b_norm = (shift > 0) ? ShiftLeftBits(b, shift) : b;
      TrimZeros(a_norm);
      TrimZeros(b_norm);

      SizeT n = (SizeT)b_norm.size();
      SizeT na = (SizeT)a_norm.size();

      // High-precision reciprocal needed whenever the divide will pit R against a chunk of
      // size ≥ 2n. For single-block ratio-2 (na ≥ 2n), standard Newton precision (~half-bits
      // accurate at large n due to integer rounding) makes Q_est off by O(n) — fixup loop
      // can't catch up. The extra refinement iter at full precision drops Q error to ≤ 1.
      bool need_high_precision = (na >= 2 * n);

      vector<DataT> R = ApproxReciprocal(b_norm, need_high_precision);
      return DivideNormalizedWithReciprocal(a, b, base, shift, a_norm, b_norm, R, computeRemainder);
    }

    static vector<DataT> Divide(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }
  };
}

#endif
