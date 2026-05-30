#ifndef NEWTON_DIVISION
#define NEWTON_DIVISION

#include <algorithm>
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

namespace BigMath
{
  // Newton-Raphson division.
  // Precomputes an n-limb approximate reciprocal R of the normalized divisor D
  // such that R*D ≈ B^(2n). For na ≤ 2n+1: one Q = (a*R) >> 2n + small fixup.
  // For na > 2n+1: blockwise — process top first_chunk ∈ [n+1, 2n] limbs, then slide down
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
      vector<DataT> &QB = scratch.v2;

      // Q ≈ (chunk * R) >> (2n limbs)
      CR = Multiply(chunk, R, CurrentBase);
      if (CR.size() > 2 * n)
        Q.assign(CR.begin() + 2 * n, CR.end());
      else
        Q.assign(1, 0);
      TrimZerosToOne(Q);

      QB = Multiply(Q, b_norm, CurrentBase);

      const int FIXUP_LIMIT = 8;
      static const vector<DataT> one{1};

      int iters = 0;
      while (Compare(QB, chunk) > 0)
      {
        if (++iters > FIXUP_LIMIT)
          return {{}, {}, false};
        Q = Subtract(Q, one, CurrentBase);
        QB = Subtract(QB, b_norm, CurrentBase);
      }
      vector<DataT> rem = Subtract(chunk, QB, CurrentBase);

      iters = 0;
      while (Compare(rem, b_norm) >= 0)
      {
        if (++iters > FIXUP_LIMIT)
          return {{}, {}, false};
        rem = Subtract(rem, b_norm, CurrentBase);
        Q = Add(Q, one, CurrentBase);
      }

      TrimZerosToOne(Q);

      return {Q, rem, true};
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
      //   Single-block (na ≤ 2n+1): one reciprocal-divide on the whole a.
      //   Blockwise (na > 2n+1): top chunk ∈ [n+1, 2n] limbs, then slide.
      bool blockwise = (na > 2 * n + 1);
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
