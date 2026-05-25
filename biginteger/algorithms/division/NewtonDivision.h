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
    // Bit-shift left by `bits` in [0,31].
    static vector<DataT> ShiftLeftBits(vector<DataT> const &v, Int bits)
    {
      if (bits == 0 || IsZero(v))
        return v;
      vector<DataT> out(v.size() + 1, 0);
      for (SizeT i = 0; i < v.size(); ++i)
      {
        ULong cur = (ULong)v[i] << bits;
        out[i] |= (DataT)(cur & 0xFFFFFFFFULL);
        out[i + 1] |= (DataT)(cur >> 32);
      }
      TrimZeros(out);
      if (out.empty())
        out.push_back(0);
      return out;
    }

    // Bit-shift right by `bits` in [0,31].
    static vector<DataT> ShiftRightBits(vector<DataT> const &v, Int bits)
    {
      if (bits == 0 || IsZero(v))
        return v;
      vector<DataT> out(v.size(), 0);
      for (Int i = 0; i < (Int)v.size(); ++i)
      {
        ULong cur = v[i];
        if (i + 1 < (Int)v.size())
          cur |= ((ULong)v[i + 1]) << 32;
        out[i] = (DataT)((cur >> bits) & 0xFFFFFFFFULL);
      }
      TrimZeros(out);
      if (out.empty())
        out.push_back(0);
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

      // n == 1 base case: direct 128/32 div.
      if (n == 1)
      {
        ULong128 num = (ULong128)1 << 64; // B^2
        ULong128 R = num / D[0];
        vector<DataT> result;
        result.push_back((DataT)(R & 0xFFFFFFFFULL));
        result.push_back((DataT)((R >> 32) & 0xFFFFFFFFULL));
        if (R >> 64)
          result.push_back((DataT)(R >> 64));
        TrimZeros(result);
        if (result.empty())
          result.push_back(0);
        return result;
      }

      // Bootstrap from top 2 limbs: R_seed ≈ B^4 / D_top.
      ULong128 D_top = ((ULong128)D[n - 1] << 32) | D[n - 2];
      ULong128 R_seed = ((ULong128)-1) / D_top; // floor((2^128 - 1) / D_top)

      vector<DataT> R;
      R.push_back((DataT)(R_seed & 0xFFFFFFFFULL));
      R.push_back((DataT)((R_seed >> 32) & 0xFFFFFFFFULL));
      if (R_seed >> 64)
        R.push_back((DataT)(R_seed >> 64));

      SizeT cur_n = 2;
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
        vector<DataT> R_pad(R.size() + extend, 0);
        std::memcpy(R_pad.data() + extend, R.data(), R.size() * sizeof(DataT));

        // D_new = top new_n limbs of D.
        vector<DataT> D_new(D.begin() + (n - new_n), D.end());

        // T = D_new * R_pad
        vector<DataT> T = Multiply(D_new, R_pad, Base2_32);

        // diff = 2*B^(2*new_n) - T  (should be ≥ 0 for valid seeds).
        vector<DataT> two_S(2 * new_n + 1, 0);
        two_S[2 * new_n] = 2;

        vector<DataT> diff;
        if (Compare(two_S, T) >= 0)
        {
          diff = Subtract(two_S, T, Base2_32);
        }
        else
        {
          // R was over-estimate; shouldn't normally trigger. Clamp to 0.
          diff = vector<DataT>{0};
        }

        // R_new = (R_pad * diff) >> (32 * 2 * new_n)
        vector<DataT> RD = Multiply(R_pad, diff, Base2_32);
        if (RD.size() > 2 * new_n)
        {
          R.assign(RD.begin() + 2 * new_n, RD.end());
        }
        else
        {
          R = vector<DataT>{0};
        }
        TrimZeros(R);
        if (R.empty())
          R.push_back(0);

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

      // Q ≈ (chunk * R) >> (2n limbs)
      vector<DataT> CR = Multiply(chunk, R, Base2_32);
      vector<DataT> Q;
      if (CR.size() > 2 * n)
        Q.assign(CR.begin() + 2 * n, CR.end());
      else
        Q = vector<DataT>{0};
      TrimZeros(Q);
      if (Q.empty())
        Q.push_back(0);

      vector<DataT> QB = Multiply(Q, b_norm, Base2_32);

      const int FIXUP_LIMIT = 8;
      vector<DataT> one{1};

      int iters = 0;
      while (Compare(QB, chunk) > 0)
      {
        if (++iters > FIXUP_LIMIT)
          return {{}, {}, false};
        Q = Subtract(Q, one, Base2_32);
        QB = Subtract(QB, b_norm, Base2_32);
      }
      vector<DataT> rem = Subtract(chunk, QB, Base2_32);

      iters = 0;
      while (Compare(rem, b_norm) >= 0)
      {
        if (++iters > FIXUP_LIMIT)
          return {{}, {}, false};
        rem = Subtract(rem, b_norm, Base2_32);
        Q = Add(Q, one, Base2_32);
      }

      TrimZeros(Q);
      if (Q.empty())
        Q.push_back(0);

      return {Q, rem, true};
    }

  public:
    static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base,
        bool computeRemainder = true)
    {
      if (IsZero(b))
        throw invalid_argument("Division by zero");
      if (IsZero(a))
        return {vector<DataT>{0}, vector<DataT>{0}};

      Int cmp = Compare(a, b);
      if (cmp < 0)
        return {vector<DataT>{0}, computeRemainder ? a : vector<DataT>()};
      if (cmp == 0)
        return {vector<DataT>{1}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

      // Restricted to Base2_32 multi-limb divisor.
      if (base != Base2_32 || b.size() <= 1)
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      // Normalize: shift so top bit of b's top limb is set.
      DataT b_top = b.back();
      Int shift = 0;
      while ((b_top & 0x80000000U) == 0)
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

      // Two paths:
      //   Single-block (na ≤ 2n+1): one reciprocal-divide on the whole a.
      //   Blockwise (na > 2n+1): top chunk ∈ [n+1, 2n] limbs, then slide.
      bool blockwise = (na > 2 * n + 1);

      // High-precision reciprocal needed whenever the divide will pit R against a chunk of
      // size ≥ 2n. For single-block ratio-2 (na ≥ 2n), standard Newton precision (~half-bits
      // accurate at large n due to integer rounding) makes Q_est off by O(n) — fixup loop
      // can't catch up. The extra refinement iter at full precision drops Q error to ≤ 1.
      bool need_high_precision = (na >= 2 * n);

      vector<DataT> R = ApproxReciprocal(b_norm, need_high_precision);

      if (!blockwise)
      {
        DivideResult res = DivideChunk(a_norm, b_norm, R);
        if (!res.ok)
          return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

        vector<DataT> rem_final;
        if (computeRemainder)
        {
          rem_final = (shift > 0) ? ShiftRightBits(res.rem, shift) : res.rem;
          TrimZeros(rem_final);
          if (rem_final.empty())
            rem_final.push_back(0);
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
      vector<vector<DataT>> q_pieces;
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
        vector<DataT> next_chunk(next_chunk_size, 0);
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
      TrimZeros(Q);
      if (Q.empty())
        Q.push_back(0);

      vector<DataT> rem_final;
      if (computeRemainder)
      {
        rem_final = (shift > 0) ? ShiftRightBits(rem, shift) : rem;
        TrimZeros(rem_final);
        if (rem_final.empty())
          rem_final.push_back(0);
      }

      return {Q, rem_final};
    }

    static vector<DataT> Divide(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }
  };
}

#endif
