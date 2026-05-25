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
  // such that R * D ≈ B^(2n). Then quotient Q ≈ (A * R) >> 2n with a small ±k
  // fixup loop. Asymptotic cost O(M(n)), one log factor better than Burnikel-Ziegler.
  //
  // Current scope: handles na ≤ 2n (square-ish divisions). For na > 2n falls back
  // to FastDivision so the dispatcher in algorithms/Division.h doesn't need extra logic.
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
    static vector<DataT> ApproxReciprocal(vector<DataT> const &D)
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
      while (cur_n < n)
      {
        SizeT new_n = std::min((SizeT)(2 * cur_n), n);

        // Pad R as seed at new_n precision: shift up by (new_n - cur_n) limbs.
        SizeT extend = new_n - cur_n;
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

      // Single-reciprocal Newton handles na ≤ 2n only.
      // For larger na fall back to FastDivision (which the dispatcher would have picked anyway
      // for non-skewed cases) — keeps the algorithm self-contained at this stage.
      if (na > 2 * n)
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      // Approximate reciprocal of normalized divisor.
      vector<DataT> R = ApproxReciprocal(b_norm);

      // Q ≈ (a_norm * R) >> (2n limbs) — top (na - n + 1)-ish limbs of the product.
      vector<DataT> AR = Multiply(a_norm, R, Base2_32);
      vector<DataT> Q;
      if (AR.size() > 2 * n)
        Q.assign(AR.begin() + 2 * n, AR.end());
      else
        Q = vector<DataT>{0};
      TrimZeros(Q);
      if (Q.empty())
        Q.push_back(0);

      // QB = Q * b_norm. Then rem = a_norm - QB with ±k fixup.
      vector<DataT> QB = Multiply(Q, b_norm, Base2_32);

      // Decrement Q while QB > a_norm. Typically 0-1 iters (Newton precision gives off-by-≤2).
      vector<DataT> one{1};
      while (Compare(QB, a_norm) > 0)
      {
        Q = Subtract(Q, one, Base2_32);
        QB = Subtract(QB, b_norm, Base2_32);
      }
      vector<DataT> rem = Subtract(a_norm, QB, Base2_32);

      // Increment Q while rem >= b_norm. Typically 0-2 iters.
      while (Compare(rem, b_norm) >= 0)
      {
        rem = Subtract(rem, b_norm, Base2_32);
        Q = Add(Q, one, Base2_32);
      }

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
