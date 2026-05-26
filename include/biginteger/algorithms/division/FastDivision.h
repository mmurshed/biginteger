#ifndef FAST_DIVISION
#define FAST_DIVISION

#include <algorithm>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;

#include "../../common/Comparator.h"
#include "../../common/Util.h"
#include "../multiplication/ClassicMultiplication.h"
#include "ClassicDivision.h"

// Möller-Granlund 3-by-2 preinverted-reciprocal qhat in the Knuth Algorithm D
// inner loop, replacing the per-iteration 128/64 division. Default on; set to 0
// for A/B comparison or emergency rollback.
#ifndef BIGMATH_FASTDIV_USE_MG_QHAT
#define BIGMATH_FASTDIV_USE_MG_QHAT 1
#endif

namespace BigMath
{
  class FastDivision
  {
  private:
    static DataT Digit(ULong128 value, BaseT base)
    {
      return base == Base2_32 ? (DataT)(value & 0xFFFFFFFFULL) : (DataT)(value % base);
    }

    static ULong Carry(ULong128 value, BaseT base)
    {
      return base == Base2_32 ? (ULong)(value >> 32) : (ULong)(value / base);
    }

    static bool SubtractMul(vector<DataT> &u, vector<DataT> const &v, ULong qhat, SizeT j, BaseT base)
    {
      ULong borrow = 0;
      SizeT n = (SizeT)v.size();

      for (SizeT i = 0; i < n; ++i)
      {
        ULong128 product = (ULong128)v[i] * qhat + borrow;
        DataT low = Digit(product, base);
        borrow = Carry(product, base);

        if (u[j + i] < low)
        {
          u[j + i] = (DataT)((ULong)u[j + i] + (ULong)base - low);
          ++borrow;
        }
        else
        {
          u[j + i] = (DataT)(u[j + i] - low);
        }
      }

      if (u[j + n] < borrow)
      {
        u[j + n] = (DataT)((ULong)u[j + n] + (ULong)base - borrow);
        return true;
      }

      u[j + n] = (DataT)(u[j + n] - borrow);
      return false;
    }

    static void AddBack(vector<DataT> &u, vector<DataT> const &v, SizeT j, BaseT base)
    {
      ULong carry = 0;
      SizeT n = (SizeT)v.size();

      for (SizeT i = 0; i < n; ++i)
      {
        ULong sum = (ULong)u[j + i] + v[i] + carry;
        u[j + i] = Digit(sum, base);
        carry = Carry(sum, base);
      }

      u[j + n] = Digit((ULong)u[j + n] + carry, base);
    }

    // Möller-Granlund 3-by-2 reciprocal for B = 2^32. Returns
    //   v = floor((B^3 - 1) / d) - B,  where d = d1*B + d0, d1's top bit set.
    // Computed once per FastDivision call, amortized over the j-loop.
    static inline DataT Reciprocal3by2_Base32(DataT d1, DataT d0)
    {
      ULong d = ((ULong)d1 << 32) | d0;
      ULong128 num = ((ULong128)1 << 96) - 1; // 2^96 - 1
      return (DataT)((ULong)(num / d) - ((ULong)1 << 32));
    }

    // M-G Algorithm 4 (3/2 division) for B = 2^32.
    // Returns qhat = quotient digit, exact or +1 vs true digit. The caller's
    // SubtractMul+AddBack handles the +1 case unchanged.
    static inline DataT MGQhat_Base32(DataT u2, DataT u1, DataT u0,
                                       DataT d1, DataT d0, DataT v)
    {
      // (q1, q0) = v*u2 + (u2:u1), discarding any carry-out (2-word view).
      ULong q = (ULong)u2 * v + (((ULong)u2 << 32) | u1);
      DataT q1 = (DataT)(q >> 32);
      DataT q0 = (DataT)(q & 0xFFFFFFFFULL);

      // 2-word remainder via modular subtraction: r = (u1:u0) - q1*(d1:d0) mod 2^64.
      DataT r1_tmp = (DataT)((ULong)u1 - (ULong)q1 * d1);
      ULong r = ((ULong)r1_tmp << 32) | u0;
      ULong d_pair = ((ULong)d1 << 32) | d0;
      r -= d_pair;
      r -= (ULong)q1 * d0;
      ++q1;

      // First correction: if r-high >= saved q0, q1 was over-bumped; undo.
      if ((DataT)(r >> 32) >= q0)
      {
        --q1;
        r += d_pair; // 64-bit add; discard carry-out
      }

      // Second correction: residual r >= d.
      if (r >= d_pair)
        ++q1;

      return q1;
    }

    static vector<DataT> MultiplyByScalar(span<const DataT> a, DataT d, BaseT base)
    {
      if (d == 0 || IsZero(a))
        return vector<DataT>{0};

      vector<DataT> out(a.size() + 1, 0);
      if (base == Base2_32)
      {
        ULong128 carry = 0;
        for (SizeT i = 0; i < a.size(); ++i)
        {
          ULong128 product = (ULong128)a[i] * d + carry;
          out[i] = (DataT)(product & 0xFFFFFFFFULL);
          carry = product >> 32;
        }
        out[a.size()] = (DataT)carry;
      }
      else
      {
        ULong carry = 0;
        for (SizeT i = 0; i < a.size(); ++i)
        {
          ULong product = (ULong)a[i] * d + carry;
          out[i] = (DataT)(product % base);
          carry = product / base;
        }
        out[a.size()] = (DataT)carry;
      }

      TrimZerosToOne(out);
      return out;
    }

    static vector<DataT> DivideByScalar(span<const DataT> a, DataT d, BaseT base, DataT *remainder = nullptr)
    {
      if (d == 0)
        throw invalid_argument("Division by zero");

      vector<DataT> q(a.size(), 0);

      if (base == Base2_32)
      {
        ULong rem = 0;
        for (Int i = (Int)a.size() - 1; i >= 0; --i)
        {
          ULong128 cur = ((ULong128)rem << 32) | a[i];
          q[i] = (DataT)(cur / d);
          rem = (ULong)(cur % d);
        }
        if (remainder)
          *remainder = (DataT)rem;
      }
      else
      {
        ULong rem = 0;
        for (Int i = (Int)a.size() - 1; i >= 0; --i)
        {
          ULong cur = rem * (ULong)base + a[i];
          q[i] = (DataT)(cur / d);
          rem = cur % d;
        }
        if (remainder)
          *remainder = (DataT)rem;
      }

      TrimZerosToOne(q);
      return q;
    }

  public:
    // Span-based primary entry. Zero-copy for read-only inputs; internal u/v/q/r vectors hold mutable state.
    static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
        span<const DataT> a,
        span<const DataT> b,
        BaseT base,
        bool computeRemainder = true)
    {
      if (IsZero(b))
        throw invalid_argument("Division by zero");

      if (IsZero(a))
        return {vector<DataT>{0}, vector<DataT>{0}};

      Int cmp = Compare(a, b);
      if (cmp < 0)
        return {vector<DataT>{0}, computeRemainder ? vector<DataT>(a.begin(), a.end()) : vector<DataT>()};
      if (cmp == 0)
        return {vector<DataT>{1}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

      if (b.size() == 1)
      {
        DataT rem = 0;
        vector<DataT> q = DivideByScalar(a, b[0], base, computeRemainder ? &rem : nullptr);
        return {q, computeRemainder ? vector<DataT>{rem} : vector<DataT>()};
      }

      DataT d = (DataT)(base / (b[b.size() - 1] + 1));

      vector<DataT> u = d > 1
                            ? MultiplyByScalar(a, d, base)
                            : vector<DataT>(a.begin(), a.end());
      vector<DataT> v = d > 1
                            ? MultiplyByScalar(b, d, base)
                            : vector<DataT>(b.begin(), b.end());

      TrimZeros(u);
      TrimZeros(v);
      SizeT n = (SizeT)v.size();
      SizeT m = (SizeT)(u.size() - n);
      u.push_back(0);

      vector<DataT> q(m + 1, 0);

      bool useMG = false;
      DataT mg_v = 0;
#if BIGMATH_FASTDIV_USE_MG_QHAT
      if (base == Base2_32 && n >= 2)
      {
        useMG = true;
        mg_v = Reciprocal3by2_Base32((DataT)v[n - 1], (DataT)v[n - 2]);
      }
#endif

      for (Int j = (Int)m; j >= 0; --j)
      {
        ULong qhat;
        if (useMG)
        {
          qhat = (ULong)MGQhat_Base32(
              (DataT)u[j + n], (DataT)u[j + n - 1], (DataT)u[j + n - 2],
              (DataT)v[n - 1], (DataT)v[n - 2], mg_v);
        }
        else
        {
          ULong128 numerator = (ULong128)u[j + n] * (ULong)base + u[j + n - 1];
          qhat = (ULong)(numerator / v[n - 1]);
          ULong rhat = (ULong)(numerator % v[n - 1]);

          while (qhat == (ULong)base ||
                 (n > 1 && (ULong128)qhat * v[n - 2] >
                               (ULong128)rhat * (ULong)base + u[j + n - 2]))
          {
            --qhat;
            rhat += v[n - 1];
            if (rhat >= (ULong)base)
              break;
          }
        }

        if (SubtractMul(u, v, qhat, (SizeT)j, base))
        {
          --qhat;
          AddBack(u, v, (SizeT)j, base);
        }

        q[j] = (DataT)qhat;
      }

      TrimZerosToOne(q);

      vector<DataT> r;
      if (computeRemainder)
      {
        r.assign(u.begin(), u.begin() + n);
        TrimZerosToOne(r);
        if (d > 1)
          r = DivideByScalar(span<const DataT>(r), d, base);
        TrimZerosToOne(r);
      }

      return {q, r};
    }

    static vector<DataT> Divide(span<const DataT> a, span<const DataT> b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }

    // Vector overloads — thin wrappers for backward compatibility.
    static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base,
        bool computeRemainder = true)
    {
      return DivideAndRemainder(span<const DataT>(a), span<const DataT>(b), base, computeRemainder);
    }

    static vector<DataT> Divide(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }
  };
}

#endif
