#ifndef FAST_DIVISION
#define FAST_DIVISION

#include <algorithm>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;

#include "../../common/BitShifts.h"
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
      if (base == Base2_32) return (DataT)(value & 0xFFFFFFFFULL);
      if (base == Base2_64) return (DataT)(value & 0xFFFFFFFFFFFFFFFFULL);
      return (DataT)(value % base);
    }

    static ULong Carry(ULong128 value, BaseT base)
    {
      if (base == Base2_32) return (ULong)(value >> 32);
      if (base == Base2_64) return (ULong)(value >> 64);
      return (ULong)(value / base);
    }

    static bool SubtractMul(vector<DataT> &u, vector<DataT> const &v, ULong qhat, SizeT j, BaseT base)
    {
      ULong borrow = 0;
      SizeT n = (SizeT)v.size();

      if (base == Base2_64)
      {
        // Base2_64: the (ULong)base = 0 sentinel makes the historical
        // `u[j+i] + base - low` underflow-correction unusable. Unsigned
        // modular subtraction wraps correctly without it: when u[j+i] < low,
        // `(ULong)u[j+i] - low` = u[j+i] + 2^64 - low (mod 2^64), which is
        // exactly what the +base correction does for power-of-two B.
        for (SizeT i = 0; i < n; ++i)
        {
          ULong128 product = (ULong128)v[i] * qhat + borrow;
          DataT low = (DataT)(product & 0xFFFFFFFFFFFFFFFFULL);
          borrow = (ULong)(product >> 64);

          if (u[j + i] < low)
          {
            u[j + i] = (DataT)((ULong)u[j + i] - low); // wraps mod 2^64
            ++borrow;
          }
          else
          {
            u[j + i] = (DataT)(u[j + i] - low);
          }
        }

        if (u[j + n] < borrow)
        {
          u[j + n] = (DataT)((ULong)u[j + n] - borrow); // wraps mod 2^64
          return true;
        }
        u[j + n] = (DataT)(u[j + n] - borrow);
        return false;
      }

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
      SizeT n = (SizeT)v.size();

      if (base == Base2_64)
      {
        // ULong128 sum to capture the carry on full 64-bit limbs.
        ULong128 carry = 0;
        for (SizeT i = 0; i < n; ++i)
        {
          ULong128 sum = (ULong128)u[j + i] + v[i] + carry;
          u[j + i] = (DataT)(sum & 0xFFFFFFFFFFFFFFFFULL);
          carry = sum >> 64;
        }
        u[j + n] = (DataT)(((ULong128)u[j + n] + carry) & 0xFFFFFFFFFFFFFFFFULL);
        return;
      }

      ULong carry = 0;
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

    // Möller-Granlund 3-by-2 reciprocal for B = 2^64. Uses Algorithm 6 of
    // M-G 2010 to derive the 3/2 reciprocal from a 2/1 reciprocal of d1 —
    // avoiding the otherwise-needed 192/128 software division.
    // Precondition: d1's top bit set (divisor already normalized).
    static inline ULong Reciprocal3by2_Base64(ULong d1, ULong d0)
    {
      // 2-by-1 reciprocal of d1: v = floor((2^128 - 1) / d1) - 2^64.
      ULong v = (ULong)((~(ULong128)0) / d1);

      // Algorithm 6 (M-G 2010): refine v to 3-by-2 reciprocal of (d1, d0).
      ULong p = d1 * v + d0;
      if (p < d0)
      {
        --v;
        if (p >= d1)
        {
          --v;
          p -= d1;
        }
        p -= d1;
      }
      ULong128 t = (ULong128)v * d0;
      ULong t1 = (ULong)(t >> 64);
      ULong t0 = (ULong)t;
      p += t1;
      if (p < t1)
      {
        --v;
        if (p > d1 || (p == d1 && t0 >= d0))
          --v;
      }
      return v;
    }

    // M-G Algorithm 4 (3/2 division) for B = 2^64. Mirrors the Base32 path
    // at 128-bit widths via ULong128.
    static inline ULong MGQhat_Base64(ULong u2, ULong u1, ULong u0,
                                       ULong d1, ULong d0, ULong v)
    {
      // (q1, q0) = v * u2 + (u2, u1), discarding any 129th-bit carry-out.
      ULong128 q = (ULong128)v * u2 + (((ULong128)u2 << 64) | u1);
      ULong q1 = (ULong)(q >> 64);
      ULong q0 = (ULong)q;

      // 2-word remainder via 128-bit modular subtraction:
      //   r = (u1:u0) - q1*(d1:d0) - (d1:d0) mod 2^128.
      ULong r1_tmp = u1 - q1 * d1;
      ULong128 r = ((ULong128)r1_tmp << 64) | u0;
      ULong128 d_pair = ((ULong128)d1 << 64) | d0;
      r -= d_pair;
      r -= (ULong128)q1 * d0;
      ++q1;

      // First correction: if r-high >= saved q0, q1 was over-bumped; undo.
      if ((ULong)(r >> 64) >= q0)
      {
        --q1;
        r += d_pair; // 128-bit add; discard carry-out
      }

      // Second correction: residual r >= d.
      if (r >= d_pair)
        ++q1;

      return q1;
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
      else if (base == Base2_64)
      {
        ULong carry = 0;
        for (SizeT i = 0; i < a.size(); ++i)
        {
          ULong128 product = (ULong128)a[i] * d + carry;
          out[i] = (DataT)(product & 0xFFFFFFFFFFFFFFFFULL);
          carry = (ULong)(product >> 64);
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
      else if (base == Base2_64)
      {
        ULong rem = 0;
        for (Int i = (Int)a.size() - 1; i >= 0; --i)
        {
          ULong128 cur = ((ULong128)rem << 64) | a[i];
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

      bool pow2base = (base == Base2_32 || base == Base2_64);
      int limbBits = (base == Base2_64) ? 64 : 32;
      int shift = 0;
      DataT d = 1;
      vector<DataT> u, v;
      if (pow2base)
      {
        DataT btop = b[b.size() - 1];
        shift = (limbBits == 64) ? __builtin_clzll((unsigned long long)btop)
                                 : __builtin_clz((unsigned)btop);
        u = ShiftLeftBits(a, shift, limbBits);
        v = ShiftLeftBits(b, shift, limbBits);
      }
      else
      {
        d = (DataT)(base / (b[b.size() - 1] + 1));
        u = d > 1 ? MultiplyByScalar(a, d, base)
                  : vector<DataT>(a.begin(), a.end());
        v = d > 1 ? MultiplyByScalar(b, d, base)
                  : vector<DataT>(b.begin(), b.end());
      }

      TrimZeros(u);
      TrimZeros(v);
      SizeT n = (SizeT)v.size();
      SizeT m = (SizeT)(u.size() - n);
      u.push_back(0);

      vector<DataT> q(m + 1, 0);

      bool useMG32 = false;
      bool useMG64 = false;
      DataT mg_v32 = 0;
      ULong mg_v64 = 0;
#if BIGMATH_FASTDIV_USE_MG_QHAT
      if (n >= 2)
      {
        if (base == Base2_32)
        {
          useMG32 = true;
          mg_v32 = Reciprocal3by2_Base32((DataT)v[n - 1], (DataT)v[n - 2]);
        }
        else if (base == Base2_64)
        {
          useMG64 = true;
          mg_v64 = Reciprocal3by2_Base64((ULong)v[n - 1], (ULong)v[n - 2]);
        }
      }
#endif

      for (Int j = (Int)m; j >= 0; --j)
      {
        ULong qhat;
        // MG 3/2 precondition is (u2:u1) < (d1:d0). The Knuth loop invariant
        // (window < d*B) only gives <= when the divisor has nonzero limbs
        // below the top two; on equality MGQhat's modular arithmetic wraps
        // and returns garbage. The true digit in that case is exactly B-1:
        // (B-1)*d <= window < B*d. GMP's sbpi1_div_qr special-cases this
        // identically. (For n == 2 equality cannot occur: remainder < d.)
        bool mgTopEqual =
            (useMG32 || useMG64) &&
            u[j + n] == v[n - 1] && u[j + n - 1] == v[n - 2];
        if (mgTopEqual)
        {
          qhat = (base == Base2_64) ? ~(ULong)0 : (ULong)0xFFFFFFFFULL;
        }
        else if (useMG32)
        {
          qhat = (ULong)MGQhat_Base32(
              (DataT)u[j + n], (DataT)u[j + n - 1], (DataT)u[j + n - 2],
              (DataT)v[n - 1], (DataT)v[n - 2], mg_v32);
        }
        else if (useMG64)
        {
          qhat = MGQhat_Base64(
              (ULong)u[j + n], (ULong)u[j + n - 1], (ULong)u[j + n - 2],
              (ULong)v[n - 1], (ULong)v[n - 2], mg_v64);
        }
        else if (base == Base2_64)
        {
          // Knuth qhat for 64-bit limbs. The (u_hi:u_lo) / v_top divide is a
          // 128/64 division; rhat overflow is detected via __builtin_add_overflow
          // since B = 2^64 doesn't fit in ULong.
          ULong128 numerator = ((ULong128)u[j + n] << 64) | u[j + n - 1];
          ULong128 qhat128 = numerator / v[n - 1];

          // Cap qhat at B-1 = 2^64-1: when u[j+n] == v[n-1] (legal under
          // Algorithm D's invariant) the raw estimate is >= B and a narrowing
          // cast would silently wrap to a tiny digit. Clamp and recompute rhat
          // against the clamped digit; the correction loop and the SubtractMul
          // fixup handle the rest. This is the `qhat == base` check of the
          // generic path, done in 128 bits.
          if (qhat128 >= ((ULong128)1 << 64))
            qhat128 = (((ULong128)1 << 64) - 1);
          ULong128 rhat128 = numerator - qhat128 * v[n - 1];
          qhat = (ULong)qhat128;

          // Once rhat >= B the mul test is unconditionally false, which also
          // covers the old __builtin_add_overflow early-exit.
          while (n > 1 && rhat128 < ((ULong128)1 << 64) &&
                 (ULong128)qhat * v[n - 2] > ((rhat128 << 64) | u[j + n - 2]))
          {
            --qhat;
            rhat128 += v[n - 1];
          }
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
        if (pow2base)
          ShiftRightBitsInPlace(r, shift, limbBits);
        else if (d > 1)
        {
          TrimZerosToOne(r);
          r = DivideByScalar(span<const DataT>(r), d, base);
        }
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
