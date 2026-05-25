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

      for (Int j = (Int)m; j >= 0; --j)
      {
        ULong128 numerator = (ULong128)u[j + n] * (ULong)base + u[j + n - 1];
        ULong qhat = (ULong)(numerator / v[n - 1]);
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
