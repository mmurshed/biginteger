#ifndef FAST_DIVISION
#define FAST_DIVISION

#include <algorithm>
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
    static DataT Digit(unsigned __int128 value, BaseT base)
    {
      return base == Base2_32 ? (DataT)(value & 0xFFFFFFFFULL) : (DataT)(value % base);
    }

    static ULong Carry(unsigned __int128 value, BaseT base)
    {
      return base == Base2_32 ? (ULong)(value >> 32) : (ULong)(value / base);
    }

    static bool SubtractMul(vector<DataT> &u, vector<DataT> const &v, ULong qhat, SizeT j, BaseT base)
    {
      ULong borrow = 0;
      SizeT n = (SizeT)v.size();

      for (SizeT i = 0; i < n; ++i)
      {
        unsigned __int128 product = (unsigned __int128)v[i] * qhat + borrow;
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

      if (b.size() == 1)
        return ClassicDivision::DivideAndRemainder(a, b[0], base);

      SizeT n = (SizeT)b.size();
      SizeT m = (SizeT)(a.size() - n);
      DataT d = (DataT)(base / (b.back() + 1));

      vector<DataT> u = d > 1 ? ClassicMultiplication::Multiply(a, d, base) : a;
      vector<DataT> v = d > 1 ? ClassicMultiplication::Multiply(b, d, base) : b;

      TrimZeros(u);
      TrimZeros(v);
      n = (SizeT)v.size();
      m = (SizeT)(u.size() - n);
      u.push_back(0);

      vector<DataT> q(m + 1, 0);

      for (Int j = (Int)m; j >= 0; --j)
      {
        unsigned __int128 numerator = (unsigned __int128)u[j + n] * (ULong)base + u[j + n - 1];
        ULong qhat = (ULong)(numerator / v[n - 1]);
        ULong rhat = (ULong)(numerator % v[n - 1]);

        while (qhat == (ULong)base ||
               (n > 1 && (unsigned __int128)qhat * v[n - 2] >
                             (unsigned __int128)rhat * (ULong)base + u[j + n - 2]))
        {
          --qhat;
          rhat += v[n - 1];
          if (rhat >= (ULong)base)
            break;
        }

        if (SubtractMul(u, v, qhat, j, base))
        {
          --qhat;
          AddBack(u, v, j, base);
        }

        q[j] = (DataT)qhat;
      }

      TrimZeros(q);
      if (q.empty())
        q.push_back(0);

      vector<DataT> r;
      if (computeRemainder)
      {
        r.assign(u.begin(), u.begin() + n);
        TrimZeros(r);
        if (d > 1)
          r = ClassicDivision::Divide(r, d, base);
        TrimZeros(r);
        if (r.empty())
          r.push_back(0);
      }

      return {q, r};
    }

    static vector<DataT> Divide(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }
  };
}

#endif
