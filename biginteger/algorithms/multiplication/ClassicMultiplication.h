/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef CLASSIC_MULTIPLICATION
#define CLASSIC_MULTIPLICATION

#include <vector>
#include <string>
using namespace std;

#include "../../common/Util.h"

namespace BigMath
{
  // All operations are unsigned
  class ClassicMultiplication
  {
  private:
    static DataT LowDigit(ULong value, BaseT base)
    {
      return base == Base2_32 ? (DataT)(value & 0xFFFFFFFFULL) : (DataT)(value % base);
    }

    static ULong NextCarry(ULong value, BaseT base)
    {
      return base == Base2_32 ? value >> 32 : value / base;
    }

    static void SetOrPush(vector<DataT> &a, SizeT pos, DataT value)
    {
      if (pos < a.size())
        a.at(pos) = value;
      else
        a.push_back(value);
    }

  public:
    static void MultiplyTo(
        vector<DataT> &a,
        DataT b,
        BaseT base)
    {
      MultiplyTo(
          a, 0, a.size() - 1,
          b,
          base);
    }

    static SizeT MultiplyTo(
        vector<DataT> &a, SizeT aStart, SizeT aEnd,
        DataT b,
        BaseT base)
    {
      if (b == 0 ||                // a times 0
          IsZero(a, aStart, aEnd)) // 0 times b
      {
        SetBit(a, aStart, aEnd, 0);
        return 0;
      }

      // Base-2^32 fast path: __int128 carry handles full 64-bit b without overflow,
      // shift/mask replaces runtime mod/div.
      if (base == Base2_32)
      {
        unsigned __int128 carry = 0;
        for (SizeT j = aStart; j <= aEnd; ++j)
        {
          unsigned __int128 p = (unsigned __int128)a.at(j) * b + carry;
          SetOrPush(a, j, (DataT)(p & 0xFFFFFFFFULL));
          carry = p >> 32;
        }
        SizeT j = aEnd + 1;
        while (carry)
        {
          SetOrPush(a, j, (DataT)(carry & 0xFFFFFFFFULL));
          carry >>= 32;
          ++j;
        }
        return j;
      }

      ULong carry = 0;

      for (Int j = aStart; j <= aEnd; j++)
      {
        ULong multiply = a.at(j);
        multiply *= b;
        multiply += carry;

        SetOrPush(a, j, LowDigit(multiply, base));
        carry = NextCarry(multiply, base);
      }

      Int j = aEnd + 1;
      while (carry > 0)
      {
        SetOrPush(a, j, LowDigit(carry, base));
        carry = NextCarry(carry, base);
        j++;
      }

      return j;
    }

    static vector<DataT> Multiply(
        vector<DataT> const &a,
        DataT b,
        BaseT base)
    {
      vector<DataT> w(a.size());
      Multiply(
          a, 0, a.size() - 1,
          b,
          w, 0, w.size() - 1,
          base);
      TrimZeros(w);
      return w;
    }

    static SizeT Multiply(
        vector<DataT> const &a, SizeT aStart, SizeT aEnd,
        ULong b,
        vector<DataT> &w, SizeT wStart, SizeT wEnd,
        BaseT base)
    {
      if (b == 0 ||                // a times 0
          IsZero(a, aStart, aEnd)) // 0 times b
      {
        SetBit(w, wStart, wEnd, 0);
        return 0;
      }

      SizeT len = Len(aStart, aEnd);
      SizeT wPos = wStart;

      ULong carry = 0;

      for (Int j = 0; j < len; j++)
      {
        ULong multiply = 0;
        SizeT aPos = aStart + j;
        if (aPos <= aEnd)
          multiply = a.at(aPos);
        multiply *= b;
        multiply += carry;

        wPos = wStart + j;
        SetOrPush(w, wPos, LowDigit(multiply, base));
        carry = NextCarry(multiply, base);
      }

      wPos = wStart + len;
      while (carry > 0)
      {
        SetOrPush(w, wPos, LowDigit(carry, base));
        carry = NextCarry(carry, base);
        wPos++;
      }

      return wPos;
    }

  public:
    static void MultiplyTo(
        vector<DataT> &a,
        vector<DataT> const &b,
        BaseT base)
    {
      a = Multiply(a, b, base);
    }

    static vector<DataT> Multiply(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base)
    {
      if (IsZero(a) || IsZero(b)) // 0 times
        return vector<DataT>();

        // If b is a single digit, use the scalar multiplication
      if (b.size() == 1)
        return Multiply(a, b[0], base);
      // If a is a single digit, use the scalar multiplication
      if (a.size() == 1)
        return Multiply(b, a[0], base);

      SizeT size = (SizeT)(a.size() + b.size() + 1);
      vector<DataT> result(size);

      Multiply(a, 0, (SizeT)a.size() - 1, b, 0, (SizeT)b.size() - 1, result, 0, base);

      TrimZeros(result);
      return result;
    }

    // Implentation of multiplication by paper-pencil method
    // Classical algorithm
    // Runtime O(n^2), Space O(n)
    static SizeT Multiply(
        vector<DataT> const &a, SizeT aStart, SizeT aEnd,
        vector<DataT> const &b, SizeT bStart, SizeT bEnd,
        vector<DataT> &result, SizeT rStart,
        BaseT base)
    {
      aEnd = min(aEnd, (SizeT)(a.size() - 1));
      bEnd = min(bEnd, (SizeT)(b.size() - 1));

      SizeT k = rStart;
      SizeT lenA = aEnd - aStart + 1;

      for (SizeT i = bStart; i <= bEnd; i++)
      {
        ULong carry = 0;
        SizeT jStart = rStart + (i - bStart);
        for (SizeT j = aStart; j <= aEnd; j++)
        {
          SizeT k = jStart + (j - aStart);
          ULong multiply = a.at(j);
          multiply *= b.at(i);
          multiply += result.at(k);
          multiply += carry;

          result.at(k) = LowDigit(multiply, base);
          carry = NextCarry(multiply, base);
        }
        k = jStart + lenA;
        result.at(k) = (DataT)carry;
      }
      return k;
    }
  };
}

#endif
