/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef CLASSIC_MULTIPLICATION
#define CLASSIC_MULTIPLICATION

#include <vector>
#include <span>
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

    // Base2_64 paths use ULong128 accumulators directly (a 64-bit * 64-bit
    // product overflows ULong), so they don't go through LowDigit/NextCarry.
    static DataT LowDigit128(ULong128 value)
    {
      return (DataT)(value & 0xFFFFFFFFFFFFFFFFULL);
    }

    static ULong NextCarry128(ULong128 value)
    {
      return (ULong)(value >> 64);
    }

    static void SetOrPush(vector<DataT> &a, SizeT pos, DataT value)
    {
      if (pos < a.size())
        a[pos] = value;
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

      // Base-2^32 fast path: ULong128 carry handles full 64-bit b without overflow,
      // shift/mask replaces runtime mod/div.
      if (base == Base2_32)
      {
        ULong128 carry = 0;
        for (SizeT j = aStart; j <= aEnd; ++j)
        {
          ULong128 p = (ULong128)a[j] * b + carry;
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

      // Base-2^64 fast path: 64-bit * 64-bit gives a 128-bit product. Carry is
      // at most 2^64-1 (the high half of the previous product), so a single
      // push suffices after the main loop.
      if (base == Base2_64)
      {
        ULong carry = 0;
        for (SizeT j = aStart; j <= aEnd; ++j)
        {
          ULong128 p = (ULong128)a[j] * b + carry;
          SetOrPush(a, j, LowDigit128(p));
          carry = NextCarry128(p);
        }
        SizeT j = aEnd + 1;
        if (carry)
        {
          SetOrPush(a, j, (DataT)carry);
          ++j;
        }
        return j;
      }

      ULong carry = 0;

      for (Int j = aStart; j <= aEnd; j++)
      {
        ULong multiply = a[j];
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

    // Span overload: multiply a span by scalar b into a fresh vector.
    // Mirrors the vector-input fast path (ULong accumulator).
    // Caller responsible for ensuring b is small enough to avoid uint64 overflow
    // in a[i] * b + carry (true for normalization scalars in FastDivision).
    static vector<DataT> Multiply(
        span<const DataT> a,
        DataT b,
        BaseT base)
    {
      if (b == 0 || a.empty())
        return vector<DataT>();
      vector<DataT> w(a.size());

      if (base == Base2_64)
      {
        ULong carry = 0;
        for (SizeT j = 0; j < a.size(); ++j)
        {
          ULong128 p = (ULong128)a[j] * b + carry;
          w[j] = LowDigit128(p);
          carry = NextCarry128(p);
        }
        if (carry)
          w.push_back((DataT)carry);
        TrimZeros(w);
        return w;
      }

      ULong carry = 0;
      for (SizeT j = 0; j < a.size(); ++j)
      {
        ULong multiply = a[j];
        multiply *= b;
        multiply += carry;
        w[j] = LowDigit(multiply, base);
        carry = NextCarry(multiply, base);
      }
      while (carry > 0)
      {
        w.push_back(LowDigit(carry, base));
        carry = NextCarry(carry, base);
      }
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

      if (base == Base2_64)
      {
        ULong carry = 0;
        for (SizeT j = 0; j < len; ++j)
        {
          ULong128 av = 0;
          SizeT aPos = aStart + j;
          if (aPos <= aEnd)
            av = a[aPos];
          ULong128 p = av * b + carry;
          wPos = wStart + j;
          SetOrPush(w, wPos, LowDigit128(p));
          carry = NextCarry128(p);
        }
        wPos = wStart + len;
        if (carry)
        {
          SetOrPush(w, wPos, (DataT)carry);
          ++wPos;
        }
        return wPos;
      }

      ULong carry = 0;

      for (Int j = 0; j < len; j++)
      {
        ULong multiply = 0;
        SizeT aPos = aStart + j;
        if (aPos <= aEnd)
          multiply = a[aPos];
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

      if (base == Base2_64)
      {
        // 64×64→128 product plus 64-bit accumulator (existing result limb) plus
        // 64-bit carry: max sum ≤ (2^64-1)² + (2^64-1) + (2^64-1) = 2^128 - 1.
        // Fits ULong128.
        for (SizeT i = bStart; i <= bEnd; i++)
        {
          ULong carry = 0;
          SizeT jStart = rStart + (i - bStart);
          for (SizeT j = aStart; j <= aEnd; j++)
          {
            SizeT kk = jStart + (j - aStart);
            ULong128 p = (ULong128)a[j] * b[i] + result[kk] + carry;
            result[kk] = LowDigit128(p);
            carry = NextCarry128(p);
          }
          k = jStart + lenA;
          result[k] = (DataT)carry;
        }
        return k;
      }

      for (SizeT i = bStart; i <= bEnd; i++)
      {
        ULong carry = 0;
        SizeT jStart = rStart + (i - bStart);
        for (SizeT j = aStart; j <= aEnd; j++)
        {
          SizeT k = jStart + (j - aStart);
          ULong multiply = a[j];
          multiply *= b[i];
          multiply += result[k];
          multiply += carry;

          result[k] = LowDigit(multiply, base);
          carry = NextCarry(multiply, base);
        }
        k = jStart + lenA;
        result[k] = (DataT)carry;
      }
      return k;
    }
  };
}

#endif
