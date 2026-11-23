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

#include "../../common/LimbKernels.h"
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

    static void SetOrPush(std::vector<DataT> &a, SizeT pos, DataT value)
    {
      if (pos < a.size())
        a[pos] = value;
      else
        a.push_back(value);
    }

  public:
    static void MultiplyTo(
        std::vector<DataT> &a,
        DataT b,
        BaseT base)
    {
      MultiplyTo(
          a, 0, a.size() - 1,
          b,
          base);
    }

    static SizeT MultiplyTo(
        std::vector<DataT> &a, SizeT aStart, SizeT aEnd,
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

      for (SizeT j = aStart; j <= aEnd; j++)
      {
        ULong multiply = a[j];
        multiply *= b;
        multiply += carry;

        SetOrPush(a, j, LowDigit(multiply, base));
        carry = NextCarry(multiply, base);
      }

      SizeT j = aEnd + 1;
      while (carry > 0)
      {
        SetOrPush(a, j, LowDigit(carry, base));
        carry = NextCarry(carry, base);
        j++;
      }

      return j;
    }

    static std::vector<DataT> Multiply(
        std::vector<DataT> const &a, DataT b, BaseT base)
    {
      return Multiply(std::span<const DataT>(a), b, base);
    }

    // Multiply a limb view by a scalar into a fresh result.
    static std::vector<DataT> Multiply(std::span<const DataT> a, DataT b, BaseT base)
    {
      if (b == 0 || a.empty())
        return {};
      if (base == Base2_64)
      {
        const SizeT n = (SizeT)a.size();
        std::vector<DataT> w(n + 1);
        w[n] = LimbKernel::Mul1(w.data(), a.data(), n, b);
        TrimZeros(w);
        return w;
      }

      std::vector<DataT> w(a.size());
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
        std::vector<DataT> const &a, SizeT aStart, SizeT aEnd,
        ULong b,
        std::vector<DataT> &w, SizeT wStart, SizeT wEnd,
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

      for (SizeT j = 0; j < len; j++)
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
        std::vector<DataT> &a,
        std::vector<DataT> const &b,
        BaseT base)
    {
      a = Multiply(a, b, base);
    }

    // Both inputs are nonempty. The disjoint output has a.size() + b.size()
    // limbs and is fully overwritten.
    static void MultiplyInto(std::span<const DataT> a, std::span<const DataT> b,
                             DataT *out, BaseT base)
    {
      if (base == Base2_64)
      {
        LimbKernel::MulBasecase(out, a.data(), (SizeT)a.size(), b.data(), (SizeT)b.size());
        return;
      }

      std::fill_n(out, a.size() + b.size(), 0);
      AccumulateRows(a, b, out, base);
    }

    static std::vector<DataT> MultiplyLimbs(std::span<const DataT> a,
                                           std::span<const DataT> b, BaseT base)
    {
      if (IsZero(a) || IsZero(b))
        return {};
      if (b.size() == 1)
        return Multiply(a, b[0], base);
      if (a.size() == 1)
        return Multiply(b, a[0], base);

      std::vector<DataT> result(a.size() + b.size());
      MultiplyInto(a, b, result.data(), base);
      TrimZeros(result);
      return result;
    }

    static std::vector<DataT> Multiply(
        std::vector<DataT> const &a, std::vector<DataT> const &b, BaseT base)
    {
      return MultiplyLimbs(a, b, base);
    }

    // Accumulating variant for callers with an existing result range.
    static SizeT Multiply(
        std::vector<DataT> const &a, SizeT aStart, SizeT aEnd,
        std::vector<DataT> const &b, SizeT bStart, SizeT bEnd,
        std::vector<DataT> &result, SizeT rStart, BaseT base)
    {
      if (a.empty() || b.empty())
        return rStart;
      aEnd = std::min(aEnd, (SizeT)(a.size() - 1));
      bEnd = std::min(bEnd, (SizeT)(b.size() - 1));
      if (bStart > bEnd)
        return rStart;
      const SizeT la = aEnd - aStart + 1, lb = bEnd - bStart + 1;
      AccumulateRows(std::span<const DataT>(a).subspan(aStart, la),
                     std::span<const DataT>(b).subspan(bStart, lb), result.data() + rStart, base);
      return rStart + la + lb - 1;
    }

  private:
    static void AccumulateRows(std::span<const DataT> a, std::span<const DataT> b,
                               DataT *out, BaseT base)
    {
      const SizeT la = (SizeT)a.size(), lb = (SizeT)b.size();
      if (base == Base2_64)
      {
        for (SizeT i = 0; i < lb; ++i)
          out[i + la] = LimbKernel::AddMul1(out + i, a.data(), la, b[i]);
        return;
      }

      for (SizeT i = 0; i < lb; ++i)
      {
        ULong carry = 0;
        for (SizeT j = 0; j < la; ++j)
        {
          const ULong product = a[j] * b[i] + out[i + j] + carry;
          out[i + j] = LowDigit(product, base);
          carry = NextCarry(product, base);
        }
        out[i + la] = (DataT)carry;
      }
    }
  };
}

#endif
