/**
 * BigMath: Multi-precision addition implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/Addition.h"

#include <algorithm>

namespace BigMath
{
  void AddTo(std::vector<DataT> &a, SizeT aStart, SizeT aEnd,
             ULong b, BaseT base)
  {
    if (base == Base2_32)
    {
      ULong128 acc = b;
      SizeT i = aStart;
      while (acc)
      {
        if (i <= aEnd && i < a.size())
        {
          acc += a[i];
          a[i] = (DataT)(acc & 0xFFFFFFFFULL);
        }
        else
        {
          a.push_back((DataT)(acc & 0xFFFFFFFFULL));
        }
        acc >>= 32;
        ++i;
      }
      return;
    }

    ULong sum = b;
    ULong carry = 0;
    if (aEnd >= aStart)
    {
      sum += a[aStart];
      a[aStart] = (DataT)(sum % base);
    }
    else
    {
      a.push_back((DataT)(sum % base));
    }

    carry = sum / base;

    SizeT aPos = aStart + 1;
    while (carry > 0)
    {
      sum = carry;
      if (aEnd >= aPos)
      {
        sum += a[aPos];
        a[aPos] = (DataT)(sum % base);
      }
      else
        a.push_back((DataT)(sum % base));
      carry = sum / base;
      aPos++;
    }

    while (carry > 0)
    {
      if (aEnd >= aPos)
        a[aPos] = (DataT)(carry % base);
      else
        a.push_back((DataT)(carry % base));
      carry = carry / base;
      aPos++;
    }
  }

  // Paper-pencil add: result[rStart..rStart+size-1] = a[aStart..aEnd] + b[bStart..bEnd].
  // O(n) time, O(1) extra space.
  void Add(std::vector<DataT> const &a, SizeT aStart, SizeT aEnd,
           std::vector<DataT> const &b, SizeT bStart, SizeT bEnd,
           std::vector<DataT> &result, SizeT rStart,
           BaseT base)
  {
    aEnd = std::min(aEnd, (SizeT)(a.size() - 1));
    bEnd = std::min(bEnd, (SizeT)(b.size() - 1));

    Int size = std::max(Len(aStart, aEnd), Len(bStart, bEnd));
    Long carry = 0;

    for (Int i = 0; i < size; i++)
    {
      Long digitOps = 0;

      Int aPos = i + aStart;
      if (aPos <= aEnd && aPos < (Int)a.size())
        digitOps = a[aPos];

      digitOps += carry;

      Int bPos = i + bStart;
      if (bPos <= bEnd && bPos < (Int)b.size())
        digitOps += b[bPos];

      Int rPos = rStart + i;
      if (rPos < (Int)result.size())
        result[rPos] = (DataT)(digitOps % base);
      carry = digitOps / base;
    }
    Int rPos = rStart + size;
    if (carry > 0 && rPos < (Int)result.size())
      result[rPos] += carry;
  }

  void AddTo(std::vector<DataT> &a, std::vector<DataT> const &b, BaseT base)
  {
    Add(a, 0, a.size() - 1,
        b, 0, b.size() - 1,
        a, 0,
        base);
  }

  void AddTo(std::vector<DataT> &a, SizeT aStart, SizeT aEnd,
             std::vector<DataT> const &b, SizeT bStart, SizeT bEnd,
             BaseT base)
  {
    Add(a, aStart, aEnd,
        b, bStart, bEnd,
        a, aStart,
        base);
  }

  std::vector<DataT> Add(std::vector<DataT> const &a,
                         std::vector<DataT> const &b,
                         BaseT base)
  {
    SizeT size = (SizeT)std::max(a.size(), b.size()) + 1;
    std::vector<DataT> result(size);
    result[size - 1] = 0;

    Add(a, 0, (SizeT)a.size() - 1,
        b, 0, (SizeT)b.size() - 1,
        result, 0,
        base);

    TrimZeros(result);
    return result;
  }

  void AddTo(std::vector<DataT> &a, ULong b, BaseT base)
  {
    if (a.empty())
    {
      if (b == 0)
        return;
      if (base == Base2_32)
      {
        while (b > 0)
        {
          a.push_back((DataT)(b & 0xFFFFFFFFULL));
          b >>= 32;
        }
      }
      else
      {
        while (b > 0)
        {
          a.push_back((DataT)(b % base));
          b /= base;
        }
      }
      return;
    }
    AddTo(a, 0, a.size() - 1, b, base);
  }
}
