/**
 * BigMath: Multi-precision subtraction implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/Subtraction.h"

#include <algorithm>

namespace BigMath
{
  void SubtractFrom(std::vector<DataT> &a, SizeT aStart, SizeT aEnd,
                    ULong b, BaseT base)
  {
    if (b == 0)
      return;

    if (aEnd >= aStart)
    {
      if (a[aStart] >= b)
      {
        a[aStart] -= static_cast<DataT>(b);
        return;
      }
      else
      {
        b -= a[aStart];
        a[aStart] = static_cast<DataT>((base - b) % base);
        b = 1;
      }
    }
    else
    {
      while (a.size() <= aStart)
        a.push_back(0);

      if (b > 0)
      {
        a[aStart] = static_cast<DataT>((base - b) % base);
        b = 1;
      }
    }

    SizeT aPos = aStart + 1;
    while (b > 0)
    {
      if (aPos < a.size())
      {
        if (a[aPos] >= 1)
        {
          a[aPos] -= 1;
          b = 0;
        }
        else
        {
          a[aPos] = static_cast<DataT>(base - 1);
          b = 1;
        }
      }
      else
      {
        a.push_back(static_cast<DataT>(base - 1));
        b = 1;
      }
      aPos++;
    }

    while (!a.empty() && a.back() == 0)
      a.pop_back();
  }

  void SubtractFrom(std::vector<DataT> &a, ULong b, BaseT base)
  {
    SubtractFrom(a, 0, a.size() - 1, b, base);
  }

  // Paper-pencil subtract; assumes a >= b. O(n) time, O(1) extra space.
  void Subtract(std::vector<DataT> const &a, SizeT aStart, SizeT aEnd,
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

      Int aPos = aStart + i;
      if (aPos <= aEnd && aPos < (Int)a.size())
        digitOps = a[aPos];

      digitOps -= carry;

      Int bPos = bStart + i;
      if (bPos <= bEnd && bPos < (Int)b.size())
        digitOps -= b[bPos];

      carry = 0;
      if (digitOps < 0)
      {
        digitOps += base;
        carry = 1;
      }

      Int rPos = rStart + i;
      if (rPos < (Int)result.size())
        result[rPos] = (DataT)digitOps;
    }
  }

  void SubtractFrom(std::vector<DataT> &a, std::vector<DataT> const &b, BaseT base)
  {
    Subtract(a, 0, a.size() - 1,
             b, 0, b.size() - 1,
             a, 0,
             base);
  }

  void SubtractFrom(std::vector<DataT> &a, SizeT aStart, SizeT aEnd,
                    std::vector<DataT> const &b, SizeT bStart, SizeT bEnd,
                    BaseT base)
  {
    Subtract(a, aStart, aEnd,
             b, bStart, bEnd,
             a, aStart,
             base);
  }

  std::vector<DataT> Subtract(std::vector<DataT> const &a,
                              std::vector<DataT> const &b,
                              BaseT base)
  {
    SizeT size = (SizeT)std::max(a.size(), b.size()) + 1;
    std::vector<DataT> result(size);

    Subtract(a, 0, (SizeT)a.size() - 1,
             b, 0, (SizeT)b.size() - 1,
             result, 0, base);

    TrimZeros(result);
    return result;
  }

  void Subtract(std::vector<DataT> &result,
                std::vector<DataT> const &a,
                std::vector<DataT> const &b,
                BaseT base)
  {
    SizeT need = (SizeT)std::max(a.size(), b.size()) + 1;
    if (result.size() < need)
      result.assign(need, 0);
    else
      std::fill(result.begin(), result.end(), (DataT)0);

    Subtract(a, 0, (SizeT)a.size() - 1,
             b, 0, (SizeT)b.size() - 1,
             result, 0, base);

    TrimZeros(result);
  }
}
