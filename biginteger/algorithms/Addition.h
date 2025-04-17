/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef ADDITION
#define ADDITION

#include <vector>
#include <string>
using namespace std;

#include "../common/Util.h"

namespace BigMath
{
  // All operations are unsigned
  void AddTo(
      vector<DataT> &a, SizeT aStart, SizeT aEnd,
      ULong b,
      BaseT base)
  {
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

  // Implentation of addition by paper-pencil method
  // Runtime O(n), Space O(n)
  // Adds two multi-precision numbers a and b in little-endian order.
  // The result is stored in result, starting at rStart.
  void Add(
      vector<DataT> const &a, SizeT aStart, SizeT aEnd,
      vector<DataT> const &b, SizeT bStart, SizeT bEnd,
      vector<DataT> &result, SizeT rStart,
      BaseT base)
  {
    aEnd = min(aEnd, (SizeT)(a.size() - 1));

    bEnd = min(bEnd, (SizeT)(b.size() - 1));

    Int size = max(Len(aStart, aEnd), Len(bStart, bEnd));

    Long carry = 0;

    for (Int i = 0; i < size; i++)
    {
      Long digitOps = 0;

      Int aPos = i + aStart;
      if (aPos <= aEnd && aPos < a.size())
        digitOps = a.at(aPos);

      digitOps += carry;

      Int bPos = i + bStart;
      if (bPos <= bEnd && bPos < b.size())
        digitOps += b.at(bPos);

      Int rPos = rStart + i;
      if (rPos < result.size())
        result.at(rPos) = (DataT)(digitOps % base);
      carry = digitOps / base;
    }
    Int rPos = rStart + size;
    if (carry > 0 && rPos < result.size())
      result.at(rPos) += carry;
  }

  void AddTo(
      vector<DataT> &a,
      vector<DataT> const &b,
      BaseT base)
  {

    Add(
        a, 0, a.size() - 1,
        b, 0, b.size() - 1,
        a, 0,
        base);
  }

  void AddToAt(
      vector<DataT> &a,
      vector<DataT> const &b,
      SizeT aStart,
      BaseT base)
  {

    Add(
        a, aStart, a.size() - 1,
        b, 0, b.size() - 1,
        a, aStart,
        base);
  }

  void AddTo(
      vector<DataT> &a, SizeT aStart, SizeT aEnd,
      vector<DataT> const &b, SizeT bStart, SizeT bEnd,
      BaseT base)
  {

    Add(
        a, aStart, aEnd,
        b, bStart, bEnd,
        a, aStart,
        base);
  }

  vector<DataT> Add(
      vector<DataT> const &a,
      vector<DataT> const &b,
      BaseT base)
  {
    SizeT size = (SizeT)max(a.size(), b.size()) + 1;
    vector<DataT> result(size);

    result[size - 1] = 0;

    Add(
        a, 0, (SizeT)a.size() - 1,
        b, 0, (SizeT)b.size() - 1,
        result, 0,
        base);

    TrimZeros(result);

    return result;
  }

  void AddTo(
      vector<DataT> &a,
      ULong b,
      BaseT base)
  {
    AddTo(a, 0, a.size() - 1, b, base);
  }

}

#endif
