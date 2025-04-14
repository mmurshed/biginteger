/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef SUBTRACTION
#define SUBTRACTION

#include <vector>
#include <string>
using namespace std;

#include "../common/Util.h"

namespace BigMath
{
  void SubtractFrom(
      vector<DataT> &a, SizeT aStart, SizeT aEnd,
      ULong b,
      BaseT base)
  {
    if (b == 0)
      return;

    // Subtract `b` from a[aStart]
    if (aEnd >= aStart)
    {
      if (a[aStart] >= b)
      {
        a[aStart] -= static_cast<DataT>(b);
        return; // Done — no borrow needed
      }
      else
      {
        b -= a[aStart]; // Remaining amount to subtract
        a[aStart] = static_cast<DataT>((base - b) % base);
        b = 1; // One borrow to subtract from next limb
      }
    }
    else
    {
      // a is too short, pad with zero up to aStart
      while (a.size() <= aStart)
        a.push_back(0);

      if (b > 0)
      {
        a[aStart] = static_cast<DataT>((base - b) % base);
        b = 1; // One borrow for the next
      }
    }

    // Propagate the borrow
    SizeT aPos = aStart + 1;
    while (b > 0)
    {
      if (aPos < a.size())
      {
        if (a[aPos] >= 1)
        {
          a[aPos] -= 1;
          b = 0; // Borrow resolved
        }
        else
        {
          a[aPos] = static_cast<DataT>(base - 1); // Continue the borrow
          b = 1;
        }
      }
      else
      {
        // If we need to borrow beyond current digits
        a.push_back(static_cast<DataT>(base - 1));
        b = 1;
      }

      aPos++;
    }

    // Trim any unnecessary trailing zeros (optional)
    while (!a.empty() && a.back() == 0)
      a.pop_back();
  }

  void SubtractFrom(
      vector<DataT> &a,
      ULong b,
      BaseT base)
  {
    SubtractFrom(a, 0, a.size() - 1, b, base);
  }

  // Implentation of subtraction by paper-pencil method
  // Assumption: a > b
  // Runtime O(n), Space O(n)
  void Subtract(
      vector<DataT> const &a, SizeT aStart, SizeT aEnd,
      vector<DataT> const &b, SizeT bStart, SizeT bEnd,
      vector<DataT> &result, SizeT rStart,
      BaseT base)
  {
    aEnd = min(
        aEnd,
        (SizeT)(a.size() - 1));

    bEnd = min(
        bEnd,
        (SizeT)(b.size() - 1));

    Int size = max(
        Len(aStart, aEnd),
        Len(bStart, bEnd));

    Long carry = 0;

    for (Int i = 0; i < size; i++)
    {
      Long digitOps = 0;

      Int aPos = aStart + i;
      if (aPos <= aEnd && aPos < a.size())
        digitOps = a.at(aPos);

      digitOps -= carry;

      Int bPos = bStart + i;
      if (bPos <= bEnd && bPos < b.size())
        digitOps -= b.at(bPos);

      carry = 0;
      if (digitOps < 0)
      {
        digitOps += base;
        carry = 1;
      }

      Int rPos = rStart + i;
      if (rPos < result.size())
        result.at(rPos) = (DataT)digitOps;
    }
  }

  void SubtractFrom(
      vector<DataT> &a,
      vector<DataT> const &b,
      BaseT base)
  {
    Subtract(
        a, 0, a.size() - 1,
        b, 0, b.size() - 1,
        a, 0,
        base);
  }

  void SubtractFrom(
      vector<DataT> &a, SizeT aStart, SizeT aEnd,
      vector<DataT> const &b, SizeT bStart, SizeT bEnd,
      BaseT base)
  {
    Subtract(
        a, aStart, aEnd,
        b, bStart, bEnd,
        a, aStart,
        base);
  }

  vector<DataT> Subtract(
      vector<DataT> const &a,
      vector<DataT> const &b,
      BaseT base)
  {
    SizeT size = (SizeT)max(a.size(), b.size()) + 1;
    vector<DataT> result(size);

    Subtract(
        a, 0, (SizeT)a.size() - 1,
        b, 0, (SizeT)b.size() - 1,
        result, 0, base);

    TrimZeros(result);

    return result;
  }

}

#endif
