/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef CLASSIC_SUBTRACTION
#define CLASSIC_SUBTRACTION

#include <vector>
#include <string>
using namespace std;

#include "../../BigIntegerUtil.h"

namespace BigMath
{
  class ClassicSubtraction
  {
  public:
    static void SubtractFrom(
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

    static vector<DataT> Subtract(
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

      BigIntegerUtil::TrimZeros(result);

      return result;
    }

    // Implentation of subtraction by paper-pencil method
    // Assumption: a > b
    // Runtime O(n), Space O(n)
    static void Subtract(
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
          BigIntegerUtil::Len(aStart, aEnd),
          BigIntegerUtil::Len(bStart, bEnd));

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
  };
}

#endif
