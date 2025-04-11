/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef CLASSIC_ADDITION
#define CLASSIC_ADDITION

#include <vector>
#include <string>
using namespace std;

#include "../BigIntegerUtil.h"

namespace BigMath
{
  // All operations are unsigned
  class ClassicAddition 
  {
    public:
    static void AddTo(
      vector<DataT> & a,
      ULong b,
      ULong base)
      {
        AddTo(a, 0, a.size() - 1, b, base);
      }

    static void AddTo(
      vector<DataT> & a, SizeT aStart, SizeT aEnd, 
      ULong b,
      ULong base)
    {
      ULong sum = b;
      ULong carry = 0;
      if(aEnd >= aStart)
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
      while(carry > 0)
      {
        sum = carry;
        if(aEnd >= aPos)
        {
          sum += a[aPos];
          a[aPos] = (DataT)(sum % base);
        }
        else
          a.push_back((DataT)(sum % base));
        carry = sum / base;
        aPos++;
      }

      while(carry > 0)
      {
        if(aEnd >= aPos)
          a[aPos] = (DataT)(carry % base);
        else
          a.push_back((DataT)(carry % base));
        carry = carry / base;
        aPos++;
      }

      BigIntegerUtil::TrimZeros(a);
    }    

    public:
    static void AddTo(
      vector<DataT>& a, SizeT aStart, SizeT aEnd, 
      vector<DataT> const& b, SizeT bStart, SizeT bEnd,
      ULong base)
    {
      Add(
          a, aStart, aEnd,
          b, bStart, bEnd,
          a, aStart,
          base);
    }

    static vector<DataT> Add(
      vector<DataT> const& a,
      vector<DataT> const& b,
      ULong base)
    {
      SizeT size = (SizeT)max(a.size(),  b.size()) + 1;
      vector<DataT> result(size);

      result[size - 1] = 0;
      
      Add(
          a, 0, (SizeT)a.size() - 1,
          b, 0, (SizeT)b.size() - 1,
          result, 0,
          base);
      
      return result;
    }

    // Implentation of addition by paper-pencil method
    // Runtime O(n), Space O(n)
    static void Add(
      vector<DataT> const& a, SizeT aStart, SizeT aEnd, 
      vector<DataT> const& b, SizeT bStart, SizeT bEnd, 
      vector<DataT>& result, SizeT rStart,
      ULong base)
    {
      // bool aZero = BigIntegerUtil::IsZero(a, aStart, aEnd);
      // bool bZero = BigIntegerUtil::IsZero(b, bStart, bEnd);
      // if(aZero && bZero)
      //   return;
      // else if(aZero)
      // {
      //   BigIntegerUtil::Copy(b, bStart, result, rStart);
      //   return;
      // }

      aEnd = min(
        aEnd,
        (SizeT) (a.size() - 1));

      bEnd = min(
        bEnd,
        (SizeT) (b.size() - 1));

      Int size = max(
        BigIntegerUtil::Len(aStart, aEnd),
        BigIntegerUtil::Len(bStart, bEnd));

      Long carry = 0;

      for(Int i = 0; i < size; i++)
      {
        Long digitOps = 0;

        Int aPos = i + aStart;
        if(aPos <= aEnd && aPos < a.size())
          digitOps = a.at(aPos);

        digitOps += carry;
        
        Int bPos = i + bStart;
        if(bPos <= bEnd && bPos < b.size())
          digitOps += b.at(bPos);

        Int rPos = rStart + i;
        if(rPos < result.size()) 
          result.at(rPos) = (DataT)(digitOps % base);
        carry = digitOps / base;
      }
      Int rPos = rStart + size;
      if(carry > 0 && rPos < result.size())
        result.at(rPos) += carry;

      BigIntegerUtil::TrimZeros(result);
    }
   };
}

#endif
