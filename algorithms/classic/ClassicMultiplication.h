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

#include "../BigIntegerUtil.h"

namespace BigMath
{
  // All operations are unsigned
  class ClassicMultiplication 
  {
    public:
    static void MultiplyTo(
      vector<DataT>& a,
      DataT b,
      ULong base)
    {
        MultiplyTo(
          a, 0, a.size() - 1,
          b,
          base);
    }    

    static SizeT MultiplyTo(
      vector<DataT>& a, SizeT aStart, SizeT aEnd,
      DataT b,
      ULong base)
    {
      if(b == 0 || // a times 0
        BigIntegerUtil::IsZero(a, aStart, aEnd)) // 0 times b
      {
        BigIntegerUtil::SetBit(a, aStart, aEnd, 0);
        return 0;
      }

      ULong carry = 0;

      for(Int j = aStart; j <= aEnd; j++)
      {
        ULong multiply = a.at(j);
        multiply *= b;
        multiply += carry;
        
        if(j < a.size())
        {
          a.at(j) = (DataT)(multiply % base);
        }
        else
        {
          a.push_back( (DataT)(multiply % base) );
        }
          
        carry = multiply / base;
      }

      Int j = aEnd + 1;
      while(carry > 0)
      {
        if(j < a.size())
        {
          a.at(j) = (DataT)(carry % base);
        }
        else
        {
          a.push_back( (DataT)(carry % base) );
        }
        
        carry = carry / base;
        j++;
      }

      BigIntegerUtil::TrimZeros(a);

      return j;
    }    

    static vector<DataT> Multiply(
      vector<DataT> const& a,
      DataT b,
      ULong base)
      {
        vector<DataT> w(a.size());
        Multiply(
          a, 0, a.size() - 1,
          b,
          w, 0, w.size() - 1,
          base);
        BigIntegerUtil::TrimZeros(w);
        return w;
      }

    static SizeT  Multiply(
      vector<DataT> const& a, SizeT aStart, SizeT aEnd,
      DataT b,
      vector<DataT>& w, SizeT wStart, SizeT wEnd,
      ULong base)
    {
      if(b == 0 || // a times 0
        BigIntegerUtil::IsZero(a, aStart, aEnd)) // 0 times b
      {
        BigIntegerUtil::SetBit(w, wStart, wEnd, 0);
        return 0;
      }

      SizeT len = BigIntegerUtil::Len(aStart, aEnd);
      SizeT wPos = wStart;

      ULong carry = 0;
      for(Int j = 0 ; j < len; j++)
      {
        ULong multiply = 0;
        SizeT aPos = aStart + j;
        if(aPos <= aEnd)
          multiply = a.at(aPos);
        multiply *= b;
        multiply += carry;
        
        wPos = wStart + j;
        if(wPos < w.size())
        {
          w.at(wPos) = (DataT)(multiply % base);
        }
        else
        {
          w.push_back( (DataT)(multiply % base) );
        }
          
        carry = multiply / base;
      }

      wPos = wStart + len;
      while(carry > 0)
      {
        if(wPos < w.size())
        {
          w.at(wPos) = (DataT)(carry % base);
        }
        else
        {
          w.push_back( (DataT)(carry % base) );
        }
        
        carry = carry / base;
        wPos++;
      }

      return wPos;
    }

    public:
    static void MultiplyTo(
      vector<DataT> & a,
      vector<DataT> const& b,
      ULong base)
    {
      a = Multiply(a, b, base);
    }

    static vector<DataT> Multiply(
      vector<DataT> const& a,
      vector<DataT> const& b,
      ULong base)
    {
      if(BigIntegerUtil::IsZero(a) || BigIntegerUtil::IsZero(b)) // 0 times
        return vector<DataT>();

      SizeT size = (SizeT)(a.size() + b.size() + 1);
      vector<DataT> result(size);

      Multiply(a, 0, (SizeT)a.size() - 1, b, 0, (SizeT)b.size() - 1, result, 0, base);

      BigIntegerUtil::TrimZeros(result);
      return result;      
    }

    // Implentation of multiplication by paper-pencil method
    // Classical algorithm
    // Runtime O(n^2), Space O(n)
    static SizeT Multiply(
      vector<DataT> const& a, SizeT aStart, SizeT aEnd, 
      vector<DataT> const& b, SizeT bStart, SizeT bEnd, 
      vector<DataT>& result, SizeT rStart,
      ULong base)
    {

      aEnd = min(aEnd, (SizeT) (a.size() - 1));
      bEnd = min(bEnd, (SizeT) (b.size() - 1));

      SizeT k = rStart;
      SizeT lenA = aEnd - aStart + 1;
      
      for(SizeT i = bStart; i <= bEnd; i++)
      {
        ULong carry = 0;
        SizeT jStart = rStart + (i - bStart);
        for(SizeT j = aStart; j <= aEnd; j++)
        {
          SizeT k = jStart + (j - aStart);
          ULong multiply = a.at(j);
          multiply *= b.at(i);
          multiply += result.at(k);
          multiply += carry;
          
          result.at(k) = (DataT)(multiply % base);
          carry = multiply / base;
        }
        k = jStart + lenA;
        result.at(k) = (DataT)carry;
      }
      return k;
    }
   };
}

#endif
