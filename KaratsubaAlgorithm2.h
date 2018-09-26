/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef KARATSUBA_ALGORITHM2_H
#define KARATSUBA_ALGORITHM2_H

#include <vector>
using namespace std;

#include "BigInteger.h"
#include "BigIntegerUtil.h"
#include "ClassicalAlgorithms.h"

namespace BigMath
{
  class KaratsubaAlgorithm2
  {
    private:
    static vector<DataT>& MultiplyUnsignedR(vector<DataT> const& x, vector<DataT> const& y, ULong base)
    {    
      // bottom of the recursion
      if (x.size() <= 4 && y.size() <= 4)
      {
        return ClassicalAlgorithms::MultiplyUnsigned(x, y, base);
      }

      // x = x1 * B^m + x2 
      // y = y1 * B^m + y2
      // xy = (x1 * B^m + x2)(y1 * B^m + y2)
      //    = (x1 * y1) * B^2m + (x1 * y2 + x2 * y1) * B^m + x2 * y2
      //    =     a     * B^2m +          b          * B^m +    c
      // Where,
      // a = x1 * y1
      // b = x1 * y2 + x2 * y1
      // c = x2 * y2

      // However, b canbe calcluated from
      // b = (x1 + x2)(y1 + y2) - a - c
      
      SizeT half = x.size() / 2;

      // Split x in half for x1 and x2
      vector<DataT>& x1 = *new vector<DataT>(x.begin(), x.begin() + half);
      vector<DataT>& x2 = *new vector<DataT>(x.begin() + half, x.end());

      half = y.size() / 2;

      // Split y in half for y1 and y2
      vector<DataT>& y1 = *new vector<DataT>(y.begin(), y.begin() + half);
      vector<DataT>& y2 = *new vector<DataT>(y.begin() + half, y.end());

      // a = x1 * y1
      vector<DataT>& a = MultiplyUnsignedR(x1, y1, base);
      
      // c = x2 * y2
      vector<DataT>& c = MultiplyUnsignedR(x2, y2, base);

      // b = (x1 + x2)(y1 + y2) - a - c
      //   = (x1 + x2)(y1 + y2) - (a + c)
      //   =     p   *    q   -      r

      // p = x1 + x2
      vector<DataT>& p = ClassicalAlgorithms::AddUnsigned(x1, x2, base);
      // q = y1 + y2
      vector<DataT>& q = ClassicalAlgorithms::AddUnsigned(y1, y2, base);
      // r = a + c
      vector<DataT>& r = ClassicalAlgorithms::AddUnsigned(a, c, base);
      
      // pq = p * q = (x1 + x2)(y1 + y2)
      vector<DataT>& pq = MultiplyUnsignedR(p, q, base);

      // b = pq - r
      vector<DataT>& b =  ClassicalAlgorithms::SubtractUnsigned(pq, r, base);

      // a_hat = a * B^2m
      // A faster way to do this is to insert zeros at the end of the number
      // 856 * 10^3 = 856000
      // This is true for any base
      vector<DataT>& a_hat = ClassicalAlgorithms::ShiftLeft(a, 2 * half);

      // b_hat = b * B^m
      vector<DataT>& b_hat = ClassicalAlgorithms::ShiftLeft(b, half);

      // for(SizeT i = 0; i < half; i++) {
      //   ClassicalAlgorithms::MultiplyToUnsigned(b, base, base);
      // }

      // g = a * B^2m + b * B^m
      //   =  a_hat   +  b_hat
      vector<DataT>& g = ClassicalAlgorithms::AddUnsigned(a_hat, b_hat, base);
      
      // xy = a * B^2m + b * B^m + c
      vector<DataT>& xy = ClassicalAlgorithms::AddUnsigned(g, c, base);

      return xy;
    }

    static vector<DataT>& MultiplyUnsigned(vector<DataT> const& a, vector<DataT> const& b, ULong base)
    {
      vector<DataT> aa(a);
      vector<DataT> bb(b);
      SizeT size = max(a.size(), b.size());

      while(aa.size() < size)
        aa.push_back(0);

      while(bb.size() < size)
        bb.push_back(0);

      vector<DataT> & result = MultiplyUnsignedR(aa, bb, base);

      return result;
    }

    static BigInteger& MultiplyUnsigned(BigInteger const& a, BigInteger const& b)
    {
      vector<DataT>& result = MultiplyUnsigned(a.GetInteger(), b.GetInteger(), BigInteger::Base());
      return *new BigInteger(result, false);
    }

public: 
    static BigInteger& Multiply(BigInteger const& a, BigInteger const& b)
    {
      if(a.IsZero() || b.IsZero())
        return *new BigInteger(); // 0 times anything is zero

      BigInteger& result = MultiplyUnsigned(a, b);
      if(a.IsNegative() != b.IsNegative())
        result.SetSign(true);

      result.TrimZeros();
      
      return result;
    } 
   };
}

#endif