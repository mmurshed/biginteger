/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef KARATSUBA_ALGORITHM_H
#define KARATSUBA_ALGORITHM_H

#include <vector>
using namespace std;

#include "BigInteger.h"
#include "BigIntegerUtil.h"
#include "ClassicalAlgorithms.h"

namespace BigMath
{
  class KaratsubaAlgorithm
  {
    private:
    static vector<DataT>& MultiplyUnsignedR(vector<DataT>& x, vector<DataT>& y, ULong base)
    {
      SizeT size = max(x.size(), y.size());

      // bottom of the recursion
      if (size <= 4)
      {
        return ClassicalAlgorithms::MultiplyUnsigned(x, y, base);
      }

      while(x.size() < size)
        x.push_back(0);

      while(y.size() < size)
        y.push_back(0);

      // x = x1 * B^m + x2 
      // y = y1 * B^m + y2
      // xy = (x1 * B^m + x2)(y1 * B^m + y2)
      //    = (x1 * y1) * B^2m + (x1 * y2 + x2 * y1) * B^m + x2 * y2
      //    =     a     * B^2m +          b          * B^m +    c
      // Where,
      // a = x1 * y1
      // b = x1 * y2 + x2 * y1
      // c = x2 * y2

      // However, b can be calcluated from
      // b = (x1 + x2)(y1 + y2) - a - c
      
      SizeT half = size / 2;

      // Split x in half for x1 and x2
      // Vector contains the number is reverse order
      // So first half is x2 and second half is x1
      vector<DataT>& x2 = *new vector<DataT>(x.begin(), x.begin() + half + 1);
      vector<DataT>& x1 = *new vector<DataT>(x.begin() + half + 1, x.end());

      // Split y in half for y1 and y2
      // Vector contains the number is reverse order
      // So first half is y2 and second half is y1
      vector<DataT>& y2 = *new vector<DataT>(y.begin(), y.begin() + half + 1);
      vector<DataT>& y1 = *new vector<DataT>(y.begin() + half + 1, y.end());

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
      vector<DataT>& b = *new vector<DataT>(0);
      bool isNegative = false;
      Int cmp = ClassicalAlgorithms::UnsignedCompareTo(pq, r);
      if(cmp > 0)
        b =  ClassicalAlgorithms::SubtractUnsigned(pq, r, base);
      else
      {
        isNegative = true;
        b =  ClassicalAlgorithms::SubtractUnsigned(r, pq, base);
      }


      // a_hat = a * B^2m
      // A faster way to do this is to insert zeros at the end of the number
      // i.e. Shift left.
      // 856 * 10^3 = 856000
      // This is true for any base
      vector<DataT>& a_hat = ClassicalAlgorithms::ShiftLeft(a, 2 * half);

      // b_hat = b * B^m
      vector<DataT>& b_hat = ClassicalAlgorithms::ShiftLeft(b, half);

      // g = a * B^2m +  c
      //   =  a_hat   +  c
      vector<DataT>& g = ClassicalAlgorithms::AddUnsigned(a_hat, c, base);
      
      // xy = a * B^2m +/- b * B^m + c
      //    =    g     +/- b_hat
      vector<DataT>& xy = *new vector<DataT>(0);
      if(isNegative)
        xy = ClassicalAlgorithms::SubtractUnsigned(g, b_hat, base);
      else
        xy = ClassicalAlgorithms::AddUnsigned(g, b_hat, base);

      BigIntegerUtil::TrimZeros(xy);

      return xy;
    }

    static vector<DataT>& MultiplyUnsigned(vector<DataT> const& a, vector<DataT> const& b, ULong base)
    {
      SizeT size = max(a.size(), b.size());
      vector<DataT> x(a);
      vector<DataT> y(b);

      vector<DataT> & result = MultiplyUnsignedR(x, y, base);

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
      
      return result;
    } 
   };
}

#endif