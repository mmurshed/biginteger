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

    // Copies the digits a_0, ... ,a_n-1 into r_0, ..., r_n-1.
    static void CopyTo(vector<DataT>& r, SizeT rStart, vector<DataT> const& a, SizeT aStart, SizeT n)
    {
      for(SizeT i = 0; i < n; i++)
      {
        r[rStart + i] = a[aStart + i];
      }
    }

    // Sets l Elements of the array r to zero.
    static void Clear(vector<DataT>& r, SizeT rStart, SizeT n)
    {
      for(SizeT i = 0; i < n; i++)
      {
        r[rStart + i] = 0;
      }
    }


    // Return the product of a and b to c.
    // Assumption: la >= lb > [la/2].
    // c must be la + lb + 1 in size
    // The array w is used as work array (temporary storage)
    // Karatsuba() adds the product of a and b to the old contents of c.
    // The output array c is assumed cleared.
    // Recursively this precondition is established in lines 16 and 21.
    // This formulation of Karatsuba's algorithm avoids negative intermediate results. 
    // Such results can occur in the formulation in [3], where the sums 
    // a_h - a_l and b_l - b_h are formed instead of a_h + a_l and b_l + b_h.

    // The final result has, of course, at most l(a) + l(b) digits. 
    // The intermediate result after line 18 could cause a carry into the next digit.
    // Therefore we must allocate and clear one digit more. This is recursively done in lines 15 and 20.

    static void Karatsuba(
      vector<DataT> const& a, SizeT aStart, SizeT aEnd, 
      vector<DataT> const& b, SizeT bStart, SizeT bEnd, 
      vector<DataT>& c, SizeT cStart, 
      vector<DataT>& w, SizeT wStart, 
      ULong base)
      {
        SizeT la = aEnd - aStart + 1;
        SizeT lb = bEnd - bStart + 1;
        
        if(la <= 4)
        {
          // Use naive method
          ClassicalAlgorithms::MultiplyUnsigned(a, aStart, aEnd, b, bStart, bEnd, c, cStart, base);
          return;
        }
        
        SizeT m = (la + 1) / 2;

        // Save a_l + a_h into w_0,...,w_m
        ClassicalAlgorithms::AddUnsigned(a, aStart, aStart + m - 1, a, aStart + m, aEnd, w, wStart, base);
  
        // Save b_l + b_h into w_m+1, ... w_2m+1
        ClassicalAlgorithms::AddUnsigned(b, bStart, bStart + m - 1, b, bStart + m, bEnd, w, wStart + m + 1, base);

        // Compute (a_l + a_h)(b_l + b_h) into c_m, ..., c_3m+1
        // Compute w_0,...,w_m times w_m+1,...,w_2m+1 into c_m, ..., c_3m+1
        Karatsuba(w, wStart, wStart + m, w, wStart + m + 1, wStart + 2 * m + 1, c, cStart + m, w, wStart + 2 * (m + 1), base);
  
        // Space needed for a_h * b_h
        SizeT lt = (la - m) + (lb - m) + 1;
        
        // Clear result array
        Clear(w, wStart, wStart + lt);

        // Compute a_h * b_h into w_0, ... ,w_(la+lb-2m-1)
        Karatsuba(a, aStart + m, aEnd, b, bStart + m, bEnd, w, wStart, w, wStart + lt, base);
  
        // Add a_h * b_h * B^2m
        ClassicalAlgorithms::AddToUnsigned(c, cStart + 2 * m, cStart + la + lb + 2 * m, w, wStart, wStart + la + lb + 2 * m, base);
  
        // Subtract a_h * b_h * B^m
        ClassicalAlgorithms::SubtractFromUnsigned(c, cStart + m, cStart + la + lb + m, w, wStart, wStart + la + lb + m, base);
        
        // Space needed for a_l * b_l
        lt = 2 * m + 1;
        
        // Clear result array
        Clear(w, wStart, wStart + lt);
 
        // Compute a_l * b_l into w_0, ..., w_2m-1
        Karatsuba(a, aStart, aStart + m, b, bStart, bStart + m, w, wStart, w, wStart + lt, base);
        
        // Add a_l * b_l
        ClassicalAlgorithms::AddToUnsigned(c, cStart, cStart + m, w, wStart, wStart + 2 * m - 1, base);
        
        // Subtract a_l * b_l * B^m
        ClassicalAlgorithms::SubtractFromUnsigned(c, cStart, cStart + m, w, wStart, wStart + 2 * m - 1, base);
      }



// The procedure addto(digit *r, digit *a, int l), adds the digits a_0, ..., a_l-1 to the old contents of r_0, ..., r_l-1, 
// propagating the carry as usual. Any carry out of the last place rl-1 is propagated further into rl, rl+l, .... 
// This is the reason for setting the next higher digit to zero in lines 9 and 12.

// The procedure subfrom(digit *r, digit *a, int l) subtracts a from the old contents of r, 
// propagating the borrow as usual. Any borrow out of the last place r_l-1 is propagated further. 
// The way this is called in line 19 there will never be a negative result since line 18 guarantees
// that at least one higher element of c is nonzero and can therefore absorb the borrow.


    // static vector<DataT> *x1;
    // static vector<DataT> *x2;
    // static vector<DataT> *y1;
    // static vector<DataT> *y2;
    // static vector<DataT> *a;
    // static vector<DataT> *b;
    // static vector<DataT> *c;
    // static vector<DataT> *p;
    // static vector<DataT> *q;
    // static vector<DataT> *r;
    // static vector<DataT> *pq;
    // static vector<DataT> *xy;

    // static void Init(SizeT size)
    // {
    //   ++size;
    //   x1 = new vector<DataT>(size);
    //   x2 = new vector<DataT>(size);
    //   y1 = new vector<DataT>(size);
    //   y2 = new vector<DataT>(size);

    //   a = new vector<DataT>(size);
    //   b = new vector<DataT>(size);
    //   c = new vector<DataT>(size);

    //   p = new vector<DataT>(size);
    //   q = new vector<DataT>(size);
    //   r = new vector<DataT>(size);

    //   size *= 2;
    //   pq = new vector<DataT>(size);
    //   xy = new vector<DataT>(size);
    // }

    // static void MultiplyUnsignedR(
    //   vector<DataT> const& x, SizeT xStart, SizeT xEnd, 
    //   vector<DataT> const& y, SizeT yStart, SizeT yEnd, 
    //   vector<DataT>& result, SizeT rStart, ULong base)
    // {
    //   if(BigIntegerUtil::IsZero(x, xStart, xEnd) || BigIntegerUtil::IsZero(y, yStart, yEnd))
    //     return *new vector<DataT>();
      
    //   SizeT size = max(x.size(), y.size());

    //   // bottom of the recursion
    //   if (size <= 4)
    //   {
    //     ClassicalAlgorithms::MultiplyUnsigned(x, y, base);
    //     return;
    //   }

    //   while(x.size() < size)
    //     x.push_back(0);

    //   while(y.size() < size)
    //     y.push_back(0);

    //   // x = x1 * B^m + x2 
    //   // y = y1 * B^m + y2
    //   // xy = (x1 * B^m + x2)(y1 * B^m + y2)
    //   //    = (x1 * y1) * B^2m + (x1 * y2 + x2 * y1) * B^m + x2 * y2
    //   //    =     a     * B^2m +          b          * B^m +    c
    //   // Where,
    //   // a = x1 * y1
    //   // b = x1 * y2 + x2 * y1
    //   // c = x2 * y2

    //   // However, b can be calcluated from
    //   // b = (x1 + x2)(y1 + y2) - a - c
      
    //   SizeT half = size / 2;

    //   // Split x in half for x1 and x2
    //   // Vector contains the number is reverse order
    //   // So first half is x2 and second half is x1
    //   vector<DataT>& x2 = *new vector<DataT>(x.begin(), x.begin() + half);
    //   vector<DataT>& x1 = *new vector<DataT>(x.begin() + half, x.end());

    //   // Split y in half for y1 and y2
    //   // Vector contains the number is reverse order
    //   // So first half is y2 and second half is y1
    //   vector<DataT>& y2 = *new vector<DataT>(y.begin(), y.begin() + half);
    //   vector<DataT>& y1 = *new vector<DataT>(y.begin() + half, y.end());

    //   // a = x1 * y1
    //   vector<DataT>& a = MultiplyUnsignedR(x1, y1, base);
      
    //   // c = x2 * y2
    //   vector<DataT>& c = MultiplyUnsignedR(x2, y2, base);

    //   // b = (x1 + x2)(y1 + y2) - a - c
    //   //   = (x1 + x2)(y1 + y2) - (a + c)
    //   //   =     p   *    q   -      r

    //   // p = x1 + x2
    //   vector<DataT>& p = ClassicalAlgorithms::AddUnsigned(x1, x2, base);
    //   // q = y1 + y2
    //   vector<DataT>& q = ClassicalAlgorithms::AddUnsigned(y1, y2, base);
    //   // r = a + c
    //   vector<DataT>& r = ClassicalAlgorithms::AddUnsigned(a, c, base);
      
    //   // pq = p * q = (x1 + x2)(y1 + y2)
    //   vector<DataT>& pq = MultiplyUnsignedR(p, q, base);

    //   // b = pq - r
    //   vector<DataT>& b = *new vector<DataT>(0);
    //   bool isNegative = false;
    //   Int cmp = ClassicalAlgorithms::UnsignedCompareTo(pq, r);
    //   if(cmp > 0)
    //     b =  ClassicalAlgorithms::SubtractUnsigned(pq, r, base);
    //   else
    //   {
    //     isNegative = true;
    //     b =  ClassicalAlgorithms::SubtractUnsigned(r, pq, base);
    //   }


    //   // a_hat = a * B^2m
    //   // A faster way to do this is to insert zeros at the end of the number
    //   // i.e. Shift left.
    //   // 856 * 10^3 = 856000
    //   // This is true for any base
    //   vector<DataT>& a_hat = ClassicalAlgorithms::ShiftLeft(a, 2 * half);

    //   // b_hat = b * B^m
    //   vector<DataT>& b_hat = ClassicalAlgorithms::ShiftLeft(b, half);

    //   // g = a * B^2m +  c
    //   //   =  a_hat   +  c
    //   vector<DataT>& g = ClassicalAlgorithms::AddUnsigned(a_hat, c, base);
      
    //   // xy = a * B^2m +/- b * B^m + c
    //   //    =    g     +/- b_hat
    //   vector<DataT>& xy = *new vector<DataT>(0);
    //   if(isNegative)
    //     xy = ClassicalAlgorithms::SubtractUnsigned(g, b_hat, base);
    //   else
    //     xy = ClassicalAlgorithms::AddUnsigned(g, b_hat, base);

    //   BigIntegerUtil::TrimZeros(xy);

    //   return xy;
    // }

    static vector<DataT>& MultiplyUnsigned(vector<DataT> const& a, vector<DataT> const& b, ULong base)
    {
      SizeT size = max(a.size(), b.size());

      size = 2 * size + 1;

      vector<DataT> c(size);
      vector<DataT> w(size);

      Clear(c, 0, size);
      Clear(w, 0, size);

      Karatsuba(a, 0, a.size(), b, 0, b.size(), c, 0, w, 0, base);

      return c;
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