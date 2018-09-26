/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef KARATSUBA_ALGORITHM_H
#define KARATSUBA_ALGORITHM_H

#include "BigInteger.h"
#include "BigIntegerUtil.h"
#include "ClassicalAlgorithms.h"

namespace BigMath
{
  class KaratsubaAlgorithm
  {
    private:
    static const int KARAT_CUTOFF = 4;

    static void gradeSchool(Long *a, Long *b, Long *ret, SizeT d)
    {
        for(SizeT i = 0; i < 2 * d; i++)
          ret[i] = 0;
        for(SizeT i = 0; i < d; i++)
        {
            for(SizeT j = 0; j < d; j++)
              ret[i + j] += a[i] * b[j];
        }
    }

    // Karatsuba Algorithm
    // http://www.cburch.com/proj/karat/karat.txt

    // ret must have space for 6d digits.
    // the result will be in only the first 2d digits
    // my use of the space in ret is pretty creative.
    // | ar*br  | al*bl  | asum*bsum | lower-recursion space | asum | bsum |
    //  d digits d digits  d digits     3d digits              d/2    d/2
    static void karatsuba(Long *a, Long *b, Long *ret, SizeT d, ULong base)
    {
        Long *ar = &a[0]; // low-order half of a
        Long *al = &a[d/2]; // high-order half of a
        Long *br = &b[0]; // low-order half of b
        Long *bl = &b[d/2]; // high-order half of b
        Long *asum = &ret[d * 5]; // sum of a's halves
        Long *bsum = &ret[d * 5 + d/2]; // sum of b's halves
        Long *x1 = &ret[d * 0]; // ar*br's location
        Long *x2 = &ret[d * 1]; // al*bl's location
        Long *x3 = &ret[d * 2]; // asum*bsum's location

        // when d is small, we're better off just reverting to
        // grade-school multiplication, since it's faster at this point.
        if(d <= KARAT_CUTOFF)
        {
            gradeSchool(a, b, ret, d);
            return;
        }

        // compute asum and bsum
        for(SizeT i = 0; i < d / 2; i++)
        {
            asum[i] = al[i] + ar[i];
            bsum[i] = bl[i] + br[i];
        }

        // do recursive calls (I have to be careful about the order,
        // since the scratch space for the recursion on x1 includes
        // the space used for x2 and x3)
        karatsuba(ar, br, x1, d/2, base);
        karatsuba(al, bl, x2, d/2, base);
        karatsuba(asum, bsum, x3, d/2, base);

        // combine recursive steps
        for(SizeT i = 0; i < d; i++)
          x3[i] = x3[i] - x1[i] - x2[i];
        for(SizeT i = 0; i < d; i++)
          ret[i + d/2] += x3[i];
    }

    static void doCarry(Long *a, SizeT d, ULong base)
    {
        Long c = 0;
        for(SizeT i = 0; i < d; i++)
        {
            a[i] += c;
            if(a[i] < 0)
            {
                c = -(-(a[i] + 1) / base + 1);
            }
            else
            {
                c = a[i] / base;
            }
            a[i] -= c * base;
        }

        // Overflow
        if(c != 0) 
        {} 
    }

    static vector<DataT>& MultiplyUnsigned(vector<DataT> const& a, vector<DataT> const& b, ULong base)
    {
      SizeT size = max(a.size(), b.size());
      SizeT rsize = size * 6;
      
      // Copy a
      Long* aa  = new Long[size];
      for(SizeT i = 0; i < a.size(); i++)
        aa[i] = a[i];
      for(SizeT i = a.size(); i < size; i++)
        aa[i] = 0;

      // Copy b
      Long* bb  = new Long[size];
      for(SizeT i = 0; i < b.size(); i++)
        bb[i] = b[i];
      for(SizeT i = b.size(); i < size; i++)
        bb[i] = 0;

      SizeT d = 1;
      for(d = 1; d < size; d *= 2);

      Long* r  = new Long[rsize];

      karatsuba(aa, bb, r, d, base); // compute product w/o regard to carry
      doCarry(r, 2 * d, base); // now do any carrying

      vector<DataT>& result  = *new vector<DataT>(rsize);
      for(SizeT i = 0; i < rsize; i++)
      {
        if(r[i] == 0)
          break;
        result[i] = r[i];
      }

      // delete [] r;
      // delete [] aa;
      // delete [] bb;

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