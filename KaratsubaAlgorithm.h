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

    // Algorithm from paper
    // "Storage Allocation for the Karatsuba Integer Multiplication Algorithm"
    // by Roman E. Maeder
    //
    // Return the product of a and b to c.
    // Assumption: la >= lb > [la/2].
    // c must be la + lb + 1 in size
    // The array w is used as work array (temporary storage)
    // Karatsuba() adds the product of a and b to the old contents of c.
    // The output array c is assumed cleared.
    // Recursively this precondition is established in lines 16 and 21.
    // This formulation of Karatsuba's algorithm avoids negative intermediate results. 
    // Such results can occur in the formulation in [3], where the sums 
    // ah - al and bl - bh are formed instead of ah + al and bl + bh.

    // The final result has, of course, at most l(a) + l(b) digits. 
    // The intermediate result after line 18 could cause a carry into the next digit.
    // Therefore we must allocate and clear one digit more. This is recursively done in lines 15 and 20.

    static void MultiplyUnsignedRecursive(
      vector<DataT> const& a, SizeT aStart, SizeT aEnd, 
      vector<DataT> const& b, SizeT bStart, SizeT bEnd, 
      vector<DataT>& c, SizeT cStart, 
      vector<DataT>& w, SizeT wStart, 
      ULong base)
      {
        SizeT la = aEnd - aStart + 1;
        SizeT lb = bEnd - bStart + 1;
        
        if(la <= 5)
        {
          // Use naive method 
          ClassicalAlgorithms::MultiplyUnsigned(
            a, aStart, aEnd, 
            b, bStart, bEnd, 
            c, cStart,
            base);
          return;
        }
        
        SizeT m = (la + 1) / 2;

        // c = al * bl + ( (al + ah)(bl + bh) - al * bl - ah * b_h) * B^m + ah * bh*B^2m
        // c = al * bl + (al + ah)(bl + bh) * B^m - al * bl * B^m - ah * bh * B^m + ah * bh * B^2m
        // c = (al + ah)(bl + bh) * B^m
        // c += ah * bh * B^2m
        // c -= ah * bh * B^m
        // c += al * bl
        // c -= al * bl * B^m

        // Clear carry digit
        w[wStart + m] = 0;
        // Save al + ah into w_0,...,w_m
        ClassicalAlgorithms::AddUnsigned(
          a, aStart, aStart + m - 1, // al
          a, aStart + m, aEnd, // ah
          w, wStart, // wl = al + ah
          base);
  
        // Clear carry digit
        w[wStart + m + m + 1] = 0;
      
        // Save bl + bh into w_m+1, ... w_2m+1
        ClassicalAlgorithms::AddUnsigned(
          b, bStart, bStart + m - 1, // bl
          b, bStart + m, bEnd, // bh
          w, wStart + m + 1, // wh = bl + bh
          base);

        // Compute (al + ah)(bl + bh) into c_m, ..., c_3m+1
        // Compute w_0,...,w_m times w_m+1,...,w_2m+1 into c_m, ..., c_3m+1
        // c = (al + ah)(bl + bh) * B^m
        MultiplyUnsignedRecursive(
          w, wStart, wStart + m, // wl
          w, wStart + m + 1, wStart + m + m + 1, // wh
          c, cStart + m, // c = wl * wh
          w, wStart + 2 * (m + 1), // work array
          base);
  
        // Space needed for ah * bh
        SizeT lt = (la - m) + (lb - m) + 1;
        
        // Clear result array
        BigIntegerUtil::SetBit(w, wStart, wStart + lt);

        // Compute a_h * b_h into w_0, ... ,w_(la+lb-2m-1)
        // w = ah * bh
        MultiplyUnsignedRecursive(
          a, aStart + m, aEnd, // ah
          b, bStart + m, bEnd,  // bh
          w, wStart, // w = ah * bh
          w, wStart + lt, // work array
          base);
  
        SizeT wEnd = wStart + (la - m) + (lb - m) - 1;
        // Add ah * bh * B^2m
        // c += ah * bh * B^2m
        ClassicalAlgorithms::AddToUnsigned(
          c, cStart + m + m, cStart + la + lb - 1, // c
          w, wStart, wEnd, // ah * bh
          base);
  
        // Subtract ah * bh * B^m
        // c -= ah * bh * B^m
        ClassicalAlgorithms::SubtractFromUnsigned(
          c, cStart + m, c.size() - 1, // c
          w, wStart, wEnd, // w
          base);
        
        // Space needed for al * bl
        lt = m + m + 1;
        
        // Clear result array
        BigIntegerUtil::SetBit(w, wStart, wStart + lt - 1);
 
        // Compute al * bl into w_0, ..., w_2m-1
        // w = al * bl
        MultiplyUnsignedRecursive(
          a, aStart, aStart + m - 1, // al 
          b, bStart, bStart + m - 1, // bl
          w, wStart, // w = al * bl
          w, wStart + lt, // work array
          base);
        
        wEnd = wStart + m + m - 1;

        // Add al * bl        
        // c += al * bl
        ClassicalAlgorithms::AddToUnsigned(
          c, cStart, cStart + m + m - 1, // c
          w, wStart, wEnd, // w
          base);
        
        // Subtract al * bl * B^m
        // c -= al * bl * B^m
        ClassicalAlgorithms::SubtractFromUnsigned(
          c, cStart + m, c.size() - 1, // c
          w, wStart, wEnd, // w
          base);
      }

    static vector<DataT>& MultiplyUnsigned(vector<DataT> const& a, vector<DataT> const& b, ULong base)
    {
      SizeT size = max(a.size(), b.size());

      size = 2 * size + 1;

      vector<DataT>& c = *new vector<DataT>(size);
      vector<DataT> w(size);

      BigIntegerUtil::SetBit(c, 0, size - 1);
      BigIntegerUtil::SetBit(w, 0, size - 1);

      MultiplyUnsignedRecursive(a, 0, a.size() - 1, b, 0, b.size() - 1, c, 0, w, 0, base);
      
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