/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef KARATSUBA_MULTIPLICATION
#define KARATSUBA_MULTIPLICATION

#include <vector>
using namespace std;

#include "../../BigInteger.h"
#include "../../BigIntegerUtil.h"
#include "../addition/ClassicAddition.h"
#include "../subtraction/ClassicSubtraction.h"
#include "../multiplication/ClassicMultiplication.h"

namespace BigMath
{
  class KaratsubaMultiplication
  {
  private:
    // Cut off deterioriate over 32
    static const SizeT KARATSUBA_THRESHOLD = 32;

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

    static void MultiplyRecursive(
        vector<DataT> const &a, SizeT aStart, SizeT aEnd,
        vector<DataT> const &b, SizeT bStart, SizeT bEnd,
        vector<DataT> &c, SizeT cStart,
        vector<DataT> &w, SizeT wStart,
        BaseT base)
    {
      Int la = BigIntegerUtil::Len(aStart, aEnd);
      Int lb = BigIntegerUtil::Len(bStart, bEnd);

      if (la <= KARATSUBA_THRESHOLD)
      {
        // Use naive method
        ClassicMultiplication::Multiply(
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
      w.at(wStart + m) = 0;

      // Save al + ah into w_0,...,w_m
      ClassicAddition::Add(
          a, aStart, aStart + m - 1, // al
          a, aStart + m, aEnd,       // ah
          w, wStart,                 // wl = al + ah
          base);

      // Clear carry digit
      w.at(wStart + m + m + 1) = 0;

      // Save bl + bh into w_m+1, ... w_2m+1
      ClassicAddition::Add(
          b, bStart, bStart + m - 1, // bl
          b, bStart + m, bEnd,       // bh
          w, wStart + m + 1,         // wh = bl + bh
          base);

      // Compute (al + ah)(bl + bh) into c_m, ..., c_3m+1
      // Compute w_0,...,w_m times w_m+1,...,w_2m+1 into c_m, ..., c_3m+1
      // c = (al + ah)(bl + bh) * B^m
      MultiplyRecursive(
          w, wStart, wStart + m,                 // wl
          w, wStart + m + 1, wStart + m + m + 1, // wh
          c, cStart + m,                         // c = wl * wh
          w, wStart + 2 * (m + 1),               // work array
          base);

      // Space needed for ah * bh
      SizeT lt = (la - m) + (lb - m) + 1;

      // Clear result array
      BigIntegerUtil::SetBit(w, wStart, wStart + lt - 1);

      // Compute a_h * b_h into w_0, ... ,w_(la+lb-2m-1)
      // w = ah * bh
      MultiplyRecursive(
          a, aStart + m, aEnd, // ah
          b, bStart + m, bEnd, // bh
          w, wStart,           // w = ah * bh
          w, wStart + lt,      // work array
          base);

      SizeT wEnd = wStart + (la - m) + (lb - m) - 1;
      // Add ah * bh * B^2m
      // c += ah * bh * B^2m
      ClassicAddition::AddTo(
          c, cStart + m + m, cStart + la + lb - 1, // c
          w, wStart, wEnd,                         // ah * bh
          base);

      // Subtract ah * bh * B^m
      // c -= ah * bh * B^m
      ClassicSubtraction::SubtractFrom(
          c, cStart + m, (SizeT)c.size() - 1, // c
          w, wStart, wEnd,                    // w
          base);

      // Space needed for al * bl
      lt = m + m + 1;

      // Clear result array
      BigIntegerUtil::SetBit(w, wStart, wStart + lt - 1);

      // Compute al * bl into w_0, ..., w_2m-1
      // w = al * bl
      MultiplyRecursive(
          a, aStart, aStart + m - 1, // al
          b, bStart, bStart + m - 1, // bl
          w, wStart,                 // w = al * bl
          w, wStart + lt,            // work array
          base);

      wEnd = wStart + m + m - 1;

      // Add al * bl
      // c += al * bl
      ClassicAddition::AddTo(
          c, cStart, cStart + m + m - 1, // c
          w, wStart, wEnd,               // w
          base);

      // Subtract al * bl * B^m
      // c -= al * bl * B^m
      ClassicSubtraction::SubtractFrom(
          c, cStart + m, (SizeT)c.size() - 1, // c
          w, wStart, wEnd,                    // w
          base);
    }

  public:
    static vector<DataT> Multiply(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base)
    {
      SizeT n = (SizeT)max(a.size(), b.size());
      SizeT size = 3 * n;

      vector<DataT> c(size, 0);
      vector<DataT> w(size, 0);

      vector<DataT> x(n, 0);
      vector<DataT> y(n, 0);

      BigIntegerUtil::Copy(a, x);
      BigIntegerUtil::Copy(b, y);

      MultiplyRecursive(
          x, 0, (SizeT)x.size() - 1, // b
          y, 0, (SizeT)y.size() - 1, // a
          c, 0,                      // c = b * a
          w, 0,                      // work array
          base);

      BigIntegerUtil::TrimZeros(c);

      return c;
    }
  };
}

#endif
