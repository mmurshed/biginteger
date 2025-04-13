/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef COMMON_ALGORITHMS
#define COMMON_ALGORITHMS

#include <vector>
#include <string>
using namespace std;

#include "../BigIntegerUtil.h"
#include "addition/ClassicAddition.h"
#include "multiplication/ClassicMultiplication.h"
#include "divison/KnuthDivision.h"

namespace BigMath
{
  // All operations are unsigned
  class CommonAlgorithms
  {

  public:
    static vector<DataT> ConvertBase(
        vector<DataT> const &bigIntB1,
        BaseT base1,
        BaseT base2)
    {
      return ConvertBase(
          bigIntB1, 0, bigIntB1.size() - 1,
          base1,
          base2);
    }

    // Convert an integer from base1 to base2
    // Donald E. Knuth 4.4
    static vector<DataT> ConvertBase(
        vector<DataT> const &bigIntB1, SizeT start, SizeT end,
        BaseT base1,
        BaseT base2)
    {
      if (base1 == base2)
      {
        vector<DataT> bigIntB2(bigIntB1);
        return bigIntB2;
      }

      vector<DataT> bigIntB2(1, 0);

      for (Int i = (Int)end; i >= (Int)start; i--)
      {
        ClassicMultiplication::MultiplyTo(
            bigIntB2,
            base1,
            base2);
        ClassicAddition::AddTo(
            bigIntB2,
            bigIntB1[i],
            base2);
      }

      return bigIntB2;
    }

    // Returns an integer by shifting n places
    // Equivalent to a * B^n

    /*
  •	bigInt[0] = coefficient of B^0  (least significant digit),
  •	bigInt[1] = coefficient of B^1,
  •	…
  •	bigInt[n - 1] = coefficient of B^{n-1} (most significant digit).

Hence, “shifting left by shift digits” means multiplying by B^shift.
In little‐endian form, you want to add shift zeros at the start (indices 0 to shift - 1)
so that the old digit at bigInt[i] (coefficient of B^i) now appears at index i + shift (coefficient of B^{i + shift}).

An easy way to see this is with an example. Suppose
  •	bigInt = [4, 3, 2, 1] (little‐endian).
  •	That represents 4 * B^0 + 3 * B^1 + 2 * B^2 + 1 * B^3.

Shifting left by 2 means the result should be:
4 * B^2 + 3 * B^3 + 2 * B^4 + 1 * B^5,
which in little‐endian form is
[0,0,4,3,2,1].
(Zero for B^0, zero for B^1, then 4 at B^2, 3 at B^3, etc.)
    */
    static vector<DataT> ShiftLeft(
        vector<DataT> const &bigInt,
        SizeT shift)
    {
      if (shift == 0 || BigIntegerUtil::IsZero(bigInt))
        return bigInt;

      SizeT size = (SizeT)bigInt.size() + shift;
      vector<DataT> result(size, 0);
      // Copy
      for (DataT i = 0; i < bigInt.size(); i++)
      {
        result[i + shift] = bigInt[i];
      }

      return result;
    }

    static vector<DataT> ShiftRight(const vector<DataT> &bigInt, SizeT shift)
    {
      if (shift == 0 || BigIntegerUtil::IsZero(bigInt))
        return bigInt;

      // If shifting by more digits than exist, return 0.
      if (shift >= bigInt.size())
      {
        return vector<DataT>{0};
      }

      // The new size is reduced by the shift amount.
      SizeT newSize = bigInt.size() - shift;
      vector<DataT> result(newSize, 0);

      // Copy the higher-order digits from bigInt.
      // Since bigInt is little-endian, we start copying from index "shift".
      for (SizeT i = shift; i < bigInt.size(); i++)
      {
        result[i - shift] = bigInt[i];
      }

      return result;
    }
  };
}

#endif
