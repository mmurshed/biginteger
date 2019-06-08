/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef COMMON_ALGORITHMS
#define COMMON_ALGORITHMS

#include <vector>
#include <string>
using namespace std;

#include "BigIntegerUtil.h"

namespace BigMath
{
  // All operations are unsigned
  class CommonAlgorithms 
  {

    public:
    static vector<DataT> ConvertBase(
      vector<DataT> const& bigIntB1,
      ULong base1,
      ULong base2)
      {
        return ConvertBase(
          bigIntB1, 0, bigIntB1.size() - 1,
          base1,
          base2);
      }

    // Convert an integer from base1 to base2
    // Donald E. Knuth 4.4
    static vector<DataT> ConvertBase(
      vector<DataT> const& bigIntB1, SizeT start, SizeT end,
      ULong base1,
      ULong base2)
    {
      if(base1 == base2)
      {
        vector<DataT> bigIntB2(bigIntB1);
        return bigIntB2;
      }

      vector<DataT> bigIntB2(1, 0);

      for(Int i = (Int)end; i >= (Int)start; i--)
      {
        MultiplyTo(
          bigIntB2,
          base1,
          base2);
        AddTo(
          bigIntB2,
          bigIntB1[i],
          base2);
      }

      return bigIntB2;
    }

    // Returns an integer by shifting n places
    // Equivalent to a * B^n
    static vector<DataT> ShiftLeft(
      vector<DataT> const& bigInt,
      SizeT shift)
    {
      SizeT size = (SizeT)bigInt.size() + shift;
      vector<DataT> result(size);

      Int j = size - 1;

      // Copy
      for(Int i = (Int)bigInt.size() - 1; i >= 0; i--, j--)
      {
        result[j] = bigInt[i];
      }
      
      // Insert zeros
      for(;j >= 0; j--)
      {
        result[j] = 0;
      }

      return result;
    }
   };
}

#endif
