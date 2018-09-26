/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIG_INTEGER_UTIL_H
#define BIG_INTEGER_UTIL_H

#include <vector>
using namespace std;

#include "BigInteger.h"

namespace BigMath
{
  typedef long long Long;
  typedef unsigned long long ULong;
  
  typedef int Int;
  typedef unsigned int UInt;

  typedef UInt DataT;
  typedef UInt SizeT;

  class BigIntegerUtil
  {
    public:
    // The Base Used
    static const ULong Base10n = 100000000lu;
    // Number of digits in `BASE'
    static const SizeT Base10n_Digit = 8;

    // Base 2^32
    static const ULong Base2_32 = 4294967296ul; // 2^32

    public:
    // Trims Leading Zeros
    static void TrimZeros(vector<DataT>& aInt)
    {
      while(aInt.size() > 0 && aInt[aInt.size() - 1] == 0)
      {
        aInt.pop_back();
      }
    }
    
    static bool IsZero(vector<DataT> const& bigInt)
    {
      return bigInt.size() == 0 || (bigInt.size() == 1 && bigInt[0] == 0);
    }
   };
}

#endif