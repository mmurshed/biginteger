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

    // Base 10
    static const ULong Base10 = 10;
    static const ULong Base100 = 100;

    public:
    // Trims Leading Zeros
    // Runtime O(n), Space O(1)
    static SizeT TrimZeros(vector<DataT>& a)
    {
      SizeT size = a.size();
      Int i = size;
      while(i > 0 && a[i - 1] == 0)
      {
        a.pop_back();
        i--;
      }
      return size - i; // return how many zeros are removed
    }

    // Runtime O(n), Space O(1)
    static SizeT FindNonZeroByte(vector<DataT> const& a)
    {
      Int i = a.size();
      while(i > 0 && a[i - 1] == 0)
        i--;
      return i;
    }
    
    // Runtime O(n), Space O(1)
    static bool IsZero(vector<DataT> const& a)
    {
      bool zero = a.size() == 0 || (a.size() == 1 && a[0] == 0);
      if(!zero)
      {
        Int i = a.size();
        while(i > 0 && a[i - 1] == 0)
          i--;
        zero = (i == 0);
      }
      return zero;
    }
   };
}

#endif