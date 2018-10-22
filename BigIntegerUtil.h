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
    // Base 10
    static const ULong Base10 = 10;
    // Base 100
    static const ULong Base100 = 100;
    // The Base Used
    static const ULong Base100M = 100000000lu;
    // Number of digits in `BASE'
    static const SizeT Base100M_Zeroes = 8;

    // Base 2^32
    static const ULong Base2_32 = 4294967296ul; // 2^32

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
    static SizeT FindNonZeroByte(vector<DataT> const& a, Int start = -1, Int end = -1)
    {
      start = (start == -1 ? 0 : start);
      Int i = (end == -1 ? a.size() : end + 1);
      while(i > start && a[i - 1] == 0)
        i--;
      return i;
    }
    
    // Runtime O(n), Space O(1)
    static bool IsZero(vector<DataT> const& a, Int start = -1, Int end = -1)
    {
      bool zero = a.size() == 0 || (a.size() == 1 && a[0] == 0);
      return zero ? zero : FindNonZeroByte(a, start, end) == 0;
    }

    // Sets some elements of the array to zero.
    static void SetBit(vector<DataT>& r, SizeT rStart, SizeT rEnd, DataT bit = 0)
    {
      rEnd = min(rEnd, (SizeT)(r.size() - 1));
      for(SizeT i = rStart; i <= rEnd; i++)
        r.at(i) = bit;
    }

    static void MakeSameSize(vector<DataT>& u, vector<DataT>& v)
    {
      // Make both same size
      while(v.size() < u.size())
        v.push_back(0);
      while(u.size() < v.size())
        u.push_back(0);
    }
   };
}

#endif