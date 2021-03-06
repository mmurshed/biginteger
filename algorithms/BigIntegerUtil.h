/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIG_INTEGER_UTIL_H
#define BIG_INTEGER_UTIL_H

#include <vector>
using namespace std;

#include "../BigInteger.h"

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
    // Base 10^8
    static const ULong Base100M = 100000000lu;
    // Number of digits in `BASE'
    static const SizeT Base100M_Zeroes = 8;

    // Base 2
    static const ULong Base2 = 2;
    // Base 2^16
    static const ULong Base2_16 = 65536ul;
    // Base 2^32
    static const ULong Base2_32 = 4294967296ul; // 2^32

    public:
    // Trims Leading Zeros
    // Runtime O(n), Space O(1)
    static SizeT TrimZeros(vector<DataT>& a, Int sizeTo = 0)
    {
      SizeT size = (SizeT)a.size();
      Int i = size;
      while(i > 0 && a[i - 1] == 0 && a.size() > sizeTo)
      {
        a.pop_back();
        i--;
      }
      return size - i; // return how many zeros are removed
    }

    // Runtime O(n), Space O(1)
    static SizeT FindNonZeroByte(vector<DataT> const& a, Int start = 0, Int end = -1)
    {
      Int i = (end == -1 ? (Int)a.size() : end + 1);
      while(i > start && a[i - 1] == 0)
        i--;
      return i;
    }
    
    // Runtime O(n), Space O(1)
    static bool IsZero(vector<DataT> const& a, Int start = 0, Int end = -1)
    {
      bool zero = a.size() == 0 || (a.size() == 1 && a[0] == 0);
      return zero ? zero : FindNonZeroByte(a, start, end) == start;
    }

    // Sets some elements of the array to zero.
    static void SetBit(vector<DataT>& r, SizeT rStart, SizeT rEnd, DataT bit = 0)
    {
      rEnd = min(rEnd, (SizeT)(r.size() - 1));
      for(SizeT i = rStart; i <= rEnd; i++)
        r.at(i) = bit;
    }

    // Copy
    static void Copy(vector<DataT> const& p, vector<DataT>& q)
    {
        Copy(p, 0, p.size() - 1, q, 0, q.size() - 1);
    }

    static void Copy(vector<DataT> const& p, SizeT pStart, SizeT pEnd, vector<DataT>& q, SizeT qStart, SizeT qEnd)
    {
      if(p.size() == 0)
        return;
      for(SizeT i = pStart, j = qStart; i <= pEnd && j <= qEnd; i++, j++)
        q.at(j) = p.at(i);
    }

    static void MakeSameSize(vector<DataT>& u, vector<DataT>& v)
    {
      // Make both same size
      Resize(v, u.size());
      Resize(u, v.size());
    }

    static void Resize(vector<DataT>& u, SizeT n)
    {
      // Make to size n
      while(u.size() < n)
        u.push_back(0);
    }

    static Int Len(SizeT start, SizeT end)
    {
      return (Int) end - (Int) start + 1;
    }
   };
}

#endif
