/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_UTIL
#define BIGINTEGER_UTIL

#include <vector>
using namespace std;

#include "Constants.h"

namespace BigMath
{
  // Trims Leading Zeros
  // Runtime O(n), Space O(1)
  SizeT TrimZeros(vector<DataT> &a, Int sizeTo = 0)
  {
    SizeT size = (SizeT)a.size();
    Int i = size;
    while (i > 0 && a[i - 1] == 0 && a.size() > sizeTo)
    {
      a.pop_back();
      i--;
    }
    return size - i; // return how many zeros are removed
  }

  // Runtime O(n), Space O(1)
  SizeT FindNonZeroByte(vector<DataT> const &a, Int start = 0, Int end = -1)
  {
    Int i = (end == -1 ? (Int)a.size() : end + 1);
    while (i > start && a[i - 1] == 0)
      i--;
    return i;
  }

  // Runtime O(n), Space O(1)
  bool IsZero(vector<DataT> const &a, Int start = 0, Int end = -1)
  {
    bool zero = a.size() == 0 || (a.size() == 1 && a[0] == 0);
    return zero ? zero : FindNonZeroByte(a, start, end) == start;
  }

  // Sets some elements of the array to zero.
  void SetBit(vector<DataT> &r, SizeT rStart, SizeT rEnd, DataT bit = 0)
  {
    if (r.size() == 0)
      return;
    rEnd = min(rEnd, (SizeT)(r.size() - 1));
    for (SizeT i = rStart; i <= rEnd; i++)
      r.at(i) = bit;
  }

  // Copy
  void Copy(vector<DataT> const &p, SizeT pStart, SizeT pEnd, vector<DataT> &q, SizeT qStart, SizeT qEnd)
  {
    if (p.size() == 0)
      return;
    for (SizeT i = pStart, j = qStart; i <= pEnd && j <= qEnd; i++, j++)
      q.at(j) = p.at(i);
  }

  void Copy(vector<DataT> const &p, vector<DataT> &q)
  {
    Copy(p, 0, p.size() - 1, q, 0, q.size() - 1);
  }

  void Resize(vector<DataT> &u, SizeT n)
  {
    // Make to size n
    while (u.size() < n)
      u.push_back(0);
  }

  void MakeSameSize(vector<DataT> &u, vector<DataT> &v)
  {
    // Make both same size
    Resize(v, u.size());
    Resize(u, v.size());
  }

  Int Len(SizeT start, SizeT end)
  {
    return (Int)end - (Int)start + 1;
  }
}

#endif
