/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_UTIL
#define BIGINTEGER_UTIL

#include <cstring>
#include <span>
#include <vector>
using namespace std;

#include "Constants.h"

namespace BigMath
{
  // compute the greatest common divisor (GCD)
  DataT gcd(DataT a, DataT b)
  {
    while (b != 0)
    {
      DataT temp = b;
      b = a % b;
      a = temp;
    }
    return a;
  }

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

  // Canonical trim for values that cross algorithm/API boundaries: zero is {0}.
  // Internal kernels may still use an empty vector as a temporary zero to avoid work.
  inline SizeT TrimZerosToOne(vector<DataT> &a)
  {
    SizeT removed = TrimZeros(a);
    if (a.empty())
      a.push_back(0);
    return removed;
  }

  inline vector<DataT> CanonicalZero()
  {
    return vector<DataT>{0};
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

  // Span overload — used by div paths that pass subviews without materializing vectors.
  inline bool IsZero(span<const DataT> a)
  {
    if (a.empty())
      return true;
    for (DataT x : a)
      if (x)
        return false;
    return true;
  }

  // Sets some elements of the array to zero.
  void SetBit(vector<DataT> &r, SizeT rStart, SizeT rEnd, DataT bit = 0)
  {
    if (r.size() == 0)
      return;
    rEnd = min(rEnd, (SizeT)(r.size() - 1));
    for (SizeT i = rStart; i <= rEnd; i++)
      r[i] = bit;
  }

  // Copy
  void Copy(vector<DataT> const &p, SizeT pStart, SizeT pEnd, vector<DataT> &q, SizeT qStart, SizeT qEnd)
  {
    if (p.size() == 0)
      return;
    SizeT count = min(pEnd - pStart + 1, qEnd - qStart + 1);
    std::memcpy(q.data() + qStart, p.data() + pStart, count * sizeof(DataT));
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
