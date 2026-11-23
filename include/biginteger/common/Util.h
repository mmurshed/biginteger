/**
 * BigMath: Common limb-level utilities.
 *
 * All functions are header-inline: they are tiny helpers on raw limb vectors
 * called from hot paths and the compiler benefits from inlining them at every
 * call site.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_UTIL
#define BIGINTEGER_UTIL

#include <cstddef>
#include <span>
#include <vector>
#include <algorithm>

#include "Constants.h"

namespace BigMath
{
  namespace LimbDetail
  {
    // Count with high zero padding ignored, keeping one limb for a nonempty zero.
    inline SizeT TrimmedLimbCount(std::span<const DataT> limbs)
    {
      auto size = static_cast<SizeT>(limbs.size());
      while (size > 1 && limbs[size - 1] == 0)
        --size;
      return size;
    }
  }

  // GCD on single limbs (Euclidean).
  inline DataT gcd(DataT a, DataT b)
  {
    while (b != 0)
    {
      DataT temp = b;
      b = a % b;
      a = temp;
    }
    return a;
  }

  // Strips trailing zeros from `a` until either it's empty or the high limb is non-zero.
  // Returns the number of limbs removed. Pass sizeTo > 0 to leave at least that many limbs.
  inline SizeT TrimZeros(std::vector<DataT> &a, Int sizeTo = 0)
  {
    SizeT size = (SizeT)a.size();
    Int i = size;
    while (i > 0 && a[i - 1] == 0 && (Int)a.size() > sizeTo)
    {
      a.pop_back();
      i--;
    }
    return size - i;
  }

  // Canonical trim for values crossing the algorithm/API boundary: zero is {0}, not {}.
  // Internal kernels may use an empty vector as a temporary zero; public values must not.
  inline SizeT TrimZerosToOne(std::vector<DataT> &a)
  {
    SizeT removed = TrimZeros(a);
    if (a.empty())
      a.push_back(0);
    return removed;
  }

  inline Int FindNonZeroByte(std::vector<DataT> const &a, Int start = 0, Int end = -1)
  {
    Int i = (end == -1 ? (Int)a.size() : end + 1);
    while (i > start && a[i - 1] == 0)
      i--;
    return i;
  }

  inline bool IsZero(std::vector<DataT> const &a, Int start = 0, Int end = -1)
  {
    bool zero = a.size() == 0 || (a.size() == 1 && a[0] == 0);
    return zero ? zero : FindNonZeroByte(a, start, end) == start;
  }

  // Span overload used by division paths that hand in subviews without materializing vectors.
  inline bool IsZero(std::span<const DataT> a)
  {
    if (a.empty())
      return true;
    for (DataT x : a)
      if (x)
        return false;
    return true;
  }

  inline void SetBit(std::vector<DataT> &r, SizeT rStart, SizeT rEnd, DataT bit = 0)
  {
    if (r.size() == 0)
      return;
    rEnd = std::min(rEnd, (SizeT)(r.size() - 1));
    for (SizeT i = rStart; i <= rEnd; i++)
      r[i] = bit;
  }

  inline Int Len(SizeT start, SizeT end)
  {
    return (Int)end - (Int)start + 1;
  }
}

#endif
