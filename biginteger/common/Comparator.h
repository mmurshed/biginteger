/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_COMPARATOR
#define BIGINTEGER_COMPARATOR

#include <span>
#include <vector>
using namespace std;

namespace BigMath
{
  // Span overload — used by div paths that pass subviews of larger vectors.
  inline Int Compare(span<const DataT> a, span<const DataT> b)
  {
    // Trim leading zeros logically by walking down.
    SizeT aLen = (SizeT)a.size();
    while (aLen > 0 && a[aLen - 1] == 0)
      --aLen;
    SizeT bLen = (SizeT)b.size();
    while (bLen > 0 && b[bLen - 1] == 0)
      --bLen;

    if (aLen != bLen)
      return aLen > bLen ? 1 : (aLen < bLen ? -1 : 0);

    for (Int i = (Int)aLen - 1; i >= 0; --i)
      if (a[i] != b[i])
        return a[i] > b[i] ? 1 : -1;
    return 0;
  }

  Int Compare(
      vector<DataT> const &a,
      DataT b)
  {
    // a has more digits than b
    if (a.size() > 1)
      return 1;

    // Case with zero
    bool aIsZero = IsZero(a);
    bool bIsZero = (b == 0);
    if (aIsZero && bIsZero)
      return 0;
    else if (!aIsZero && bIsZero)
      return 1;
    else if (aIsZero && !bIsZero)
      return -1;

    // Both have one digits
    if (a[0] == b)
      return 0;
    return a[0] > b ? 1 : -1;
  }

  // All operations are unsigned
  // Compares this with `with' irrespective of sign
  // Returns
  // 0 if equal
  // +value if this > with
  // -value if this < with
  // Runtime O(n), Space O(1).
  // Trims leading zeros logically (no IsZero pre-walk), then compares by size,
  // then walks down for tie-break.
  Int Compare(
      vector<DataT> const &a,
      vector<DataT> const &b)
  {
    SizeT aLen = (SizeT)a.size();
    while (aLen > 0 && a[aLen - 1] == 0)
      --aLen;
    SizeT bLen = (SizeT)b.size();
    while (bLen > 0 && b[bLen - 1] == 0)
      --bLen;

    if (aLen != bLen)
      return aLen > bLen ? 1 : -1;

    for (Int i = (Int)aLen - 1; i >= 0; --i)
      if (a[i] != b[i])
        return a[i] > b[i] ? 1 : -1;
    return 0;
  }
}

#endif
