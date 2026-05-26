/**
 * BigMath: Limb-vector comparison primitives.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_COMPARATOR
#define BIGINTEGER_COMPARATOR

#include <span>
#include <vector>

#include "Constants.h"
#include "Util.h"

namespace BigMath
{
  // Span overload — used by division paths that pass subviews of larger vectors.
  inline Int Compare(std::span<const DataT> a, std::span<const DataT> b)
  {
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

  inline Int Compare(std::vector<DataT> const &a, DataT b)
  {
    if (a.size() > 1)
      return 1;

    bool aIsZero = IsZero(a);
    bool bIsZero = (b == 0);
    if (aIsZero && bIsZero)
      return 0;
    if (!aIsZero && bIsZero)
      return 1;
    if (aIsZero && !bIsZero)
      return -1;

    if (a[0] == b)
      return 0;
    return a[0] > b ? 1 : -1;
  }

  // Unsigned magnitude compare. Returns 0/+/− for ==, >, <.
  // O(n) time, O(1) space. Trims leading zeros logically by walking down.
  inline Int Compare(std::vector<DataT> const &a, std::vector<DataT> const &b)
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
