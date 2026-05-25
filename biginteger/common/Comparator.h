/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_COMPARATOR
#define BIGINTEGER_COMPARATOR

#include <vector>
using namespace std;

namespace BigMath
{
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
  // Runtime O(n), Space O(1)
  Int Compare(
      vector<DataT> const &a,
      vector<DataT> const &b)
  {
    // Case with zero
    bool aIsZero = IsZero(a);
    bool bIsZero = IsZero(b);

    if (aIsZero && bIsZero)
      return 0;
    else if (!aIsZero && bIsZero)
      return 1;
    else if (aIsZero && !bIsZero)
      return -1;

    // Different in size
    Long diff = a.size();
    diff -= b.size();
    if (diff != 0)
      return (Int)diff;

    // Both ints have same number of digits
    Int cmp = 0;
    for (Int i = (Int)a.size() - 1; i >= 0; i--)
    {
      if (a[i] != b[i])
        return a[i] > b[i] ? 1 : -1;
    }

    return cmp;
  }
}

#endif
