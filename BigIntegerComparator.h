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
  // All operations are unsigned
  class BigIntegerComparator
  {

    public:

    // Compares this with `with' irrespective of sign
    // Returns
    // 0 if equal
    // +value if this > with
    // -value if this < with
    // Runtime O(n), Space O(1)
    static Int CompareTo(
      vector<DataT> const& a,
      vector<DataT> const& b)
    {
      // Case with zero
      bool aIsZero = BigIntegerUtil::IsZero(a);
      bool bIsZero = BigIntegerUtil::IsZero(b);

      if(aIsZero && bIsZero)
        return 0;
      else if(!aIsZero && bIsZero)
        return 1;
      else if(aIsZero && !bIsZero)
        return -1;

      // Different in size
      Long diff = a.size();
      diff -= b.size();
      if(diff != 0)
        return (Int)diff;

      // Both ints have same number of digits
      Int cmp = 0;
      for(Int i = (Int)a.size() - 1; i >= 0; i--)
      {
        diff = a[i];
        diff -= b[i];
        if(diff != 0)
          return (Int)diff;
      }

      return cmp;
    }
   };
}

#endif
