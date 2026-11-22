/**
 * BigMath: BigInteger comparison operators.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_COMPARISON
#define BIGINTEGER_COMPARISON

#include <cstring>

#include "../BigInteger.h"

namespace BigMath
{
  inline bool operator==(BigInteger const &a, BigInteger const &b)
  {
    if (&a == &b)
      return true;
    if (a.IsNegative() != b.IsNegative())
      return false;

    const auto lhs = a.Limbs();
    const auto rhs = b.Limbs();
    if (lhs.size() != rhs.size())
    {
      // The size constructor can leave zero padded to the requested width.
      const auto longer = lhs.size() > rhs.size() ? lhs : rhs;
      return longer.back() == 0 && Compare(lhs, rhs) == 0;
    }
    if (lhs.back() != rhs.back() || lhs.front() != rhs.front())
      return false;

    // Check the low block first, then walk down from the high end so that a
    // difference just below the top limb doesn't require scanning the value.
    constexpr std::size_t blockLimbs = 8;
    constexpr std::size_t blockBytes = blockLimbs * sizeof(DataT);
    if (lhs.size() <= blockLimbs)
      return lhs.size() <= 2 ||
             std::memcmp(lhs.data() + 1, rhs.data() + 1, (lhs.size() - 2) * sizeof(DataT)) == 0;
    if (std::memcmp(lhs.data(), rhs.data(), blockBytes) != 0)
      return false;
    for (auto end = lhs.size(); end > blockLimbs; end -= blockLimbs)
      if (std::memcmp(lhs.data() + (end - blockLimbs), rhs.data() + (end - blockLimbs), blockBytes) != 0)
        return false;
    return true;
  }

  inline bool operator!=(BigInteger const &a, BigInteger const &b) { return !(a == b); }
  inline bool operator>=(BigInteger const &a, BigInteger const &b) { return a.CompareTo(b) >= 0; }
  inline bool operator<=(BigInteger const &a, BigInteger const &b) { return a.CompareTo(b) <= 0; }
  inline bool operator>(BigInteger const &a, BigInteger const &b)  { return a.CompareTo(b) >  0; }
  inline bool operator<(BigInteger const &a, BigInteger const &b)  { return a.CompareTo(b) <  0; }
}

#endif
