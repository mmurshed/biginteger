/**
 * BigMath: BigInteger comparison operators.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_COMPARISON
#define BIGINTEGER_COMPARISON

#include "../BigInteger.h"

namespace BigMath
{
  inline bool operator==(BigInteger const &a, BigInteger const &b) { return a.CompareTo(b) == 0; }
  inline bool operator!=(BigInteger const &a, BigInteger const &b) { return a.CompareTo(b) != 0; }
  inline bool operator>=(BigInteger const &a, BigInteger const &b) { return a.CompareTo(b) >= 0; }
  inline bool operator<=(BigInteger const &a, BigInteger const &b) { return a.CompareTo(b) <= 0; }
  inline bool operator>(BigInteger const &a, BigInteger const &b)  { return a.CompareTo(b) >  0; }
  inline bool operator<(BigInteger const &a, BigInteger const &b)  { return a.CompareTo(b) <  0; }
}

#endif
