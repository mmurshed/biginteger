/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_COMPARISON
#define BIGINTEGER_COMPARISON

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../common/Comparator.h"

namespace BigMath
{
  // Comparison operators
  bool operator==(BigInteger const &a, BigInteger const &b)
  {
    return a.CompareTo(b) == 0;
  }

  bool operator!=(BigInteger const &a, BigInteger const &b)
  {
    return a.CompareTo(b) != 0;
  }

  bool operator>=(BigInteger const &a, BigInteger const &b)
  {
    return a.CompareTo(b) >= 0;
  }

  bool operator<=(BigInteger const &a, BigInteger const &b)
  {
    return a.CompareTo(b) <= 0;
  }

  bool operator>(BigInteger const &a, BigInteger const &b)
  {
    return a.CompareTo(b) > 0;
  }

  bool operator<(BigInteger const &a, BigInteger const &b)
  {
    return a.CompareTo(b) < 0;
  }
}

#endif