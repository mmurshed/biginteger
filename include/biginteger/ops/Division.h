/**
 * BigMath: BigInteger / BigInteger and % operators.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_DIVISION
#define BIGINTEGER_DIVISION

#include <utility>

#include "../BigInteger.h"
#include "../algorithms/Division.h"

namespace BigMath
{
  std::pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, BigInteger const &b);
  BigInteger Divide(BigInteger const &a, BigInteger const &b);
  BigInteger operator/(BigInteger const &a, BigInteger const &b);
  BigInteger operator%(BigInteger const &a, BigInteger const &b);
}

#endif
