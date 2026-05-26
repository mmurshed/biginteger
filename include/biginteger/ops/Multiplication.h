/**
 * BigMath: BigInteger × BigInteger multiplication operator.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_MULTIPLICATION
#define BIGINTEGER_MULTIPLICATION

#include "../BigInteger.h"
#include "../algorithms/Multiplication.h"

namespace BigMath
{
  BigInteger Multiply(BigInteger const &a, BigInteger const &b);
  BigInteger operator*(BigInteger const &a, BigInteger const &b);
}

#endif
