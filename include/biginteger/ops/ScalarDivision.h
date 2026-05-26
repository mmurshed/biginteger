/**
 * BigMath: BigInteger / DataT and % DataT scalar operators.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_SCALAR_DIVISION
#define BIGINTEGER_SCALAR_DIVISION

#include <utility>

#include "../BigInteger.h"
#include "../algorithms/division/ClassicDivision.h"

namespace BigMath
{
  BigInteger Divide(BigInteger const &a, DataT b);
  std::pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, DataT b);
  BigInteger operator/(BigInteger const &a, DataT const &b);
  BigInteger operator%(BigInteger const &a, DataT const &b);
}

#endif
