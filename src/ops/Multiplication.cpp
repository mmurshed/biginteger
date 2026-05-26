/**
 * BigMath: BigInteger × BigInteger multiplication.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/ops/Multiplication.h"

namespace BigMath
{
  BigInteger Multiply(BigInteger const &a, BigInteger const &b)
  {
    auto result = Multiply(a.GetInteger(), b.GetInteger(), BigInteger::Base());
    return BigInteger(result, a.IsNegative() != b.IsNegative());
  }

  BigInteger operator*(BigInteger const &a, BigInteger const &b)
  {
    return Multiply(a, b);
  }
}
