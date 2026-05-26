/**
 * BigMath: BigInteger / BigInteger and % operators.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/ops/Division.h"

namespace BigMath
{
  std::pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, BigInteger const &b)
  {
    auto [qv, rv] = DivideAndRemainder(a.GetInteger(), b.GetInteger(), BigInteger::Base());
    BigInteger q(qv, a.IsNegative() != b.IsNegative());
    BigInteger r(rv, a.IsNegative() || b.IsNegative());
    return {q, r};
  }

  BigInteger Divide(BigInteger const &a, BigInteger const &b)
  {
    auto qv = Divide(a.GetInteger(), b.GetInteger(), BigInteger::Base());
    return BigInteger(qv, a.IsNegative() != b.IsNegative());
  }

  BigInteger operator/(BigInteger const &a, BigInteger const &b)
  {
    return Divide(a, b);
  }

  BigInteger operator%(BigInteger const &a, BigInteger const &b)
  {
    return DivideAndRemainder(a, b).second;
  }
}
