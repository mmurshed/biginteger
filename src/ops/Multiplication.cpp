/**
 * BigMath: BigInteger × BigInteger multiplication.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/ops/Multiplication.h"
#include "biginteger/ops/SmallOps.h"

#include <utility>

namespace BigMath
{
  BigInteger Multiply(BigInteger const &a, BigInteger const &b)
  {
    if (SmallOps::Fits(a, b))
      return SmallOps::MulMagnitudes(a, b, a.IsNegative() != b.IsNegative());
    auto result = Multiply(a.AsVector().get(), b.AsVector().get(), BigInteger::Base());
    return BigInteger(std::move(result), a.IsNegative() != b.IsNegative());
  }

  BigInteger operator*(BigInteger const &a, BigInteger const &b)
  {
    return Multiply(a, b);
  }
}
