/**
 * BigMath: BigInteger / BigInteger and % operators.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/ops/Division.h"
#include "biginteger/ops/SmallOps.h"

namespace BigMath
{
  std::pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, BigInteger const &b)
  {
    if (SmallOps::Fits(a, b))
    {
      std::pair<BigInteger, BigInteger> qr;
      SmallOps::DivRem(a, b, &qr.first, &qr.second);
      return qr;
    }
    auto [qv, rv] = detail::DivideAndRemainderLimbs(a.Limbs(), b.Limbs(), BigInteger::Base(), true);
    BigInteger q(std::move(qv), a.IsNegative() != b.IsNegative());
    // Truncated division: remainder takes the dividend's sign.
    BigInteger r(std::move(rv), a.IsNegative());
    return {std::move(q), std::move(r)};
  }

  BigInteger Divide(BigInteger const &a, BigInteger const &b)
  {
    if (SmallOps::Fits(a, b))
    {
      BigInteger q;
      SmallOps::DivRem(a, b, &q, nullptr);
      return q;
    }
    auto qv = detail::DivideAndRemainderLimbs(a.Limbs(), b.Limbs(), BigInteger::Base(), false).first;
    return BigInteger(std::move(qv), a.IsNegative() != b.IsNegative());
  }

  BigInteger operator/(BigInteger const &a, BigInteger const &b)
  {
    return Divide(a, b);
  }

  BigInteger operator%(BigInteger const &a, BigInteger const &b)
  {
    if (SmallOps::Fits(a, b))
    {
      BigInteger r;
      SmallOps::DivRem(a, b, nullptr, &r);
      return r;
    }
    return DivideAndRemainder(a, b).second;
  }
}
