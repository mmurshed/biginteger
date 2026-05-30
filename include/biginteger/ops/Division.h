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
#include "../algorithms/division/ReciprocalDivision.h"

namespace BigMath
{
  class CachedDivision
  {
  private:
    BigInteger divisor;
    ReciprocalDivision::Divider divider;

  public:
    CachedDivision(BigInteger const &b)
        : divisor(b),
          divider(b.GetInteger(), BigInteger::Base())
    {
    }

    std::pair<BigInteger, BigInteger> DivideAndRemainder(
        BigInteger const &a,
        bool computeRemainder = true) const
    {
      auto [qv, rv] = divider.DivideAndRemainder(a.GetInteger(), computeRemainder);
      BigInteger q(qv, a.IsNegative() != divisor.IsNegative());
      BigInteger r(rv, a.IsNegative() || divisor.IsNegative());
      return {q, r};
    }

    BigInteger Divide(BigInteger const &a) const
    {
      auto qv = divider.Divide(a.GetInteger());
      return BigInteger(qv, a.IsNegative() != divisor.IsNegative());
    }

    BigInteger const &Divisor() const
    {
      return divisor;
    }
  };

  inline CachedDivision CacheDivision(BigInteger const &b)
  {
    return CachedDivision(b);
  }

  std::pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, BigInteger const &b);
  BigInteger Divide(BigInteger const &a, BigInteger const &b);
  BigInteger operator/(BigInteger const &a, BigInteger const &b);
  BigInteger operator%(BigInteger const &a, BigInteger const &b);
}

#endif
