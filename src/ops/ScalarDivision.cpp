/**
 * BigMath: BigInteger / DataT and % DataT scalar operators.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/ops/ScalarDivision.h"

#include <stdexcept>

namespace BigMath
{
  BigInteger Divide(BigInteger const &a, DataT b)
  {
    if (b == 0)
      throw std::invalid_argument("Division by zero");

    if (a.Zero())
      return BigInteger();

    std::vector<DataT> q = ClassicDivision::Divide(a.GetInteger(), b, BigInteger::Base());
    BigInteger result(q, false);
    if (a.IsNegative())
      result.SetSign(true);
    return result;
  }

  std::pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, DataT b)
  {
    if (b == 0)
      throw std::invalid_argument("Division by zero");

    if (a.Zero())
      return {BigInteger(), BigInteger()};

    auto result = ClassicDivision::DivideAndRemainder(a.GetInteger(), b, BigInteger::Base());

    BigInteger q(result.first, false);
    BigInteger r(result.second, false);
    if (a.IsNegative())
    {
      q.SetSign(true);
      r.SetSign(true);
    }
    return {q, r};
  }

  BigInteger operator/(BigInteger const &a, DataT const &b) { return Divide(a, b); }
  BigInteger operator%(BigInteger const &a, DataT const &b) { return DivideAndRemainder(a, b).second; }
}
