/**
 * BigMath: BigInteger / DataT and % DataT scalar operators.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/ops/ScalarDivision.h"

#include <stdexcept>
#include "biginteger/common/Parser.h"

namespace BigMath
{
  BigInteger Divide(BigInteger const &a, DataT b)
  {
    if (b == 0)
      throw std::invalid_argument("Division by zero");

    if (a.Zero())
      return BigInteger();

    return BigInteger(ClassicDivision::Divide(a.GetInteger(), b, BigInteger::Base()), a.IsNegative());
  }

  std::pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, DataT b)
  {
    if (b == 0)
      throw std::invalid_argument("Division by zero");

    if (a.Zero())
      return {BigInteger(), BigInteger()};

    auto [q, r] = ClassicDivision::DivideAndRemainder(a.GetInteger(), b, BigInteger::Base());
    return {BigInteger(std::move(q), a.IsNegative()), BigInteger(std::move(r), a.IsNegative())};
  }

  BigInteger operator/(BigInteger const &a, DataT const &b) { return Divide(a, b); }
  BigInteger operator%(BigInteger const &a, DataT const &b)
  {
    DataT rem = ClassicDivision::Remainder(a.GetInteger(), b, BigInteger::Base());
    return BigInteger(Convert(rem), a.IsNegative());
  }
}
