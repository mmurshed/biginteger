/**
 * BigMath: BigInteger / DataT and % DataT scalar operators.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/ops/ScalarDivision.h"

#include <stdexcept>
#include "biginteger/common/Parser.h"
#include "biginteger/ops/SmallOps.h"

namespace BigMath
{
  BigInteger Divide(BigInteger const &a, DataT b)
  {
    if (b == 0)
      throw std::invalid_argument("Division by zero");

    if (a.Zero())
      return BigInteger();

    if (SmallOps::Fits(a))
    {
      DataT q[SmallOps::MaxLimbs];
      auto x = a.Limbs();
      ClassicDivision::DivRem1(q, x.data(), (SizeT)x.size(), b);
      return BigInteger::FromLimbs(std::span<const DataT>(q, x.size()), a.IsNegative());
    }
    return BigInteger(ClassicDivision::Divide(a.AsVector().get(), b, BigInteger::Base()), a.IsNegative());
  }

  std::pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, DataT b)
  {
    if (b == 0)
      throw std::invalid_argument("Division by zero");

    if (a.Zero())
      return {BigInteger(), BigInteger()};

    if (SmallOps::Fits(a))
    {
      DataT q[SmallOps::MaxLimbs];
      auto x = a.Limbs();
      const DataT rem = ClassicDivision::DivRem1(q, x.data(), (SizeT)x.size(), b);
      return {BigInteger::FromLimbs(std::span<const DataT>(q, x.size()), a.IsNegative()),
              BigInteger::FromLimbs(std::span<const DataT>(&rem, 1), a.IsNegative())};
    }
    auto [q, r] = ClassicDivision::DivideAndRemainder(a.AsVector().get(), b, BigInteger::Base());
    return {BigInteger(std::move(q), a.IsNegative()), BigInteger(std::move(r), a.IsNegative())};
  }

  BigInteger operator/(BigInteger const &a, DataT const &b) { return Divide(a, b); }
  BigInteger operator%(BigInteger const &a, DataT const &b)
  {
    if (b == 0)
      throw std::invalid_argument("Division by zero");
    if (SmallOps::Fits(a))
    {
      DataT q[SmallOps::MaxLimbs];
      auto x = a.Limbs();
      const DataT rem = ClassicDivision::DivRem1(q, x.data(), (SizeT)x.size(), b);
      return BigInteger::FromLimbs(std::span<const DataT>(&rem, 1), a.IsNegative());
    }
    DataT rem = ClassicDivision::Remainder(a.AsVector().get(), b, BigInteger::Base());
    return BigInteger(Convert(rem), a.IsNegative());
  }
}
