/**
 * BigMath: forced-classic multiplication entry (cross-check only).
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/ops/ClassicMultiplication.h"

#include <utility>

namespace BigMath
{
  BigInteger MultiplyClassic(BigInteger const &a, BigInteger const &b)
  {
    if (a.Zero() || b.Zero())
      return BigInteger();
    auto result = ClassicMultiplication::Multiply(a.AsVector().get(), b.AsVector().get(), BigInteger::Base());
    return BigInteger(std::move(result), a.IsNegative() != b.IsNegative());
  }
}
