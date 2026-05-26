/**
 * BigMath: forced-classic multiplication entry (cross-check only).
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/ops/ClassicMultiplication.h"

namespace BigMath
{
  BigInteger MultiplyClassic(BigInteger const &a, BigInteger const &b)
  {
    if (a.Zero() || b.Zero())
      return BigInteger();
    auto result = ClassicMultiplication::Multiply(a.GetInteger(), b.GetInteger(), BigInteger::Base());
    return BigInteger(result, a.IsNegative() != b.IsNegative());
  }
}
