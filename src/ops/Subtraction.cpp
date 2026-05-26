/**
 * BigMath: Sign-aware BigInteger subtraction.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/ops/Subtraction.h"

namespace BigMath
{
  BigInteger SubtractCompared(BigInteger const &a, BigInteger const &b)
  {
    int cmp = Compare(a.GetInteger(), b.GetInteger());
    if (cmp > 0)
      return BigInteger(
          Subtract(a.GetInteger(), b.GetInteger(), BigInteger::Base()),
          false);
    if (cmp < 0)
      return BigInteger(
          Subtract(b.GetInteger(), a.GetInteger(), BigInteger::Base()),
          true);
    return BigInteger();
  }

  BigInteger Subtract(BigInteger const &a, BigInteger const &b)
  {
    bool aZero = a.Zero();
    bool bZero = b.Zero();

    if (aZero && bZero)
      return BigInteger();

    if (aZero)
    {
      BigInteger result(b);
      result.SetSign(!b.IsNegative());
      return result;
    }

    if (bZero)
      return BigInteger(a);

    bool aNeg = a.IsNegative();
    bool bNeg = b.IsNegative();

    // Opposite signs → magnitude-add, sign from `a`.
    if (aNeg != bNeg)
      return BigInteger(
          Add(a.GetInteger(), b.GetInteger(), BigInteger::Base()),
          aNeg && !bNeg);

    return SubtractCompared(a, b);
  }

  BigInteger operator-(BigInteger const &a, BigInteger const &b)
  {
    return Subtract(a, b);
  }

  BigInteger operator-(BigInteger const &a, DataT b)
  {
    return Subtract(a, BigIntegerBuilder::From(b));
  }
}
