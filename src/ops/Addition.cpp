/**
 * BigMath: Sign-aware BigInteger addition.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/ops/Addition.h"
#include "biginteger/ops/Subtraction.h"

namespace BigMath
{
  BigInteger Add(BigInteger const &a, BigInteger const &b)
  {
    bool aZero = a.Zero();
    bool bZero = b.Zero();

    if (aZero && bZero)
      return BigInteger();
    if (aZero)
      return BigInteger(b);
    if (bZero)
      return BigInteger(a);

    bool aNeg = a.IsNegative();
    bool bNeg = b.IsNegative();

    // Opposite signs → reduce to subtraction. Branches kept distinct rather than
    // collapsed because the call argument order matters for sign correctness.
    if (aNeg && !bNeg)
      return SubtractCompared(b, a);
    if (!aNeg && bNeg)
      return SubtractCompared(a, b);

    // Same sign: magnitude-add, preserve sign.
    return BigInteger(
        Add(a.GetInteger(), b.GetInteger(), BigInteger::Base()),
        aNeg && bNeg);
  }

  BigInteger operator+(BigInteger const &a, BigInteger const &b)
  {
    return Add(a, b);
  }
}
