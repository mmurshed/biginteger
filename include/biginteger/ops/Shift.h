/**
 * BigMath: BigInteger bit-shift operators.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_SHIFT
#define BIGINTEGER_SHIFT

#include "../BigInteger.h"
#include "../algorithms/Shift.h"

namespace BigMath
{
  inline BigInteger operator<<(BigInteger const &a, SizeT b)
  {
    if (b == 0 || a.Zero())
      return a;
    return BigInteger(ShiftLeft(a.GetInteger(), b), a.IsNegative());
  }

  inline BigInteger operator>>(BigInteger const &a, SizeT b)
  {
    if (b == 0 || a.Zero())
      return a;
    return BigInteger(ShiftRight(a.GetInteger(), b), a.IsNegative());
  }
}

#endif
