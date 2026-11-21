/**
 * BigMath: BigInteger limb-shift operators — a << k multiplies by
 * Base^k (whole limbs, not bits).
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
    auto x = a.Limbs();
    if (x.size() + b <= BigInteger::InlineLimbs)
    {
      DataT t[BigInteger::InlineLimbs] = {};
      std::copy(x.begin(), x.end(), t + b);
      return BigInteger::FromLimbs(std::span<const DataT>(t, x.size() + b), a.IsNegative());
    }
    return BigInteger(ShiftLeft(a.AsVector().get(), b), a.IsNegative());
  }

  inline BigInteger operator>>(BigInteger const &a, SizeT b)
  {
    if (b == 0 || a.Zero())
      return a;
    auto x = a.Limbs();
    if (b >= x.size())
      return BigInteger();
    return BigInteger::FromLimbs(x.subspan(b), a.IsNegative());
  }
}

#endif
