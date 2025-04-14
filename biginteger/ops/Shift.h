/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_SHIFT
#define BIGINTEGER_SHIFT

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../algorithms/Shift.h"

namespace BigMath
{
  BigInteger operator<<(BigInteger const &a, SizeT b)
  {
    if (b == 0 || a.Zero())
      return a;
    vector<DataT> result = ShiftLeft(a.GetInteger(), b);
    return BigInteger(result, a.IsNegative());
  }

  BigInteger operator>>(BigInteger const &a, SizeT b)
  {
    if (b == 0 || a.Zero())
      return a;
    vector<DataT> result = ShiftRight(a.GetInteger(), b);
    return BigInteger(result, a.IsNegative());
  }
}

#endif