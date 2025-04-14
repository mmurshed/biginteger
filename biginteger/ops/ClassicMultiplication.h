/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_CLASSIC_MULTIPLICATION
#define BIGINTEGER_CLASSIC_MULTIPLICATION

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../algorithms/multiplication/ClassicMultiplication.h"

namespace BigMath
{
  BigInteger MultiplyClassic(BigInteger const &a, BigInteger const &b)
  {
    if (a.Zero() || b.Zero())
      return BigInteger(); // 0 times anything is zero

    auto result = ClassicMultiplication::Multiply(a.GetInteger(), b.GetInteger(), BigInteger::Base());
    return BigInteger(result, a.IsNegative() != b.IsNegative());
  }

  BigInteger MultiplyClassic(BigInteger const &a, DataT b)
  {
    if (b == 0 || a.Zero())
      return BigInteger(); // 0 times anything is zero

    auto result = ClassicMultiplication::Multiply(a.GetInteger(), b, BigInteger::Base());
    return BigInteger(result, a.IsNegative() != b<0);
  }
}

#endif