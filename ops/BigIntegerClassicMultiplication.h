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
#include "../algorithms/classic/ClassicMultiplication.h"

namespace BigMath
{
  class BigIntegerClassicMultiplication
  {
  public:
    static BigInteger Multiply(BigInteger const &a, BigInteger const &b)
    {
      if (a.IsZero() || b.IsZero())
        return BigInteger(); // 0 times anything is zero

      BigInteger result = ClassicMultiplication::Multiply(a.GetInteger(), b.GetInteger(), BigInteger::Base());
      if (a.IsNegative() != b.IsNegative())
        result.SetSign(true);

      return result;
    }
  };
}

#endif