/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_MULTIPLICATION
#define BIGINTEGER_MULTIPLICATION

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "MultiplicationStrategy.h"

namespace BigMath
{
  class BigIntegerMultiplication
  {
  public:
    static BigInteger Multiply(BigInteger const &a, BigInteger const &b)
    {
      if (a.IsZero() || b.IsZero())
        return BigInteger(); // 0 times anything is zero

      BigInteger result = MultiplicationStrategy::MultiplyUnsigned(a.GetInteger(), b.GetInteger(), BigInteger::Base());
      if (a.IsNegative() != b.IsNegative())
        result.SetSign(true);

      return result;
    }
  };

  // Multiplies Two BigInteger
  BigInteger operator*(BigInteger const &a, BigInteger const &b)
  {
    return BigIntegerMultiplication::Multiply(a, b);
  }
}

#endif