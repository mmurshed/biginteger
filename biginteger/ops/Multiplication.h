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
#include "../algorithms/Multiplication.h"

namespace BigMath
{
  static BigInteger Multiply(BigInteger const &a, BigInteger const &b)
  {
    auto result = Multiply(a.GetInteger(), b.GetInteger(), BigInteger::Base());
    return BigInteger(result, a.IsNegative() != b.IsNegative());
  }

  // Multiplies Two BigInteger
  BigInteger operator*(BigInteger const &a, BigInteger const &b)
  {
    return Multiply(a, b);
  }
}

#endif