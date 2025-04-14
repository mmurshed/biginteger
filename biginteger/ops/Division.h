/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_DIVISION
#define BIGINTEGER_DIVISION

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../algorithms/Division.h"

namespace BigMath
{
  pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, BigInteger const &b)
  {

    auto [qv, rv] = DivideAndRemainder(a.GetInteger(), b.GetInteger(), BigInteger::Base());
    BigInteger q = BigInteger(qv, a.IsNegative() != b.IsNegative());
    BigInteger r = BigInteger(rv, a.IsNegative() || b.IsNegative());
    return {q, r};
  }

  // Divides Two BigInteger
  BigInteger operator/(BigInteger const &a, BigInteger const &b)
  {
    return DivideAndRemainder(a, b).first;
  }

  BigInteger operator%(BigInteger const &a, BigInteger const &b)
  {
    return DivideAndRemainder(a, b).second;
  }

}

#endif