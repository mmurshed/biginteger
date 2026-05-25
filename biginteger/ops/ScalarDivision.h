/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_SCALAR_DIVISION
#define BIGINTEGER_SCALAR_DIVISION

#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../algorithms/division/ClassicDivision.h"

namespace BigMath
{
  BigInteger Divide(BigInteger const &a, DataT b)
  {
    if (b == 0)
    {
      throw invalid_argument("Division by zero");
    }

    if (a.Zero())
    {
      return BigInteger(); // case of 0
    }

    // Assume a > b
    vector<DataT> q = ClassicDivision::Divide(
        a.GetInteger(),
        b,
        BigInteger::Base());

    BigInteger result = BigInteger(q, false);

    if (a.IsNegative())
    {
      result.SetSign(true);
    }

    return result;
  }

  pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, DataT b)
  {
    if (b == 0)
    {
      throw invalid_argument("Division by zero");
    }

    if (a.Zero())
    {
      return {BigInteger(), BigInteger()}; // case of 0
    }

    // Assume a > b
    pair<vector<DataT>, vector<DataT>> result = ClassicDivision::DivideAndRemainder(
        a.GetInteger(),
        b,
        BigInteger::Base());

    BigInteger q = BigInteger(result.first, false);
    BigInteger r = BigInteger(result.second, false);

    if (a.IsNegative())
    {
      q.SetSign(true);
      r.SetSign(true);
    }

    return {q, r};
  }

  // Divides Two BigInteger
  BigInteger operator/(BigInteger const &a, DataT const &b)
  {
    return Divide(a, b);
  }

  // Divides Two BigInteger
  BigInteger operator%(BigInteger const &a, DataT const &b)
  {
    return DivideAndRemainder(a, b).second;
  }
}

#endif
