/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_SCALAR_DIVISION
#define BIGINTEGER_SCALAR_DIVISION

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../algorithms/classic/ClassicDivision.h"

namespace BigMath
{
  class BigIntegerScalarDivision
  {
  public:
    static BigInteger Divide(BigInteger const &a, DataT b)
    {
      if (a.IsZero() || b == 0)
      {
        return BigInteger(); // case of 0
      }

      // Assume a > b
      vector<DataT> q = ClassicDivision::Divide(
          a.GetInteger(),
          b,
          BigInteger::Base());

      BigInteger result = BigInteger(q, false);

      if (a.IsNegative() != b < 0)
      {
        result.SetSign(true);
      }

      return result;
    }

    static pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, DataT b)
    {
      if (a.IsZero() || b == 0)
      {
        return {BigInteger(), BigInteger()}; // case of 0
      }

      // Assume a > b
      pair< vector<DataT>, vector<DataT> > result = ClassicDivision::DivideAndRemainder(
          a.GetInteger(),
          b,
          BigInteger::Base());

      BigInteger q = BigInteger(result.first, false);
      BigInteger r = BigInteger(result.second, false);

      if (a.IsNegative() != b < 0)
      {
        q.SetSign(true);
        r.SetSign(true);
      }

      return {q, r};
    }
  };

  // Divides Two BigInteger
  BigInteger operator/(BigInteger const &a, DataT const &b)
  {
    return BigIntegerScalarDivision::Divide(a, b);
  }

  // Divides Two BigInteger
  BigInteger operator%(BigInteger const &a, DataT const &b)
  {
    return BigIntegerScalarDivision::DivideAndRemainder(a, b).second;
  }
}

#endif