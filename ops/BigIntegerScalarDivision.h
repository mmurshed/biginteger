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
  };

  // Divides Two BigInteger
  BigInteger operator/(BigInteger const &a, DataT const &b)
  {
    return BigIntegerScalarDivision::Divide(a, b);
  }
}

#endif