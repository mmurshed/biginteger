/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_SCALAR_MULTIPLICATION
#define BIGINTEGER_SCALAR_MULTIPLICATION

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../algorithms/classic/ClassicMultiplication.h"

namespace BigMath
{
  class BigIntegerScalarMultiplication
  {
  public:
    static BigInteger Multiply(BigInteger const &a, DataT b)
    {
      if (a.IsZero() || b == 0)
        return BigInteger(); // 0 times anything is zero

      vector<DataT> m = ClassicMultiplication::Multiply(a.GetInteger(), b, BigInteger::Base());
      BigInteger result = BigInteger(m, false);
      if (a.IsNegative() != b < 0)
        result.SetSign(true);

      return result;
    }
  };

  // Multiplies Two BigInteger
  BigInteger operator*(BigInteger const &a, DataT b)
  {
    return BigIntegerScalarMultiplication::Multiply(a, b);
  }

}

#endif