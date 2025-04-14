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
#include "../algorithms/multiplication/ClassicMultiplication.h"

namespace BigMath
{
  BigInteger Multiply(BigInteger const &a, DataT b)
  {
    if (b == 0 || a.Zero())
      return BigInteger(); // 0 times anything is zero

    vector<DataT> m = ClassicMultiplication::Multiply(a.GetInteger(), b, BigInteger::Base());
    return BigInteger(m, a.IsNegative() != b < 0);
  }

  // Multiplies Two BigInteger
  BigInteger operator*(BigInteger const &a, DataT b)
  {
    return Multiply(a, b);
  }

}

#endif