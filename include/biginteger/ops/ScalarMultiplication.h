/**
 * BigMath: BigInteger × scalar multiplication.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_SCALAR_MULTIPLICATION
#define BIGINTEGER_SCALAR_MULTIPLICATION

#include <vector>

#include "../BigInteger.h"
#include "../algorithms/multiplication/ClassicMultiplication.h"

namespace BigMath
{
  inline BigInteger Multiply(BigInteger const &a, DataT b)
  {
    if (b == 0 || a.Zero())
      return BigInteger();
    std::vector<DataT> m = ClassicMultiplication::Multiply(a.GetInteger(), b, BigInteger::Base());
    // DataT is unsigned; the b<0 branch was structurally dead. Sign comes solely from a.
    return BigInteger(m, a.IsNegative());
  }

  inline BigInteger operator*(BigInteger const &a, DataT b)
  {
    return Multiply(a, b);
  }
}

#endif
