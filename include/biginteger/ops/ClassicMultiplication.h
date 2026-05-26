/**
 * BigMath: BigInteger × BigInteger forced through ClassicMultiplication.
 * Reserved for cross-checking tests; production code uses the multiplication
 * dispatcher in `algorithms/Multiplication.h`.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_CLASSIC_MULTIPLICATION
#define BIGINTEGER_CLASSIC_MULTIPLICATION

#include "../BigInteger.h"
#include "../algorithms/multiplication/ClassicMultiplication.h"

namespace BigMath
{
  BigInteger MultiplyClassic(BigInteger const &a, BigInteger const &b);
}

#endif
