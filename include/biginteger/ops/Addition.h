/**
 * BigMath: BigInteger + BigInteger operator.
 * Sign-aware addition; dispatches to subtract when operand signs differ.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_ADDITION
#define BIGINTEGER_ADDITION

#include "../BigInteger.h"
#include "../algorithms/Addition.h"
#include "Subtraction.h"

namespace BigMath
{
  BigInteger Add(BigInteger const &a, BigInteger const &b);
  BigInteger operator+(BigInteger const &a, BigInteger const &b);
}

#endif
