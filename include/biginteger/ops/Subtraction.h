/**
 * BigMath: BigInteger − BigInteger operator and scalar variant.
 * Sign-aware subtraction; dispatches to add when operand signs differ.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_SUBTRACTION
#define BIGINTEGER_SUBTRACTION

#include "../BigInteger.h"
#include "../common/Builder.h"
#include "../common/Comparator.h"
#include "../algorithms/Addition.h"
#include "../algorithms/Subtraction.h"

namespace BigMath
{
  // |a| − |b| with correct sign; assumes only the magnitudes matter (callers handle signs).
  BigInteger SubtractCompared(BigInteger const &a, BigInteger const &b);
  BigInteger Subtract(BigInteger const &a, BigInteger const &b);
  BigInteger operator-(BigInteger const &a, BigInteger const &b);
  BigInteger operator-(BigInteger const &a, DataT b);
}

#endif
