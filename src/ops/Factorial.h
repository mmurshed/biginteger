#ifndef BIGMATH_FACTORIAL_DETAIL
#define BIGMATH_FACTORIAL_DETAIL

#include "biginteger/ops/NumberTheory.h"

namespace BigMath::detail
{
  // Packed product tree, also used as a benchmark baseline for prime swing.
  BigInteger FactorialProductTree(UInt n);
}

#endif
