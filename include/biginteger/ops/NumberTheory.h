#ifndef BIGINTEGER_NUMBER_THEORY
#define BIGINTEGER_NUMBER_THEORY

#include "../BigInteger.h"

namespace BigMath
{
  // Return base^exponent mod modulus in [0, modulus).
  // Throws std::invalid_argument if exponent < 0 or modulus <= 0.
  // An exponent of zero returns 1 mod modulus, including when base is zero.
  // Execution time depends on the operands.
  BigInteger PowMod(BigInteger const &base, BigInteger const &exponent,
                    BigInteger const &modulus);
}

#endif
