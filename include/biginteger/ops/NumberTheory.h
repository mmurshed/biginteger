#ifndef BIGINTEGER_NUMBER_THEORY
#define BIGINTEGER_NUMBER_THEORY

#include <utility>

#include "../BigInteger.h"

namespace BigMath
{
  // Return base^exponent mod modulus in [0, modulus).
  // Throws std::invalid_argument if exponent < 0 or modulus <= 0.
  // An exponent of zero returns 1 mod modulus, including when base is zero.
  // Execution time depends on the operands.
  BigInteger PowMod(BigInteger const &base, BigInteger const &exponent,
                    BigInteger const &modulus);

  // Return the nonnegative greatest common divisor. Gcd(0, 0) is zero.
  BigInteger Gcd(BigInteger const &a, BigInteger const &b);

  struct ExtendedGcdResult
  {
    BigInteger gcd;
    BigInteger x;
    BigInteger y;
  };

  // Return gcd, x and y such that a*x + b*y == gcd. Both inputs may be signed.
  // For two zero inputs all three result fields are zero.
  ExtendedGcdResult ExtendedGcd(BigInteger const &a, BigInteger const &b);

  // Return {root, remainder}, where root is floor(sqrt(n)) and
  // remainder is n - root*root. Throws std::invalid_argument for n < 0.
  std::pair<BigInteger, BigInteger> SqrtRem(BigInteger const &n);
}

#endif
