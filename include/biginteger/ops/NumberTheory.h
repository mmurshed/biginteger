#ifndef BIGINTEGER_NUMBER_THEORY
#define BIGINTEGER_NUMBER_THEORY

#include <utility>

#include "../BigInteger.h"
#include "ChineseRemainder.h"

namespace BigMath
{
  // Return base^exponent mod modulus in [0, modulus).
  // Throws std::invalid_argument if exponent < 0 or modulus <= 0.
  // An exponent of zero returns 1 mod modulus, including when base is zero.
  // Execution time depends on the operands.
  BigInteger PowMod(BigInteger const &base, BigInteger const &exponent,
                    BigInteger const &modulus);

  // Exact for values below 2^64. For larger values, a true result means
  // the input passed Miller–Rabin with the requested number of random bases.
  // Returns false for n < 2. Throws std::invalid_argument if rounds <= 0.
  bool IsProbablePrime(BigInteger const &n, int rounds = 32);

  // Search strictly above n, returning 2 for n < 2. Results below 2^64
  // are exact. Larger candidates use the same test as IsProbablePrime.
  // Throws std::invalid_argument if rounds <= 0.
  BigInteger NextProbablePrime(BigInteger const &n, int rounds = 32);

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

  // Return {root, n - root^degree}, rounding the real root toward zero.
  // Throws std::invalid_argument for zero degree or a negative input with even degree.
  std::pair<BigInteger, BigInteger> RootRem(BigInteger const &n, UInt degree);

  // Return n!, with 0! equal to one.
  BigInteger Factorial(UInt n);

  // Return the binomial coefficient n choose k, or zero when k > n.
  BigInteger Binomial(UInt n, UInt k);
}

#endif
