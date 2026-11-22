#include "biginteger/ops/NumberTheory.h"

#include <stdexcept>

#include "PrimeSearch.h"

namespace BigMath
{
  BigInteger NextProbablePrime(BigInteger const &n, int rounds)
  {
    if (rounds <= 0)
      throw std::invalid_argument("NextProbablePrime requires a positive round count");
    return PrimeSearch::FindNext(n, [rounds](BigInteger const &candidate) {
      return IsProbablePrime(candidate, rounds);
    });
  }
}
