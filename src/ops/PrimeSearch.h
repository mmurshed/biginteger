#ifndef BIGMATH_PRIME_SEARCH
#define BIGMATH_PRIME_SEARCH

#include <cstddef>

#include "biginteger/BigInteger.h"

namespace BigMath::PrimeSearch
{
  enum class WitnessMode { All, BaseTwo, RandomOnly };

  // The internal witness test expects a positive, odd, normalized value
  // above the word-size range, with small prime factors already screened.
  bool TestLargeCandidate(BigInteger const &mod, int rounds, WitnessMode mode);

  struct SearchStats
  {
    std::size_t screenedCandidates = 0;
  };

  BigInteger FindNext(BigInteger const &n, int rounds, SearchStats *stats = nullptr);
}

#endif
