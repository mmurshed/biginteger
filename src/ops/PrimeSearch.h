#ifndef BIGMATH_PRIME_SEARCH
#define BIGMATH_PRIME_SEARCH

#include <array>
#include <cstdint>

#include "biginteger/ops/Operations.h"

namespace BigMath::PrimeSearch
{
  inline constexpr unsigned SegmentSize = 1024;
  inline constexpr unsigned DefaultSieveLimit = 4096;
  inline constexpr unsigned LargeSieveLimit = 32768;

  inline constexpr auto SmallComposites = [] {
    std::array<bool, LargeSieveLimit> composite{};
    for (unsigned prime = 3; prime * prime < LargeSieveLimit; prime += 2)
      if (!composite[prime])
        for (unsigned multiple = prime * prime; multiple < LargeSieveLimit; multiple += 2 * prime)
          composite[multiple] = true;
    return composite;
  }();

  inline constexpr auto SmallPrimes = [] {
    constexpr auto count = [] {
      unsigned total = 0;
      for (unsigned prime = 3; prime < LargeSieveLimit; prime += 2)
        total += !SmallComposites[prime];
      return total;
    }();
    std::array<uint16_t, count> primes{};
    unsigned i = 0;
    for (unsigned prime = 3; prime < LargeSieveLimit; prime += 2)
      if (!SmallComposites[prime])
        primes[i++] = static_cast<uint16_t>(prime);
    return primes;
  }();

  inline BigInteger FromWord(DataT value)
  {
    return BigInteger::FromLimbs({&value, 1}, false);
  }

  // The predicate sees the candidates that survived the sieve, in order.
  template <class IsPrime>
  BigInteger FindNext(BigInteger const &n, IsPrime test)
  {
    const auto one = FromWord(1);
    const auto two = FromWord(2);
    const auto start = BigInteger::FromLimbs(n.Limbs(), n.IsNegative());
    if (start < two)
      return two;

    auto first = start + ((start[0] & 1) ? two : one);
    if (start.size() <= 64 / LimbBits)
    {
      // Word-size primality tests are cheap and also handle primes in the
      // sieve's own table, which must not be crossed out as multiples.
      while (!test(first))
        first = first + two;
      return first;
    }

    // At larger sizes one rejected Miller–Rabin test pays for more sieving.
    const unsigned limit = start.size() > 512 / LimbBits ? LargeSieveLimit : DefaultSieveLimit;
    const auto primeCount = static_cast<size_t>(
        std::lower_bound(SmallPrimes.begin(), SmallPrimes.end(), limit) - SmallPrimes.begin());
    std::array<uint16_t, SmallPrimes.size()> offsets;
    for (size_t i = 0; i < primeCount; ++i)
    {
      const DataT p = SmallPrimes[i];
      const DataT remainder = (first % p)[0];
      const DataT distance = remainder == 0 ? 0 : p - remainder;
      offsets[i] = static_cast<uint16_t>((distance + ((distance & 1) ? p : 0)) / 2);
    }

    const auto stride = FromWord(2 * SegmentSize);
    for (;;)
    {
      std::array<bool, SegmentSize> composite{};
      for (size_t i = 0; i < primeCount; ++i)
      {
        unsigned offset = offsets[i];
        for (; offset < SegmentSize; offset += SmallPrimes[i])
          composite[offset] = true;
        // Carry the first multiple beyond this segment into the next one.
        offsets[i] = static_cast<uint16_t>(offset - SegmentSize);
      }
      for (unsigned i = 0; i < SegmentSize; ++i)
        if (!composite[i])
        {
          auto candidate = first + FromWord(2 * i);
          if (test(candidate))
            return candidate;
        }
      first = first + stride;
    }
  }
}

#endif
