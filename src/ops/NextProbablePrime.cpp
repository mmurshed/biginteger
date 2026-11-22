#include "PrimeSearch.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <exception>
#include <stdexcept>

#include "biginteger/common/Parallel.h"
#include "biginteger/ops/Operations.h"

namespace BigMath::PrimeSearch
{
  namespace
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

    constexpr SizeT MaxBatch = 8;

    BigInteger TestBatch(std::span<const BigInteger> candidates, int rounds)
    {
      std::array<bool, MaxBatch> passes{};
      std::array<std::exception_ptr, MaxBatch> failures;
      ParallelDo(static_cast<Int>(candidates.size()), [&](Int begin, Int end) {
        try
        {
          for (Int i = begin; i < end; ++i)
            passes[i] = TestLargeCandidate(candidates[i], rounds, WitnessMode::BaseTwo);
        }
        catch (...)
        {
          // Let the pool register this chunk as finished before rethrowing.
          failures[begin] = std::current_exception();
        }
      });
      for (auto const &failure : failures)
        if (failure)
          std::rethrow_exception(failure);

      // Test in order even when later candidates finished screening first.
      for (size_t i = 0; i < candidates.size(); ++i)
        if (passes[i] && TestLargeCandidate(candidates[i], rounds, WitnessMode::RandomOnly))
          return candidates[i];
      return BigInteger();
    }
  }

  BigInteger FindNext(BigInteger const &n, int rounds, SearchStats *stats)
  {
    if (rounds <= 0)
      throw std::invalid_argument("NextProbablePrime requires a positive round count");
    const auto one = FromWord(1);
    const auto two = FromWord(2);
    const auto start = BigInteger::FromLimbs(n.Limbs(), n.IsNegative());
    if (start < two)
      return two;

    const auto record = [stats](size_t count) {
      if (stats)
        stats->screenedCandidates += count;
    };
    auto first = start + ((start[0] & 1) ? two : one);
    if (start.size() <= 64 / LimbBits)
    {
      // The word-size path also handles primes in the sieve table.
      for (;; first = first + two)
      {
        record(1);
        if (IsProbablePrime(first, rounds))
          return first;
      }
    }

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

    // Try the first survivor on its own for short gaps. Later tests can
    // overlap when the candidates are large enough to amortize a dispatch.
    const auto bits = (first.Limbs().size() - 1) * LimbBits + std::bit_width(first.Limbs().back());
    const bool allowParallel = bits >= 1024 && first.size() <= 128;
    bool firstSurvivor = true;
    SizeT batchSize = 1;
    std::vector<BigInteger> batch;
    const auto stride = FromWord(2 * SegmentSize);
    for (;;)
    {
      std::array<bool, SegmentSize> composite{};
      for (size_t i = 0; i < primeCount; ++i)
      {
        unsigned offset = offsets[i];
        for (; offset < SegmentSize; offset += SmallPrimes[i])
          composite[offset] = true;
        offsets[i] = static_cast<uint16_t>(offset - SegmentSize);
      }
      for (unsigned i = 0; i < SegmentSize; ++i)
      {
        if (composite[i])
          continue;
        auto candidate = first + FromWord(2 * i);
        if (batchSize == 1)
        {
          record(1);
          if (TestLargeCandidate(candidate, rounds, WitnessMode::All))
            return candidate;
          if (firstSurvivor && allowParallel)
          {
            batchSize = std::min(MaxBatch, ParallelNumThreads());
            batch.reserve(batchSize);
          }
          firstSurvivor = false;
        }
        else
        {
          batch.emplace_back(std::move(candidate));
          if (batch.size() == batchSize)
          {
            record(batch.size());
            auto result = TestBatch(batch, rounds);
            if (!result.Zero())
              return result;
            batch.clear();
          }
        }
      }
      first = first + stride;
    }
  }
}

namespace BigMath
{
  BigInteger NextProbablePrime(BigInteger const &n, int rounds)
  {
    return PrimeSearch::FindNext(n, rounds);
  }
}
