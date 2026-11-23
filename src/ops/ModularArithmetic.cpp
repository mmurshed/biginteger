#include "biginteger/ops/NumberTheory.h"
#include "PrimeSearch.h"

#include <array>
#include <atomic>
#include <bit>
#include <exception>
#include <random>
#include <stdexcept>

#include "biginteger/common/BitShifts.h"
#include "biginteger/common/LimbKernels.h"
#include "biginteger/common/Parallel.h"
#include "biginteger/ops/Operations.h"

namespace BigMath
{
  namespace
  {
    template <class Value, class Product>
    Value SlidingWindowPower(Value factor, std::span<const DataT> bits, std::size_t length, Product multiply)
    {
      const auto bit = [bits](std::size_t i) {
        return (bits[i / LimbBits] >> (i % LimbBits)) & 1;
      };
      unsigned width = 1;
      if (length >= 32)
      {
        const unsigned candidate = length <= 256 ? 3 : (length <= 2048 ? 4 : 5);
        const unsigned minimumSetBits = 1u << candidate;
        unsigned setBits = 0;
        for (DataT limb : bits)
        {
          setBits += std::popcount(limb);
          if (setBits >= minimumSetBits)
          {
            width = candidate;
            break;
          }
        }
      }

      if (width == 1)
      {
        auto result = factor;
        for (std::size_t i = length - 1; i > 0; --i)
        {
          multiply(result, result);
          if (bit(i - 1))
            multiply(result, factor);
        }
        return result;
      }

      // Only odd powers appear at the end of a window.
      std::array<Value, 16> powers;
      powers[0] = std::move(factor);
      auto square = powers[0];
      multiply(square, square);
      for (std::size_t i = 1; i < (std::size_t{1} << (width - 1)); ++i)
      {
        powers[i] = powers[i - 1];
        multiply(powers[i], square);
      }

      const auto window = [bit, width](std::size_t &remaining) {
        std::size_t low = remaining > width ? remaining - width : 0;
        while (!bit(low))
          ++low;
        unsigned value = 0;
        while (remaining > low)
          value = (value << 1) | bit(--remaining);
        return value >> 1;
      };
      std::size_t remaining = length;
      auto result = powers[window(remaining)];
      while (remaining > 0)
      {
        if (!bit(remaining - 1))
        {
          multiply(result, result);
          --remaining;
          continue;
        }
        const std::size_t start = remaining;
        const auto index = window(remaining);
        for (std::size_t i = remaining; i < start; ++i)
          multiply(result, result);
        multiply(result, powers[index]);
      }
      return result;
    }

    DataT InvertOddLimb(DataT odd, unsigned width = LimbBits)
    {
      DataT inverse = 1;
      // Each Newton step doubles the number of correct low bits.
      for (unsigned bits = 1; bits < width; bits *= 2)
        inverse *= 2 - odd * inverse;
      return inverse;
    }

    class Montgomery
    {
      std::span<const DataT> modulus;
      SizeT size;
      DataT inverse;
      std::vector<DataT> scratch;

      void Reduce(std::vector<DataT> &result)
      {
        // Cancel each low limb successively. Dividing by the power of the
        // radix then leaves a residue below twice the modulus.
        for (SizeT i = 0; i < size; ++i)
        {
          const DataT digit = (scratch[i] * inverse) & LimbMask;
          DataT carry = 0;
          if constexpr (LimbBits == 64)
            carry = LimbKernel::AddMul1(scratch.data() + i, modulus.data(), size, digit);
          else
          {
            for (SizeT j = 0; j < size; ++j)
            {
              const DataT sum = scratch[i + j] + modulus[j] * digit + carry;
              scratch[i + j] = sum & LimbMask;
              carry = sum >> 32;
            }
          }
          for (std::size_t j = i + size; carry != 0; ++j)
          {
            const auto sum = static_cast<ULong128>(scratch[j]) + carry;
            scratch[j] = static_cast<DataT>(sum) & LimbMask;
            carry = static_cast<DataT>(sum >> LimbBits);
          }
        }

        auto high = std::span<DataT>(scratch).subspan(size, size);
        if (scratch.back() != 0 || Compare(high, modulus) >= 0)
        {
          if constexpr (LimbBits == 64)
            LimbKernel::SubInPlace(high.data(), size, modulus.data(), size);
          else
          {
            DataT borrow = 0;
            for (SizeT i = 0; i < size; ++i)
            {
              const DataT diff = LimbBase + high[i] - modulus[i] - borrow;
              high[i] = diff & LimbMask;
              borrow = 1 - (diff >> 32);
            }
          }
        }
        std::copy(high.begin(), high.end(), result.begin());
      }

    public:
      explicit Montgomery(BigInteger const &mod)
          : modulus(mod.Limbs()), size(mod.size()), inverse(0 - InvertOddLimb(modulus.front())),
            scratch(2 * static_cast<std::size_t>(size) + 1)
      {
      }

      void Multiply(std::vector<DataT> &a, std::vector<DataT> const &b)
      {
        if constexpr (LimbBits == 64)
        {
          if (&a == &b)
            LimbKernel::SqrBasecase(scratch.data(), a.data(), size);
          else
            LimbKernel::MulBasecase(scratch.data(), a.data(), size, b.data(), size);
        }
        else
        {
          std::fill(scratch.begin(), scratch.end(), 0);
          for (SizeT i = 0; i < size; ++i)
          {
            DataT carry = 0;
            for (SizeT j = 0; j < size; ++j)
            {
              const DataT sum = scratch[i + j] + a[i] * b[j] + carry;
              scratch[i + j] = sum & LimbMask;
              carry = sum >> 32;
            }
            scratch[i + size] = carry;
          }
        }
        scratch.back() = 0;
        Reduce(a);
      }

      std::vector<DataT> PowerEncoded(std::vector<DataT> factor,
                                      std::span<const DataT> bits, std::size_t length)
      {
        return SlidingWindowPower(std::move(factor), bits, length, [this](auto &a, auto &b) {
          Multiply(a, b);
        });
      }

      BigInteger Power(BigInteger const &factor, BigInteger const &mod,
                       std::span<const DataT> bits, std::size_t length)
      {
        auto encoded = ((factor << size) % mod).GetInteger();
        encoded.resize(size);
        auto result = PowerEncoded(std::move(encoded), bits, length);
        std::copy(result.begin(), result.end(), scratch.begin());
        std::fill(scratch.begin() + size, scratch.end(), 0);
        Reduce(result);
        return BigInteger(std::move(result), false);
      }
    };

    BigInteger PowerReduced(BigInteger const &factor, BigInteger const &mod,
                            std::span<const DataT> bits, std::size_t length)
    {
      if (length == 1 || (factor.size() == 1 && factor[0] <= 1))
        return factor;
      if (mod.size() == 1)
      {
        const DataT modulus = mod[0];
        DataT result;
        if ((modulus >> 32) == 0)
          result = SlidingWindowPower(factor[0], bits, length, [modulus](DataT &a, DataT b) {
            a = (a * b) % modulus;
          });
        else
          result = SlidingWindowPower(factor[0], bits, length, [modulus](DataT &a, DataT b) {
            a = static_cast<DataT>((static_cast<ULong128>(a) * b) % modulus);
          });
        return BigInteger::FromLimbs({&result, 1}, false);
      }

      // At large sizes retain the library's asymptotically faster multiplication
      // and division dispatch. At short exponents avoid the Montgomery setup.
      if ((mod[0] & 1) && mod.size() <= 128 && length >= 16)
        return Montgomery(mod).Power(factor, mod, bits, length);

      const unsigned zeros = std::countr_zero(mod[0]);
      if (zeros > 0 && zeros < LimbBits && mod.size() <= 128 && length >= 16)
      {
        // Reconstruct the answer from its residues modulo odd and 2^zeros.
        const DataT mask = (DataT{1} << zeros) - 1;
        const auto odd = BigInteger(ShiftRightBits(mod.Limbs(), zeros, LimbBits), false);
        const auto highResidue = PowerReduced(factor % odd, odd, bits, length);
        DataT lowResidue = factor[0] & mask;
        if (lowResidue > 1)
          lowResidue = SlidingWindowPower(lowResidue, bits, length, [mask](DataT &a, DataT b) {
            a = (a * b) & mask;
          });
        const DataT coefficient = ((lowResidue - highResidue[0]) * InvertOddLimb(odd[0])) & mask;
        return highResidue + odd * coefficient;
      }

      return SlidingWindowPower(factor, bits, length, [&mod](BigInteger &a, BigInteger const &b) {
        a = (a * b) % mod;
      });
    }

    template <class Value, class Square>
    bool PassesMillerRabin(Value value, Value const &one, Value const &minusOne,
                          std::size_t twos, Square square)
    {
      if (value == one || value == minusOne)
        return true;
      for (std::size_t i = 1; i < twos; ++i)
      {
        square(value);
        if (value == minusOne)
          return true;
        if (value == one)
          return false;
      }
      return false;
    }

    bool PrimeWord(DataT n)
    {
      for (const DataT p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47})
        if (n % p == 0)
          return n == p;
      if (n < 53 * 53)
        return true;

      const unsigned twos = std::countr_zero(n - 1);
      const DataT exponent = (n - 1) >> twos;
      const bool fits32 = n <= UINT32_MAX;
      const DataT inverse = fits32 ? 0 : 0 - InvertOddLimb(n, 64);
      const DataT one = fits32 ? 1 : static_cast<DataT>((ULong128{1} << 64) % n);
      const DataT minusOne = n - one;
      const auto product = [n, fits32, inverse](DataT a, DataT b) {
        if (fits32)
          return (a * b) % n;
        // In radix 2^64 the low limbs cancel. Their sum carries exactly
        // when the low limb of a*b is nonzero.
        const ULong128 ab = static_cast<ULong128>(a) * b;
        const auto low = static_cast<DataT>(ab);
        const ULong128 multiple = static_cast<ULong128>(low * inverse) * n;
        const ULong128 reduced = (ab >> 64) + (multiple >> 64) + (low != 0);
        return static_cast<DataT>(reduced >= n ? reduced - n : reduced);
      };
      const auto trial = [&](DataT base) {
        base %= n;
        if (base == 0)
          return true;
        if (!fits32)
          base = static_cast<DataT>((static_cast<ULong128>(base) << 64) % n);
        DataT power = base;
        for (unsigned bit = std::bit_width(exponent) - 1; bit > 0; --bit)
        {
          power = product(power, power);
          if ((exponent >> (bit - 1)) & 1)
            power = product(power, base);
        }
        return PassesMillerRabin(power, one, minusOne, twos, [&](DataT &value) {
          value = product(value, value);
        });
      };
      // Jaeschke's three bases for 32 bits and Sinclair's seven for 64 bits.
      if (fits32)
        return trial(2) && trial(7) && trial(61);
      for (const DataT base : {2, 325, 9375, 28178, 450775, 9780504, 1795265022})
        if (!trial(base))
          return false;
      return true;
    }

    template <class WitnessTestFactory>
    bool CheckRandomWitnesses(BigInteger const &n, BigInteger const &two,
                              int rounds, WitnessTestFactory makeTest, bool allowParallel = false)
    {
      const auto upper = n - two;
      const unsigned topBits = std::bit_width(upper.Limbs().back());
      const DataT mask = LimbMask >> (LimbBits - topBits);
      std::vector<DataT> randomLimbs(upper.size());
      std::random_device source;
      std::uniform_int_distribution<DataT> distribution(0, LimbMask);
      const auto sample = [&] {
        BigInteger witness;
        do
        {
          for (auto &limb : randomLimbs)
            limb = distribution(source);
          randomLimbs.back() &= mask;
          witness = BigInteger::FromLimbs(randomLimbs, false);
        } while (witness < two || witness > upper);
        return witness;
      };

      // Keep cheap screens and small tests on the calling thread. The
      // randomized rounds are independent once the base-2 screen passes.
      const auto workers = allowParallel && rounds >= 8 ? ParallelNumThreads() : 1;
      if (workers <= 1)
      {
        auto test = makeTest();
        for (int i = 0; i < rounds; ++i)
          if (!test(sample()))
            return false;
        return true;
      }

      constexpr int maxBatch = 32;
      std::vector<BigInteger> witnesses;
      witnesses.reserve(std::min(rounds, maxBatch));
      while (rounds > 0)
      {
        const int count = std::min(rounds, maxBatch);
        witnesses.clear();
        for (int i = 0; i < count; ++i)
          witnesses.emplace_back(sample());
        std::atomic<bool> stop{false};
        std::array<std::exception_ptr, maxBatch> failures;
        ParallelDo(count, [&](Int begin, Int end) {
          try
          {
            auto workerTest = makeTest();
            for (Int i = begin; i < end && !stop.load(std::memory_order_relaxed); ++i)
              if (!workerTest(witnesses[i]))
              {
                stop.store(true, std::memory_order_relaxed);
                break;
              }
          }
          catch (...)
          {
            // Exceptions must not escape the worker before it reports completion.
            failures[begin] = std::current_exception();
            stop.store(true, std::memory_order_relaxed);
          }
        });
        for (auto const &failure : failures)
          if (failure)
            std::rethrow_exception(failure);
        if (stop.load(std::memory_order_relaxed))
          return false;
        rounds -= count;
      }
      return true;
    }
  }

  BigInteger PowMod(BigInteger const &base, BigInteger const &exponent,
                    BigInteger const &modulus)
  {
    if (modulus.IsNegative() || modulus.Zero() || exponent.IsNegative())
      throw std::invalid_argument("PowMod requires a positive modulus and a nonnegative exponent");

    constexpr DataT unit = 1;
    const auto one = BigInteger::FromLimbs({&unit, 1}, false);
    const auto mod = BigInteger::FromLimbs(modulus.Limbs(), false);
    if (mod == one)
      return BigInteger();
    if (exponent.Zero())
      return one;

    auto factor = BigInteger::FromLimbs(base.Limbs(), base.IsNegative()) % mod;
    if (factor.IsNegative())
      factor = factor + mod;
    const auto bits = exponent.Limbs();
    return PowerReduced(factor, mod, bits, LimbDetail::BitLength(bits));
  }

  bool IsProbablePrime(BigInteger const &n, int rounds)
  {
    if (rounds <= 0)
      throw std::invalid_argument("IsProbablePrime requires a positive round count");
    if (n.IsNegative() || n.Zero())
      return false;
    const auto bits = LimbDetail::BitLength(n.Limbs());
    if (bits <= 64)
    {
      DataT value = n[0];
      if constexpr (LimbBits == 32)
        if (n.size() > 1)
          value |= n[1] << 32;
      return value >= 2 && PrimeWord(value);
    }
    if (!(n[0] & 1))
      return false;

    // One reduction per group tests all of its prime factors.
    constexpr DataT group1 = DataT{3} * 5 * 7 * 11 * 13 * 17 * 19 * 23 * 29;
    constexpr DataT group2 = DataT{31} * 37 * 41 * 43 * 47;
    if (gcd((n % group1)[0], group1) != 1 || gcd((n % group2)[0], group2) != 1)
      return false;

    return PrimeSearch::TestLargeCandidate(
        BigInteger::FromLimbs(n.Limbs(), false), rounds, PrimeSearch::WitnessMode::All);
  }

  bool PrimeSearch::TestLargeCandidate(BigInteger const &mod, int rounds, WitnessMode mode)
  {
    constexpr DataT unit = 1, pair = 2;
    const auto one = BigInteger::FromLimbs({&unit, 1}, false);
    const auto two = BigInteger::FromLimbs({&pair, 1}, false);
    const auto minusOne = mod - one;
    auto oddLimbs = minusOne.Limbs();
    while (oddLimbs.front() == 0)
      oddLimbs = oddLimbs.subspan(1);
    const unsigned trailing = std::countr_zero(oddLimbs.front());
    const std::size_t twos = (minusOne.size() - oddLimbs.size()) * LimbBits + trailing;
    const auto exponent = ShiftRightBits(oddLimbs, trailing, LimbBits);
    const auto length = LimbDetail::BitLength(exponent);

    if (mod.size() <= 128)
    {
      const auto montOne = (one << mod.size()) % mod;
      auto encodedOne = montOne.GetInteger();
      auto encodedMinusOne = (mod - montOne).GetInteger();
      encodedOne.resize(mod.size());
      encodedMinusOne.resize(mod.size());
      const auto testEncoded = [&](Montgomery &context, std::vector<DataT> encoded) {
        auto power = context.PowerEncoded(std::move(encoded), exponent, length);
        return PassesMillerRabin(std::move(power), encodedOne, encodedMinusOne, twos, [&](auto &value) {
          context.Multiply(value, value);
        });
      };

      if (mode != WitnessMode::RandomOnly)
      {
        auto montTwo = montOne + montOne;
        if (montTwo >= mod)
          montTwo = montTwo - mod;
        auto encodedTwo = montTwo.ReleaseInteger();
        encodedTwo.resize(mod.size());
        Montgomery context(mod);
        const bool passes = testEncoded(context, std::move(encodedTwo));
        if (!passes || mode == WitnessMode::BaseTwo)
          return passes;
      }

      auto radixSquared = ((montOne << mod.size()) % mod).GetInteger();
      radixSquared.resize(mod.size());
      const auto makeTest = [&] {
        return [&, context = Montgomery(mod)](BigInteger const &witness) mutable {
          auto encoded = witness.GetInteger();
          encoded.resize(mod.size());
          context.Multiply(encoded, radixSquared);
          return testEncoded(context, std::move(encoded));
        };
      };
      return CheckRandomWitnesses(mod, two, rounds, makeTest, LimbDetail::BitLength(mod.Limbs()) >= 1024);
    }

    const auto test = [&](BigInteger const &witness) {
      auto power = PowerReduced(witness, mod, exponent, length);
      return PassesMillerRabin(std::move(power), one, minusOne, twos, [&](BigInteger &value) {
        value = (value * value) % mod;
      });
    };
    if (mode != WitnessMode::RandomOnly)
    {
      const bool passes = test(two);
      if (!passes || mode == WitnessMode::BaseTwo)
        return passes;
    }
    return CheckRandomWitnesses(mod, two, rounds, [&] { return test; });
  }
}
