#include "biginteger/ops/NumberTheory.h"

#include <array>
#include <bit>
#include <stdexcept>

#include "biginteger/common/BitShifts.h"
#include "biginteger/common/LimbKernels.h"
#include "biginteger/ops/Operations.h"

namespace BigMath
{
  namespace
  {
    size_t BitLength(std::span<const DataT> bits)
    {
      while (bits.size() > 1 && bits.back() == 0)
        bits = bits.first(bits.size() - 1);
      return (bits.size() - 1) * LimbBits + std::bit_width(bits.back());
    }

    template <class Value, class Product>
    Value SlidingWindowPower(Value factor, std::span<const DataT> bits, size_t length, Product multiply)
    {
      const auto bit = [bits](size_t i) {
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
        for (size_t i = length - 1; i > 0; --i)
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
      for (size_t i = 1; i < (size_t{1} << (width - 1)); ++i)
      {
        powers[i] = powers[i - 1];
        multiply(powers[i], square);
      }

      const auto window = [bit, width](size_t &remaining) {
        size_t low = remaining > width ? remaining - width : 0;
        while (!bit(low))
          ++low;
        unsigned value = 0;
        while (remaining > low)
          value = (value << 1) | bit(--remaining);
        return value >> 1;
      };
      size_t remaining = length;
      auto result = powers[window(remaining)];
      while (remaining > 0)
      {
        if (!bit(remaining - 1))
        {
          multiply(result, result);
          --remaining;
          continue;
        }
        const size_t start = remaining;
        const auto index = window(remaining);
        for (size_t i = remaining; i < start; ++i)
          multiply(result, result);
        multiply(result, powers[index]);
      }
      return result;
    }

    DataT InvertOddLimb(DataT odd)
    {
      DataT inverse = 1;
      // Each Newton step doubles the number of correct low bits.
      for (unsigned bits = 1; bits < LimbBits; bits *= 2)
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
          for (size_t j = i + size; carry != 0; ++j)
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
            scratch(2 * static_cast<size_t>(size) + 1)
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

      BigInteger Power(BigInteger const &factor, BigInteger const &mod,
                       std::span<const DataT> bits, size_t length)
      {
        auto encoded = ((factor << size) % mod).GetInteger();
        encoded.resize(size);
        auto result = SlidingWindowPower(std::move(encoded), bits, length, [this](auto &a, auto &b) {
          Multiply(a, b);
        });
        std::copy(result.begin(), result.end(), scratch.begin());
        std::fill(scratch.begin() + size, scratch.end(), 0);
        Reduce(result);
        return BigInteger(std::move(result), false);
      }
    };

    BigInteger PowerReduced(BigInteger const &factor, BigInteger const &mod,
                            std::span<const DataT> bits, size_t length)
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
    return PowerReduced(factor, mod, bits, BitLength(bits));
  }
}
