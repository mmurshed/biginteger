#include "biginteger/ops/NumberTheory.h"

#include <bit>
#include <stdexcept>

#include "biginteger/ops/Operations.h"

namespace BigMath
{
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
    if (factor.Zero() || factor == one)
      return factor;

    auto bits = exponent.Limbs();
    while (bits.back() == 0)
      bits = bits.first(bits.size() - 1);
    // The leading set bit accounts for the initial factor.
    BigInteger result = factor;
    unsigned width = std::bit_width(bits.back()) - 1;
    for (size_t i = bits.size(); i > 0; --i)
    {
      for (unsigned bit = width; bit > 0; --bit)
      {
        result = (result * result) % mod;
        if ((bits[i - 1] >> (bit - 1)) & 1)
          result = (result * factor) % mod;
      }
      width = LimbBits;
    }
    return result;
  }
}
