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

  BigInteger Gcd(BigInteger const &a, BigInteger const &b)
  {
    auto left = BigInteger::FromLimbs(a.Limbs(), false);
    auto right = BigInteger::FromLimbs(b.Limbs(), false);
    while (!right.Zero())
    {
      if (right.size() == 1)
      {
        const DataT reduced = (left % right[0])[0];
        const DataT result = gcd(right[0], reduced);
        return BigInteger::FromLimbs({&result, 1}, false);
      }
      auto remainder = left % right;
      left = std::move(right);
      right = std::move(remainder);
    }
    return left;
  }

  ExtendedGcdResult ExtendedGcd(BigInteger const &a, BigInteger const &b)
  {
    if (a.Zero() && b.Zero())
      return {BigInteger(), BigInteger(), BigInteger()};

    auto left = BigInteger::FromLimbs(a.Limbs(), false);
    auto right = BigInteger::FromLimbs(b.Limbs(), false);
    constexpr DataT unit = 1;
    auto x = BigInteger::FromLimbs({&unit, 1}, a.IsNegative());
    auto nextY = BigInteger::FromLimbs({&unit, 1}, b.IsNegative());
    BigInteger y, nextX;
    // The current and next coefficient pairs represent left and right as
    // linear combinations of the original signed inputs.
    while (!right.Zero())
    {
      auto [quotient, remainder] = DivideAndRemainder(left, right);
      left = std::move(right);
      right = std::move(remainder);
      if (right.Zero())
      {
        // The next coefficients already represent the last nonzero remainder.
        return {std::move(left), std::move(nextX), std::move(nextY)};
      }
      auto updatedX = x - quotient * nextX;
      auto updatedY = y - quotient * nextY;
      x = std::move(nextX);
      y = std::move(nextY);
      nextX = std::move(updatedX);
      nextY = std::move(updatedY);
    }
    return {std::move(left), std::move(x), std::move(y)};
  }
}
