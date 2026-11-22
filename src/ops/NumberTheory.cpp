#include "biginteger/ops/NumberTheory.h"

#include <bit>
#include <stdexcept>

#include "biginteger/ops/Operations.h"

namespace BigMath
{
  namespace
  {
    constexpr int64_t LehmerCoefficientLimit = int64_t{1} << 30;

    struct LehmerRow
    {
      int64_t left;
      int64_t right;

      bool Fits() const
      {
        constexpr auto limit = LehmerCoefficientLimit;
        return -limit < left && left < limit && -limit < right && right < limit;
      }

      BigInteger Apply(BigInteger const &a, BigInteger const &b) const
      {
        const auto scaled = [](BigInteger const &n, int64_t coefficient) {
          auto product = n * static_cast<DataT>(coefficient < 0 ? -coefficient : coefficient);
          product.SetSign(n.IsNegative() != (coefficient < 0));
          return product;
        };
        return scaled(a, left) + scaled(b, right);
      }
    };

    struct LehmerMatrix
    {
      LehmerRow first{1, 0};
      LehmerRow second{0, 1};

      void Apply(BigInteger &a, BigInteger &b) const
      {
        auto nextA = first.Apply(a, b);
        b = second.Apply(a, b);
        a = std::move(nextA);
      }
    };

    // Simulate Euclidean steps on a common prefix of at most 62 bits.
    // Each omitted suffix is in [0, 1) in units of that prefix. The
    // current matrix row bounds its contribution to each remainder.
    bool BuildLehmerMatrix(BigInteger const &a, BigInteger const &b, LehmerMatrix &matrix)
    {
      const auto left = a.Limbs();
      const auto right = b.Limbs();
      if (left.size() != right.size() || left.size() < 2)
        return false;

      const unsigned width = LimbBits + std::bit_width(std::max(left.back(), right.back()));
      const unsigned shift = width > 62 ? width - 62 : 0;
      const auto prefix = [shift](std::span<const DataT> limbs) {
        const auto high = static_cast<ULong128>(limbs.back()) << LimbBits;
        return static_cast<int64_t>((high | limbs[limbs.size() - 2]) >> shift);
      };
      int64_t u = prefix(left);
      int64_t v = prefix(right);
      unsigned steps = 0;
      while (v > 0 && u >= v)
      {
        const auto [a0, a1] = matrix.first;
        const auto [b0, b1] = matrix.second;
        const int64_t lowerU = u + std::min<int64_t>(0, a0) + std::min<int64_t>(0, a1);
        const int64_t upperU = u + std::max<int64_t>(0, a0) + std::max<int64_t>(0, a1);
        const int64_t lowerV = v + std::min<int64_t>(0, b0) + std::min<int64_t>(0, b1);
        const int64_t upperV = v + std::max<int64_t>(0, b0) + std::max<int64_t>(0, b1);
        const int64_t q = u / v;
        // Accept q only if every possible quotient is in [q, q+1).
        // The prefix and coefficient bounds keep each product in int64_t.
        if (q >= LehmerCoefficientLimit || lowerV <= 0 ||
            lowerU < q * upperV || upperU >= (q + 1) * lowerV)
          break;
        const LehmerRow next{a0 - q * b0, a1 - q * b1};
        if (!next.Fits())
          break;
        matrix.first = matrix.second;
        matrix.second = next;
        const int64_t remainder = u - q * v;
        u = v;
        v = remainder;
        ++steps;
      }
      return steps >= 2;
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
      LehmerMatrix matrix;
      if (BuildLehmerMatrix(left, right, matrix))
      {
        matrix.Apply(left, right);
        continue;
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
      LehmerMatrix matrix;
      if (BuildLehmerMatrix(left, right, matrix))
      {
        matrix.Apply(left, right);
        matrix.Apply(x, nextX);
        matrix.Apply(y, nextY);
        continue;
      }
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

  std::pair<BigInteger, BigInteger> SqrtRem(BigInteger const &n)
  {
    if (n.IsNegative())
      throw std::invalid_argument("SqrtRem requires a nonnegative input");
    if (n.Zero())
      return {BigInteger(), BigInteger()};

    auto limbs = n.Limbs();
    while (limbs.back() == 0)
      limbs = limbs.first(limbs.size() - 1);
    const auto bits = (limbs.size() - 1) * uint64_t{LimbBits} + std::bit_width(limbs.back());
    // Start at 2^ceil(bits / 2), above the true root. Integer Newton steps
    // decrease this upper bound until division reports n/root >= root.
    const auto seedBit = (bits + 1) / 2;
    const DataT seedLimb = DataT{1} << (seedBit % LimbBits);
    if (limbs.size() == 1)
    {
      const DataT value = limbs.front();
      DataT root = seedLimb;
      for (DataT quotient = value / root; quotient < root; quotient = value / root)
        root = (root + quotient) / 2;
      const DataT remainder = value - root * root;
      return {BigInteger::FromLimbs({&root, 1}, false),
              BigInteger::FromLimbs({&remainder, 1}, false)};
    }

    auto root = BigInteger::FromLimbs({&seedLimb, 1}, false) << static_cast<SizeT>(seedBit / LimbBits);
    for (;;)
    {
      const auto quotient = n / root;
      if (quotient >= root)
      {
        auto remainder = n - root * root;
        return {std::move(root), std::move(remainder)};
      }
      root = (root + quotient) / DataT{2};
    }
  }
}
