#include "biginteger/ops/NumberTheory.h"

#include <bit>
#include <stdexcept>

#include "biginteger/common/BitShifts.h"
#include "biginteger/common/Builder.h"
#include "biginteger/ops/Operations.h"

namespace BigMath
{
  namespace
  {
    // Both the base and exponent are positive. Zero signals that the
    // power exceeded the bound.
    BigInteger PowerAtMost(BigInteger const &base, UInt exponent, BigInteger const &limit)
    {
      const auto limitBits = LimbDetail::BitLength(limit.Limbs());
      const auto multiply = [&](BigInteger const &a, BigInteger const &b) {
        if (LimbDetail::BitLength(a.Limbs()) + LimbDetail::BitLength(b.Limbs()) - 1 > limitBits)
          return BigInteger();
        auto product = a * b;
        if (product > limit)
          return BigInteger();
        return product;
      };

      if (base > limit)
        return BigInteger();
      auto power = base;
      for (UInt bit = std::bit_floor(exponent) >> 1; bit != 0; bit >>= 1)
      {
        power = multiply(power, power);
        if (!power.Zero() && (exponent & bit))
          power = multiply(power, base);
        if (power.Zero())
          break;
      }
      return power;
    }

    BigInteger DropLowBits(BigInteger const &n, std::size_t bits)
    {
      auto top = n.Limbs().subspan(bits / LimbBits);
      if (bits % LimbBits == 0)
        return BigInteger::FromLimbs(top, false);
      return BigInteger(ShiftRightBits(top, bits % LimbBits, LimbBits), false);
    }

    ULong RootPrefix(BigInteger const &limit, UInt degree, unsigned prefixBits)
    {
      ULong low = ULong{1} << (prefixBits - 1);
      ULong high = 2 * low;
      while (high - low > 1)
      {
        const auto mid = low + (high - low) / 2;
        if (PowerAtMost(BigIntegerBuilder::From(mid), degree, limit).Zero())
          high = mid;
        else
          low = mid;
      }
      return low;
    }

    std::pair<BigInteger, BigInteger> RefineRoot(BigInteger const &n, UInt degree, BigInteger root)
    {
      for (;;)
      {
        auto denominator = PowerAtMost(root, degree - 1, n);
        auto next = root * static_cast<DataT>(degree - 1);
        if (!denominator.Zero())
          next = next + n / denominator;
        next = next / static_cast<DataT>(degree);
        // Integer Newton decreases from above until the first step that
        // would stay at the root or increase.
        if (next >= root)
          return {std::move(root), std::move(denominator)};
        root = std::move(next);
      }
    }

    std::pair<BigInteger, BigInteger> RootRemWord(ULong n, UInt degree, std::size_t rootBits)
    {
      const auto boundedPower = [n, degree](ULong base) {
        ULong128 power = base;
        for (UInt bit = std::bit_floor(degree) >> 1; bit != 0; bit >>= 1)
        {
          power *= power;
          if (power > n)
            return ULong{0};
          if (degree & bit)
          {
            power *= base;
            if (power > n)
              return ULong{0};
          }
        }
        return static_cast<ULong>(power);
      };

      ULong low = ULong{1} << (rootBits - 1);
      ULong high = 2 * low;
      ULong power = boundedPower(low);
      while (high - low > 1)
      {
        const auto mid = low + (high - low) / 2;
        if (const auto candidate = boundedPower(mid))
        {
          low = mid;
          power = candidate;
        }
        else
          high = mid;
      }
      return {BigIntegerBuilder::From(low), BigIntegerBuilder::From(n - power)};
    }

    std::pair<BigInteger, BigInteger> RootRemMagnitude(BigInteger const &n, UInt degree)
    {
      const auto bits = LimbDetail::BitLength(n.Limbs());
      if (degree >= bits)
      {
        auto one = BigIntegerBuilder::From(ULong{1});
        auto remainder = n - one;
        return {std::move(one), std::move(remainder)};
      }

      const auto rootBits = (bits - 1) / degree + 1;
      if (bits <= 64)
      {
        ULong word = n[0];
        if constexpr (LimbBits == 32)
          if (bits > 32)
            word |= n[1] << 32;
        return RootRemWord(word, degree, rootBits);
      }

      // Enough leading root bits to make the upper estimate's relative
      // error at most 1/degree, even for large degrees.
      const auto guard = std::bit_width(degree);
      const auto prefixBits = std::min<std::size_t>(rootBits, std::max<unsigned>(8, guard + 1));
      const auto droppedBits = rootBits - prefixBits;
      const auto leading = droppedBits == 0 ? BigInteger() : DropLowBits(n, degree * droppedBits);
      auto root = BigIntegerBuilder::From(RootPrefix(droppedBits == 0 ? n : leading, degree, prefixBits));
      if (droppedBits == 0)
      {
        auto remainder = n - PowerAtMost(root, degree, n);
        return {std::move(root), std::move(remainder)};
      }

      const auto one = BigIntegerBuilder::From(ULong{1});
      auto precision = prefixBits;
      for (;;)
      {
        // floor(root(n / 2^(degree*s))) is the full integer root with
        // its low s bits removed. Each exact prefix therefore bounds the
        // next precision from above after adding one and rescaling.
        const auto nextPrecision = std::min(rootBits, 2 * precision - guard);
        const auto shift = nextPrecision - precision;
        const auto upper = root + one;
        auto limbs = ShiftLeftBits(upper.Limbs(), shift % LimbBits, LimbBits);
        limbs.insert(limbs.begin(), shift / LimbBits, 0);
        BigInteger seed(std::move(limbs), false);
        if (nextPrecision == rootBits)
        {
          auto [finalRoot, denominator] = RefineRoot(n, degree, std::move(seed));
          auto remainder = n - finalRoot * denominator;
          return {std::move(finalRoot), std::move(remainder)};
        }
        auto truncated = DropLowBits(n, degree * (rootBits - nextPrecision));
        root = RefineRoot(truncated, degree, std::move(seed)).first;
        precision = nextPrecision;
      }
    }
  }

  std::pair<BigInteger, BigInteger> RootRem(BigInteger const &n, UInt degree)
  {
    if (degree == 0 || (n.IsNegative() && degree % 2 == 0))
      throw std::invalid_argument("RootRem requires a positive degree and allows negative inputs for odd degrees");
    if (n.Zero())
      return {BigInteger(), BigInteger()};
    if (degree == 1)
      return {n, BigInteger()};
    if (degree == 2)
      return SqrtRem(n);
    if (!n.IsNegative())
      return RootRemMagnitude(n, degree);

    auto [root, remainder] = RootRemMagnitude(BigInteger::FromLimbs(n.Limbs(), false), degree);
    root.SetSign(true);
    remainder.SetSign(true);
    return {std::move(root), std::move(remainder)};
  }
}
