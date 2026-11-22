#include "biginteger/ops/NumberTheory.h"

#include <bit>
#include <cmath>
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

  namespace
  {
    template <class Word>
    Word SqrtWord(Word value)
    {
      // Even when the floating-point estimate rounds up to 2^(word bits/2),
      // clamping makes the trial square fit in Word.
      constexpr Word maxRoot = (Word{1} << (sizeof(Word) * 4)) - 1;
      Word root = std::min(static_cast<Word>(std::sqrt(static_cast<double>(value))), maxRoot);
      const Word square = root * root;
      if (square > value || value - square > root + root)
      {
        // If the estimate is low, the quotient provides an upper bound.
        root = std::max(root, value / root);
        for (Word quotient = value / root; quotient < root; quotient = value / root)
          root = (root + quotient) / 2;
      }
      return root;
    }

    template <class Word>
    std::pair<BigInteger, BigInteger> SqrtRemWord(Word value)
    {
      const Word root = SqrtWord(value);
      const Word remainder = value - root * root;
      const DataT rootLimb = static_cast<DataT>(root);
      const DataT remainderLimbs[] = {
          static_cast<DataT>(remainder) & LimbMask,
          static_cast<DataT>(static_cast<ULong128>(remainder) >> LimbBits)};
      return {BigInteger::FromLimbs({&rootLimb, 1}, false),
              BigInteger::FromLimbs(remainderLimbs, false)};
    }

    std::pair<BigInteger, BigInteger> SqrtRemMagnitude(std::span<const DataT> limbs)
    {
      if (limbs.size() == 1)
        return SqrtRemWord(limbs.front());
      if (limbs.size() == 2)
      {
        const auto value = (static_cast<ULong128>(limbs.back()) << LimbBits) | limbs.front();
        if constexpr (LimbBits == 32)
          return SqrtRemWord(static_cast<DataT>(value));
        else
          return SqrtRemWord(value);
      }

      const auto bits = (limbs.size() - 1) * uint64_t{LimbBits} + std::bit_width(limbs.back());
      // Choose a limb-aligned radix B so the high part has root at least B/2.
      const auto block = static_cast<SizeT>((bits + 1) / (4 * LimbBits));
      if (block == 0)
      {
        const auto value = BigInteger::FromLimbs(limbs, false);
        // Keep 63 or 64 leading bits by discarding an even number of bits.
        // The next integer above the prefix root gives an upper bound once
        // scaled back to the original width.
        const auto shift = (bits - 63) / 2;
        const auto first = static_cast<SizeT>(2 * shift / LimbBits);
        ULong128 prefix = 0;
        for (size_t i = limbs.size(); i > first; --i)
          prefix = (prefix << LimbBits) | limbs[i - 1];
        const auto top = static_cast<DataT>(prefix >> (2 * shift % LimbBits));
        const ULong128 seed = static_cast<ULong128>(SqrtWord(top) + 1) << (shift % LimbBits);
        const DataT seedLimbs[] = {static_cast<DataT>(seed) & LimbMask,
                                   static_cast<DataT>(seed >> LimbBits)};
        auto root = BigInteger::FromLimbs(seedLimbs, false) << static_cast<SizeT>(shift / LimbBits);
        for (;;)
        {
          const auto quotient = value / root;
          if (quotient >= root)
          {
            auto remainder = value - root * root;
            return {std::move(root), std::move(remainder)};
          }
          root = (root + quotient) / DataT{2};
        }
      }

      // Zimmermann's divide-and-conquer square root. Take the root of the
      // high part and use one division to recover the low part of the root.
      auto [highRoot, highRemainder] = SqrtRemMagnitude(limbs.subspan(2 * block));
      const auto middle = BigInteger::FromLimbs(limbs.subspan(block, block), false);
      const auto low = BigInteger::FromLimbs(limbs.first(block), false);
      auto [q, r] = DivideAndRemainder((highRemainder << block) + middle, highRoot + highRoot);
      auto root = (highRoot << block) + q;
      auto remainder = (r << block) + low - q * q;
      if (remainder.IsNegative())
      {
        // The lower bound on highRoot guarantees one correction suffices.
        constexpr DataT unit = 1;
        const auto one = BigInteger::FromLimbs({&unit, 1}, false);
        remainder = remainder + (root + root - one);
        root = root - one;
      }
      return {std::move(root), std::move(remainder)};
    }
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
    return SqrtRemMagnitude(limbs);
  }
}
