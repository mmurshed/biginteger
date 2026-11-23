#include "biginteger/ops/NumberTheory.h"

#include <bit>
#include <cmath>
#include <stdexcept>

#include "biginteger/build/DispatchThresholds.h"
#include "biginteger/common/BitShifts.h"
#include "biginteger/ops/Operations.h"

namespace BigMath
{
  namespace
  {
    constexpr Long LehmerCoefficientLimit = Long{1} << 30;

    struct LehmerRow
    {
      Long left;
      Long right;

      bool Fits() const
      {
        constexpr auto limit = LehmerCoefficientLimit;
        return -limit < left && left < limit && -limit < right && right < limit;
      }

      BigInteger Apply(BigInteger const &a, BigInteger const &b) const
      {
        const auto scaled = [](BigInteger const &n, Long coefficient) {
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
        return static_cast<Long>((high | limbs[limbs.size() - 2]) >> shift);
      };
      Long u = prefix(left);
      Long v = prefix(right);
      unsigned steps = 0;
      while (v > 0 && u >= v)
      {
        const auto [a0, a1] = matrix.first;
        const auto [b0, b1] = matrix.second;
        const Long lowerU = u + std::min<Long>(0, a0) + std::min<Long>(0, a1);
        const Long upperU = u + std::max<Long>(0, a0) + std::max<Long>(0, a1);
        const Long lowerV = v + std::min<Long>(0, b0) + std::min<Long>(0, b1);
        const Long upperV = v + std::max<Long>(0, b0) + std::max<Long>(0, b1);
        const Long q = u / v;
        // Accept q only if every possible quotient is in [q, q+1).
        // The prefix and coefficient bounds keep each product in Long.
        if (q >= LehmerCoefficientLimit || lowerV <= 0 ||
            lowerU < q * upperV || upperU >= (q + 1) * lowerV)
          break;
        const LehmerRow next{a0 - q * b0, a1 - q * b1};
        if (!next.Fits())
          break;
        matrix.first = matrix.second;
        matrix.second = next;
        const Long remainder = u - q * v;
        u = v;
        v = remainder;
        ++steps;
      }
      return steps >= 2;
    }

    BigInteger FromMagnitudeWord(DataT value)
    {
      return BigInteger::FromLimbs({&value, 1}, false);
    }

    std::size_t GcdBits(BigInteger const &value)
    {
      return LimbDetail::BitLength(value.Limbs());
    }

    // Store the magnitudes of a determinant-one transform with positive
    // diagonal and nonpositive off-diagonal. The inverse is nonnegative.
    struct HalfGcdMatrix
    {
      BigInteger aa = FromMagnitudeWord(1);
      BigInteger ab;
      BigInteger ba;
      BigInteger bb = FromMagnitudeWord(1);

      bool IsIdentity() const
      {
        // With determinant one, zero off-diagonal entries imply unit diagonal.
        return ab.Zero() && ba.Zero();
      }

      static BigInteger Product(BigInteger const &value, BigInteger const &factor)
      {
        return factor.size() == 1 ? value * factor[0] : value * factor;
      }

      void Apply(BigInteger &a, BigInteger &b) const
      {
        auto first = Product(a, aa) - Product(b, ab);
        b = Product(b, bb) - Product(a, ba);
        a = std::move(first);
      }

      void Prepend(HalfGcdMatrix step)
      {
        if (step.IsIdentity())
          return;
        if (IsIdentity())
        {
          *this = std::move(step);
          return;
        }
        *this = {
            Product(aa, step.aa) + Product(ba, step.ab),
            Product(ab, step.aa) + Product(bb, step.ab),
            Product(aa, step.ba) + Product(ba, step.bb),
            Product(ab, step.ba) + Product(bb, step.bb)};
      }

      void SubtractMultiple(bool reduceFirst, BigInteger const &quotient)
      {
        if (reduceFirst)
        {
          aa = aa + Product(ba, quotient);
          ab = ab + Product(bb, quotient);
        }
        else
        {
          ba = ba + Product(aa, quotient);
          bb = bb + Product(ab, quotient);
        }
      }
    };

    bool TryBoundedLehmer(BigInteger &a, BigInteger &b, HalfGcdMatrix &matrix,
                          std::size_t stop)
    {
      const bool reversed = a < b;
      auto &larger = reversed ? b : a;
      auto &smaller = reversed ? a : b;
      LehmerMatrix batch;
      if (!BuildLehmerMatrix(larger, smaller, batch))
        return false;

      auto second = batch.second.Apply(larger, smaller);
      if (GcdBits(second) <= stop)
        return false;
      auto first = batch.first.Apply(larger, smaller);

      // Undo the implicit swap from an odd number of Euclidean steps.
      if (batch.first.left < 0)
      {
        std::swap(first, second);
        std::swap(batch.first, batch.second);
      }
      larger = std::move(first);
      smaller = std::move(second);

      // Express the batch in the caller's original coordinate order.
      if (reversed)
      {
        std::swap(batch.first, batch.second);
        std::swap(batch.first.left, batch.first.right);
        std::swap(batch.second.left, batch.second.right);
      }
      matrix.Prepend({
          FromMagnitudeWord(static_cast<DataT>(batch.first.left)),
          FromMagnitudeWord(static_cast<DataT>(-batch.first.right)),
          FromMagnitudeWord(static_cast<DataT>(-batch.second.left)),
          FromMagnitudeWord(static_cast<DataT>(batch.second.right))});
      return true;
    }

    // Reduce the larger operand without swapping the two coordinates. Both
    // operands must retain more than 'stop' bits so this transform is safe
    // to apply to longer operands with the same high bits.
    bool HalfGcdStep(BigInteger &a, BigInteger &b, HalfGcdMatrix &matrix,
                     std::size_t stop)
    {
      const bool reduceFirst = a >= b;
      auto &larger = reduceFirst ? a : b;
      auto &smaller = reduceFirst ? b : a;
      auto difference = larger - smaller;
      if (GcdBits(difference) <= stop)
        return false;
      if (difference <= smaller)
      {
        larger = std::move(difference);
        matrix.SubtractMultiple(reduceFirst, FromMagnitudeWord(1));
        return true;
      }

      auto [quotient, remainder] = DivideAndRemainder(larger, smaller);
      if (GcdBits(remainder) <= stop)
      {
        quotient = quotient - FromMagnitudeWord(1);
        remainder = remainder + smaller;
      }
      larger = std::move(remainder);
      matrix.SubtractMultiple(reduceFirst, quotient);
      return true;
    }

    BigInteger GcdHighBits(BigInteger const &value, std::size_t drop)
    {
      const auto high = value.Limbs().subspan(drop / LimbBits);
      const Int bits = drop % LimbBits;
      if (bits == 0)
        return BigInteger::FromLimbs(high, false);
      return BigInteger(ShiftRightBits(high, bits, LimbBits), false);
    }

    HalfGcdMatrix HalfGcd(BigInteger &a, BigInteger &b);

    HalfGcdMatrix HalfGcdFromHighBits(BigInteger &a, BigInteger &b, std::size_t drop)
    {
      auto highA = GcdHighBits(a, drop);
      auto highB = GcdHighBits(b, drop);
      auto matrix = HalfGcd(highA, highB);
      if (!matrix.IsIdentity())
        matrix.Apply(a, b);
      return matrix;
    }

    // Möller's difference-based half-GCD (section 6, DOI 10.1090/S0025-5718-07-02017-0).
    // Return T such that (a, b) = T * (original a, original b). When both
    // inputs exceed 'stop' bits, both outputs do too and the difference
    // has at most 'stop' bits.
    HalfGcdMatrix HalfGcd(BigInteger &a, BigInteger &b)
    {
      const auto bits = std::max(GcdBits(a), GcdBits(b));
      const auto stop = bits / 2 + 1;
      HalfGcdMatrix matrix;
      if (std::min(GcdBits(a), GcdBits(b)) <= stop)
        return matrix;
      if (bits <= BIGMATH_HALF_GCD_LEAF_THRESHOLD)
      {
        while (TryBoundedLehmer(a, b, matrix, stop) || HalfGcdStep(a, b, matrix, stop))
          ;
        return matrix;
      }

      const auto threeQuarters = bits - (bits + 3) / 4;
      if (std::min(GcdBits(a), GcdBits(b)) > threeQuarters + 2)
        matrix = HalfGcdFromHighBits(a, b, bits / 2);

      // Repair the few bits of slack introduced by applying the transform
      // to the low bits. The next recursion also starts on at most half
      // of the original bits.
      auto remainingBits = std::max(GcdBits(a), GcdBits(b));
      while (remainingBits > threeQuarters + 1)
      {
        if (!HalfGcdStep(a, b, matrix, stop))
          return matrix;
        remainingBits = std::max(GcdBits(a), GcdBits(b));
      }
      // Align the second recursive stop with this level's stop after lifting.
      if (std::min(GcdBits(a), GcdBits(b)) > stop + 2)
        matrix.Prepend(HalfGcdFromHighBits(a, b, 2 * stop - remainingBits + 1));

      while (HalfGcdStep(a, b, matrix, stop))
        ;
      return matrix;
    }

    bool UseHalfGcd(BigInteger const &a, BigInteger const &b, std::size_t threshold)
    {
      if (threshold == 0)
        return false;
      const auto minLimbs = std::min(a.size(), b.size());
      if (minLimbs < (threshold + LimbBits - 1) / LimbBits)
        return false;
      const auto aBits = GcdBits(a);
      const auto bBits = GcdBits(b);
      const auto smaller = std::min(aBits, bBits);
      return smaller >= threshold && smaller > std::max(aBits, bBits) / 2 + 1;
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
      if (UseHalfGcd(left, right, BIGMATH_HALF_GCD_THRESHOLD))
      {
        // A single high-half reduction avoids composing a full matrix
        // here. Plain GCD needs only the two reduced operands.
        HalfGcdFromHighBits(left, right, std::max(GcdBits(left), GcdBits(right)) / 2);
        if (left < right)
          std::swap(left, right);
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
      if (UseHalfGcd(left, right, BIGMATH_HALF_GCD_EXTENDED_THRESHOLD))
      {
        const auto matrix = HalfGcd(left, right);
        if (!matrix.IsIdentity())
        {
          matrix.Apply(x, nextX);
          matrix.Apply(y, nextY);
        }
        if (left < right)
        {
          std::swap(left, right);
          std::swap(x, nextX);
          std::swap(y, nextY);
        }
      }
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

      const auto bits = LimbDetail::BitLength(limbs);
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
        for (std::size_t i = limbs.size(); i > first; --i)
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
