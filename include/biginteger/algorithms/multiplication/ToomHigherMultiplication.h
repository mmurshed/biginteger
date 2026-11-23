// Higher-degree Toom products on little-endian 64-bit limb arrays.
#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

#include "../../common/Comparator.h"
#include "../../common/LimbKernels.h"
#include "../../common/Util.h"
#include "KaratsubaMultiplication.h"

namespace BigMath
{
  namespace ToomHigherDetail
  {
    using Buffer = std::span<DataT>;
    using Limbs = std::span<const DataT>;

    inline Limbs Piece(Limbs value, SizeT width, unsigned index)
    {
      const size_t offset = size_t{width} * index;
      return offset < value.size() ? value.subspan(offset, std::min<size_t>(width, value.size() - offset))
                                   : Limbs{};
    }

    inline SizeT NonzeroLength(Limbs value)
    {
      const SizeT n = LimbDetail::TrimmedLimbCount(value);
      return n == 1 && value.front() == 0 ? 0 : n;
    }

    inline void MultiplyPoint(Buffer out, Limbs a, Limbs b, DataT *scratch)
    {
      const SizeT la = NonzeroLength(a), lb = NonzeroLength(b);
      if (la && lb)
        KaratsubaMultiplication::MultiplyInto(a.data(), la, b.data(), lb, out.data(), scratch, Base2_64);
      std::fill(out.begin() + (la && lb ? la + lb : 0), out.end(), 0);
    }

    // Scale a dyadic evaluation at x = 2^power by the denominator to
    // the operand's degree. Even and odd terms are both nonnegative.
    inline void EvaluateParts(Limbs value, SizeT width, unsigned pieces, int power,
                              Buffer even, Buffer odd = {})
    {
      std::fill(even.begin(), even.end(), 0);
      std::fill(odd.begin(), odd.end(), 0);
      for (unsigned i = 0; i < pieces; ++i)
      {
        const auto chunk = Piece(value, width, i);
        if (chunk.empty())
          break;
        const unsigned shift = power >= 0 ? i * static_cast<unsigned>(power)
                                         : (pieces - 1 - i) * static_cast<unsigned>(-power);
        const auto dest = (i % 2 && !odd.empty()) ? odd : even;
        const auto shared = static_cast<SizeT>(chunk.size());
        const DataT carry = LimbKernel::AddMul1(dest.data(), chunk.data(), shared, DataT{1} << shift);
        const auto tail = dest.subspan(shared);
        [[maybe_unused]] const DataT overflow = LimbKernel::AddCarry(
            tail.data(), tail.data(), static_cast<SizeT>(tail.size()), carry);
        assert(overflow == 0);
      }
    }

    inline bool EvaluatePair(Limbs value, SizeT width, unsigned pieces, int power,
                             Buffer plus, Buffer minus, Buffer temporary)
    {
      EvaluateParts(value, width, pieces, power, plus, minus);
      const bool negative = Compare(Limbs(plus), Limbs(minus)) < 0;
      const auto n = static_cast<SizeT>(plus.size());
      [[maybe_unused]] const DataT carry = LimbKernel::AddN(temporary.data(), plus.data(), minus.data(), n);
      [[maybe_unused]] const DataT borrow = negative
          ? LimbKernel::SubN(minus.data(), minus.data(), plus.data(), n)
          : LimbKernel::SubN(minus.data(), plus.data(), minus.data(), n);
      assert(carry == 0 && borrow == 0);
      std::copy_n(temporary.data(), n, plus.data());
      return negative;
    }

    inline void Shift(Buffer value, int bits)
    {
      [[maybe_unused]] DataT discarded = 0;
      const auto n = static_cast<SizeT>(value.size());
      assert(bits > -64 && bits < 64);
      if (bits > 0)
        discarded = LimbKernel::Lshift(value.data(), value.data(), n, static_cast<unsigned>(bits));
      else if (bits < 0)
        discarded = LimbKernel::Rshift(value.data(), value.data(), n, static_cast<unsigned>(-bits));
      assert(discarded == 0);
    }

    inline void SubtractMultiple(Buffer value, Limbs term, DataT factor)
    {
      assert(value.size() == term.size());
      [[maybe_unused]] const DataT borrow = LimbKernel::SubMul1(
          value.data(), term.data(), static_cast<SizeT>(value.size()), factor);
      assert(borrow == 0);
    }

    struct ExactDivisor
    {
      int shift = 0;
      DataT odd = 1;
      DataT inverse = 1;

      constexpr ExactDivisor() = default;
      constexpr explicit ExactDivisor(DataT value)
          : shift(std::countr_zero(value)), odd(value >> shift), inverse(odd)
      {
        // Newton iteration doubles the correct low bits of the odd inverse.
        for (unsigned i = 0; i < 6; ++i)
          inverse *= 2 - odd * inverse;
      }
    };

    inline void DivideExact(Buffer value, ExactDivisor divisor)
    {
      Shift(value, -divisor.shift);
      if (divisor.odd == 1)
        return;
      DataT carry = 0;
      for (DataT &limb : value)
      {
        const DataT quotient = (limb - carry) * divisor.inverse;
        carry = static_cast<DataT>((ULong128{quotient} * divisor.odd + carry) >> 64);
        limb = quotient;
      }
      assert(carry == 0);
    }

    template <size_t Count>
    constexpr auto InterpolationDivisors()
    {
      std::array<std::array<ExactDivisor, Count>, Count> divisors{};
      for (unsigned order = 1; order < Count; ++order)
        for (unsigned i = order; i < Count; ++i)
          divisors[order][i] = ExactDivisor((DataT{1} << (2 * i)) - (DataT{1} << (2 * (i - order))));
      return divisors;
    }

    template <size_t Count>
    void Interpolate(std::array<Buffer, Count> &samples, unsigned gridShift)
    {
      // Nodes on the scaled grid are 1, 4, 16, ... . First form Newton
      // divided differences, then expand the Newton factors in place.
      // All divisions are exact for an integer-coefficient polynomial.
      // Divided differences of its nonnegative monomials and their
      // coefficients during the expansion are also nonnegative.
      constexpr auto divisors = InterpolationDivisors<Count>();
      for (size_t order = 1; order < Count; ++order)
        for (size_t i = Count; i-- > order;)
        {
          const auto n = static_cast<SizeT>(samples[i].size());
          [[maybe_unused]] const DataT borrow = LimbKernel::SubN(
              samples[i].data(), samples[i].data(), samples[i - 1].data(), n);
          assert(borrow == 0);
          DivideExact(samples[i], divisors[order][i]);
        }
      for (size_t i = Count - 1; i-- > 0;)
        for (size_t j = i; j + 1 < Count; ++j)
          SubtractMultiple(samples[j], samples[j + 1], DataT{1} << (2 * i));
      for (size_t i = 0; i < Count; ++i)
        Shift(samples[i], -static_cast<int>(gridShift * (Count - 1 - i)));
    }

    template <unsigned LeftPieces, unsigned RightPieces>
    std::vector<DataT> MultiplySplit(Limbs a, Limbs b)
    {
      constexpr unsigned degree = LeftPieces + RightPieces - 2;
      constexpr unsigned evenCount = (degree - 1) / 2;
      constexpr unsigned oddCount = degree / 2;
      constexpr unsigned gridShift = evenCount - 1;
      const SizeT k = static_cast<SizeT>(std::max((a.size() + LeftPieces - 1) / LeftPieces,
                                                  (b.size() + RightPieces - 1) / RightPieces));
      const SizeT inputWidth = k + 1;
      // At the largest dyadic point each evaluation has fewer than 32
      // extra bits. After scaling for interpolation even the largest
      // sample is below 2^(128k+80), so two extra limbs suffice.
      const SizeT rowWidth = 2 * inputWidth;
      const size_t innerScratch = std::max(KaratsubaMultiplication::ScratchLimbs(inputWidth, inputWidth),
                                           KaratsubaMultiplication::ScratchLimbs(inputWidth, inputWidth / 2));
      const size_t workSize = size_t{rowWidth} * (degree + 2) + size_t{inputWidth} * 4 + innerScratch;
      auto workspace = std::make_unique_for_overwrite<DataT[]>(workSize);
      Buffer remaining(workspace.get(), workSize);
      const auto take = [&](SizeT count) {
        const auto result = remaining.first(count);
        remaining = remaining.subspan(count);
        return result;
      };
      auto low = take(rowWidth), high = take(rowWidth);
      std::array<Buffer, evenCount> even;
      std::array<Buffer, oddCount> odd;
      for (auto &row : even)
        row = take(rowWidth);
      for (auto &row : odd)
        row = take(rowWidth);
      auto temporary = take(rowWidth);
      auto aPlus = take(inputWidth), aMinus = take(inputWidth);
      auto bPlus = take(inputWidth), bMinus = take(inputWidth);
      auto *scratch = remaining.data();
      assert(remaining.size() == innerScratch);

      MultiplyPoint(low, Piece(a, k, 0), Piece(b, k, 0), scratch);
      MultiplyPoint(high, Piece(a, k, LeftPieces - 1), Piece(b, k, RightPieces - 1), scratch);

      for (unsigned i = 0; i < evenCount; ++i)
      {
        const int power = static_cast<int>(i) - static_cast<int>(gridShift / 2);
        const bool aNegative = EvaluatePair(a, k, LeftPieces, power, aPlus, aMinus, temporary);
        const bool bNegative = EvaluatePair(b, k, RightPieces, power, bPlus, bMinus, temporary);
        MultiplyPoint(even[i], aPlus, bPlus, scratch);
        MultiplyPoint(odd[i], aMinus, bMinus, scratch);
        std::copy(even[i].begin(), even[i].end(), temporary.begin());
        auto *positive = even[i].data(), *negative = odd[i].data(), *plus = temporary.data();
        [[maybe_unused]] DataT carry, borrow;
        if (aNegative == bNegative)
        {
          carry = LimbKernel::AddN(positive, plus, negative, rowWidth);
          borrow = LimbKernel::SubN(negative, plus, negative, rowWidth);
        }
        else
        {
          borrow = LimbKernel::SubN(positive, plus, negative, rowWidth);
          carry = LimbKernel::AddN(negative, plus, negative, rowWidth);
        }
        assert(carry == 0 && borrow == 0);
        Shift(even[i], -1);
        Shift(odd[i], -1);

        // For x = v/u, the paired products were scaled by u^degree.
        // Remove the endpoints and the lowest power of x in each parity.
        // If Q is the remaining polynomial in x², store G^deg(Q)*Q(y/G)
        // at y = G*x² on the integer grid G = 2^gridShift.
        const int numerator = std::max(power, 0), denominator = std::max(-power, 0);
        SubtractMultiple(even[i], low, DataT{1} << (degree * denominator));
        SubtractMultiple(degree % 2 ? odd[i] : even[i], high, DataT{1} << (degree * numerator));
        Shift(even[i], static_cast<int>(gridShift * (evenCount - 1))
                           - static_cast<int>(degree - 2) * denominator - 2 * numerator);
        Shift(odd[i], static_cast<int>(gridShift * (oddCount - 1))
                          - static_cast<int>(degree - 1) * denominator - numerator);
      }

      Interpolate(even, gridShift);
      if constexpr (oddCount > evenCount)
      {
        // Toom-5 has one more odd coefficient. Evaluate at +4 and
        // subtract the now-known even contribution to get the last sample.
        constexpr int power = 2;
        EvaluateParts(a, k, LeftPieces, power, aPlus);
        EvaluateParts(b, k, RightPieces, power, bPlus);
        MultiplyPoint(odd.back(), aPlus, bPlus, scratch);
        SubtractMultiple(odd.back(), low, 1);
        SubtractMultiple(odd.back(), high, DataT{1} << (degree * power));
        for (unsigned i = 0; i < evenCount; ++i)
          SubtractMultiple(odd.back(), even[i], DataT{1} << (2 * (i + 1) * power));
        Shift(odd.back(), static_cast<int>(gridShift * (oddCount - 1)) - power);
      }
      Interpolate(odd, gridShift);

      std::vector<DataT> out(a.size() + b.size());
      const auto addCoefficient = [&](Limbs value, unsigned index) {
        const size_t offset = size_t{k} * index;
        const SizeT length = NonzeroLength(value);
        if (!length)
          return;
        assert(offset < out.size() && length <= out.size() - offset);
        [[maybe_unused]] const DataT carry = LimbKernel::AddInPlace(
            out.data() + offset, static_cast<SizeT>(out.size() - offset), value.data(), length);
        assert(carry == 0);
      };
      addCoefficient(low, 0);
      for (unsigned i = 0; i < oddCount; ++i)
        addCoefficient(odd[i], 2 * i + 1);
      for (unsigned i = 0; i < evenCount; ++i)
        addCoefficient(even[i], 2 * i + 2);
      addCoefficient(high, degree);
      TrimZeros(out);
      return out;
    }
  }

  class ToomHigherMultiplication
  {
  public:
    enum class Split { Toom5, Toom65, Toom85 };

    // One forced Toom level. Its point products reuse the existing recursive
    // limb multiplier. This entry point is intended for tests and tuning.
    static std::vector<DataT> Multiply(const std::vector<DataT> &a, const std::vector<DataT> &b,
                                       BaseT base, Split split)
    {
      if (IsZero(a) || IsZero(b))
        return {0};
      if (base != Base2_64)
        return KaratsubaMultiplication::Multiply(a, b, base);
      if (a.size() > static_cast<size_t>(std::numeric_limits<Int>::max()) ||
          b.size() > static_cast<size_t>(std::numeric_limits<Int>::max()) - a.size())
        throw std::length_error("Toom product exceeds the supported limb count");

      ToomHigherDetail::Limbs longer(a), shorter(b);
      if (longer.size() < shorter.size())
        std::swap(longer, shorter);
      switch (split)
      {
      case Split::Toom5: return ToomHigherDetail::MultiplySplit<5, 5>(longer, shorter);
      case Split::Toom65: return ToomHigherDetail::MultiplySplit<7, 6>(longer, shorter);
      case Split::Toom85: return ToomHigherDetail::MultiplySplit<9, 8>(longer, shorter);
      }
      throw std::invalid_argument("invalid Toom split");
    }
  };
}
