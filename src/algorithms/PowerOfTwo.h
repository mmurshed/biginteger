#ifndef BIGMATH_POWER_OF_TWO_ARITHMETIC
#define BIGMATH_POWER_OF_TWO_ARITHMETIC

#include <algorithm>
#include <span>
#include <vector>

#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/common/LimbKernels.h"

namespace BigMath::detail::PowerOfTwo
{
  inline DataT InvertOddLimb(DataT odd, unsigned width = LimbBits)
  {
    DataT inverse = 1;
    // Each Newton step doubles the number of correct low bits.
    for (unsigned bits = 1; bits < width; bits *= 2)
      inverse *= 2 - odd * inverse;
    return inverse;
  }

  inline std::span<const DataT> Low(std::span<const DataT> value, std::size_t count)
  {
    return value.first(std::min(value.size(), count));
  }

  inline std::span<const DataT> High(std::span<const DataT> value, std::size_t count)
  {
    return value.subspan(std::min(value.size(), count));
  }

  inline void AddLow(std::span<DataT> destination, std::span<const DataT> source)
  {
    if constexpr (LimbBits == 64)
      LimbKernel::AddInPlace(destination.data(), destination.size(), source.data(), source.size());
    else
    {
      DataT carry = 0;
      for (std::size_t i = 0; i < destination.size(); ++i)
      {
        const DataT sum = destination[i] + carry + (i < source.size() ? source[i] : 0);
        destination[i] = sum & LimbMask;
        carry = sum >> 32;
      }
    }
  }

  inline std::vector<DataT> MultiplyLow(std::span<const DataT> a, std::span<const DataT> b,
                                       std::size_t count)
  {
    a = Low(a, count);
    b = Low(b, count);
    // Wide products can already use the transform dispatcher. Keep
    // the recursion to the range where cross products are cheaper.
    if (a.size() + b.size() <= count || count >= 1024)
    {
      auto product = detail::MultiplyLimbs(a, b, CurrentBase);
      if (product.size() > count)
        product.resize(count);
      return product;
    }

    if (count <= 32)
    {
      std::vector<DataT> product(count);
      for (std::size_t i = 0; i < a.size(); ++i)
      {
        const auto length = std::min(b.size(), count - i);
        DataT carry = 0;
        if constexpr (LimbBits == 64)
          carry = LimbKernel::AddMul1(product.data() + i, b.data(), length, a[i]);
        else
        {
          for (std::size_t j = 0; j < length; ++j)
          {
            const DataT sum = product[i + j] + a[i] * b[j] + carry;
            product[i + j] = sum & LimbMask;
            carry = sum >> 32;
          }
        }
        for (auto j = i + length; carry != 0 && j < count; ++j)
        {
          const auto sum = static_cast<ULong128>(product[j]) + carry;
          product[j] = static_cast<DataT>(sum) & LimbMask;
          carry = static_cast<DataT>(sum >> LimbBits);
        }
      }
      return product;
    }

    // With split = ceil(count/2), the two high halves never contribute
    // to the requested limbs. Both cross products can be truncated.
    const auto split = (count + 1) / 2;
    auto product = detail::MultiplyLimbs(Low(a, split), Low(b, split), CurrentBase);
    product.resize(count);
    auto tail = std::span<DataT>(product).subspan(split);
    AddLow(tail, MultiplyLow(High(a, split), b, tail.size()));
    AddLow(tail, MultiplyLow(a, High(b, split), tail.size()));
    return product;
  }

  inline void Negate(std::span<DataT> value)
  {
    DataT borrow = 0;
    for (auto &limb : value)
    {
      const DataT negated = (DataT{0} - limb - borrow) & LimbMask;
      borrow |= limb != 0;
      limb = negated;
    }
  }

  // Invert an odd magnitude modulo B^count. The requested count is positive
  // and may exceed the number of input limbs.
  inline std::vector<DataT> InvertOddLimbs(std::span<const DataT> odd, std::size_t count)
  {
    std::vector<DataT> inverse{InvertOddLimb(odd.front()) & LimbMask};
    inverse.reserve(count);
    while (inverse.size() < count)
    {
      const auto known = inverse.size();
      const auto extra = std::min(known, count - known);
      // odd*inverse = 1 + B^known*error. The next limbs of the
      // inverse are -inverse*error, truncated to 'extra' limbs.
      const auto product = detail::MultiplyLimbs(Low(odd, known + extra), inverse, CurrentBase);
      auto correction = MultiplyLow(inverse, High(product, known), extra);
      correction.resize(extra);
      Negate(correction);
      inverse.insert(inverse.end(), correction.begin(), correction.end());
    }
    return inverse;
  }
}

#endif
