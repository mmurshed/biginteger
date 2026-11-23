/**
 * BigMath: Generic base conversion implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/BaseConversion.h"

#include <bit>
#include <span>

namespace BigMath
{
  namespace
  {
    unsigned BinaryDigitWidth(BaseT base)
    {
      if (base == Base2_64)
        return 64;

      const DataT radix = static_cast<DataT>(base);
      return base > 0 && std::has_single_bit(radix)
                 ? static_cast<unsigned>(std::countr_zero(radix))
                 : 0;
    }

    std::vector<DataT> RepackBinaryDigits(std::span<const DataT> digits,
                                          unsigned sourceBits, unsigned targetBits)
    {
      while (!digits.empty() && digits.back() == 0)
        digits = digits.first(digits.size() - 1);
      if (digits.empty())
        return {0};

      std::vector<DataT> result;
      result.reserve((digits.size() * sourceBits + targetBits - 1) / targetBits);
      const DataT mask = static_cast<DataT>((ULong128{1} << targetBits) - 1);
      ULong128 buffer = 0;
      unsigned bufferedBits = 0;

      for (DataT digit : digits)
      {
        // At most 63 bits remain from the preceding source limb.
        buffer |= static_cast<ULong128>(digit) << bufferedBits;
        bufferedBits += sourceBits;
        while (bufferedBits >= targetBits)
        {
          result.push_back(static_cast<DataT>(buffer) & mask);
          buffer >>= targetBits;
          bufferedBits -= targetBits;
        }
      }
      if (bufferedBits != 0)
        result.push_back(static_cast<DataT>(buffer));
      while (result.size() > 1 && result.back() == 0)
        result.pop_back();
      return result;
    }

    std::vector<DataT> ConvertFromBase2_64(std::span<const DataT> digits, DataT targetBase)
    {
      std::vector<DataT> result{0};
      for (std::size_t i = digits.size(); i > 0; --i)
      {
        // One Horner step multiplies the result by 2^64 and adds this limb.
        DataT carry = digits[i - 1];
        for (DataT &digit : result)
        {
          const ULong128 value = (static_cast<ULong128>(digit) << 64) | carry;
          digit = static_cast<DataT>(value % targetBase);
          carry = static_cast<DataT>(value / targetBase);
        }
        while (carry != 0)
        {
          result.push_back(carry % targetBase);
          carry /= targetBase;
        }
      }
      return result;
    }

    std::vector<DataT> ConvertDigits(std::span<const DataT> digits, BaseT base1, BaseT base2)
    {
      if (base1 == base2)
        return std::vector<DataT>(digits.begin(), digits.end());

      const unsigned sourceBits = BinaryDigitWidth(base1);
      const unsigned targetBits = BinaryDigitWidth(base2);
      if (sourceBits != 0 && targetBits != 0)
        return RepackBinaryDigits(digits, sourceBits, targetBits);
      if (base1 == Base2_64)
        return ConvertFromBase2_64(digits, static_cast<DataT>(base2));

      std::vector<DataT> result{0};
      for (std::size_t i = digits.size(); i > 0; --i)
      {
        ClassicMultiplication::MultiplyTo(result, base1, base2);
        AddTo(result, digits[i - 1], base2);
      }
      return result;
    }
  }

  std::vector<DataT> ConvertBase(std::vector<DataT> const &bigIntB1,
                                 SizeT start, SizeT end,
                                 BaseT base1,
                                 BaseT base2)
  {
    if (bigIntB1.empty() || start > end)
      return ConvertDigits({}, base1, base2);
    const auto digits = std::span<const DataT>(bigIntB1).subspan(start, std::size_t{end} - start + 1);
    return ConvertDigits(digits, base1, base2);
  }

  std::vector<DataT> ConvertBase(std::vector<DataT> const &bigIntB1,
                                 BaseT base1,
                                 BaseT base2)
  {
    return ConvertDigits(bigIntB1, base1, base2);
  }
}
