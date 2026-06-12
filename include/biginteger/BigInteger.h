/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER
#define BIGINTEGER

#include <algorithm>
#include <span>
#include <vector>

#include "common/Util.h"
#include "common/Comparator.h"

namespace BigMath
{
  class BigInteger
  {
    // Data
  private:
    // The Integer array to hold the number
    std::vector<DataT> theInteger;
    // True if the number is negative
    bool isNegative;

    // Constructor, destructor, and assignment operator
  public:
    explicit BigInteger(SizeT size = 0, bool negative = false) : theInteger(size == 0 ? 1 : size, 0), isNegative(negative)
    {
      if (isNegative && Zero())
        isNegative = false;
    }

    BigInteger(std::vector<DataT> const &aInt, bool negative) : theInteger(aInt), isNegative(negative)
    {
      TrimZerosToOne(theInteger);
      if(negative && Zero())
      {
        isNegative = false;
      }
    }

    BigInteger(std::vector<DataT>&& aInt, bool negative) : theInteger(std::move(aInt)), isNegative(negative)
    {
      TrimZerosToOne(theInteger);
      if (isNegative && Zero())
      {
        isNegative = false;
      }
    }

    // Filled with specified data
    BigInteger(SizeT size, bool negative, DataT fill) : theInteger(size), isNegative(negative)
    {
      SetBit(theInteger, 0, size - 1, fill);
      TrimZerosToOne(theInteger);
      if (isNegative && Zero())
        isNegative = false;
    }

    // Rule of five, all defaulted. The previous user-declared copy members
    // suppressed the implicit move operations, so std::move(BigInteger)
    // deep-copied the limb vector.
    BigInteger(BigInteger const &) = default;
    BigInteger(BigInteger &&) noexcept = default;
    BigInteger &operator=(BigInteger const &) = default;
    BigInteger &operator=(BigInteger &&) noexcept = default;
    ~BigInteger() = default;

    // Accessors
  public:
    std::vector<DataT> const &GetInteger() const
    {
      return theInteger;
    }

    std::vector<DataT> ReleaseInteger()
    {
      std::vector<DataT> r = std::move(theInteger);
      theInteger = {0};
      isNegative = false;
      return r;
    }

    std::vector<uint8_t> ToByteArray(bool bigEndian = true) const
    {
      if (Zero())
      {
        return {};
      }

      constexpr SizeT limbBytes = LimbBits / 8;
      std::vector<uint8_t> bytes;
      bytes.reserve(theInteger.size() * limbBytes);

      for (DataT val : theInteger)
      {
        for (SizeT b = 0; b < limbBytes; ++b)
        {
          bytes.push_back(static_cast<uint8_t>((val >> (b * 8)) & 0xFF));
        }
      }

      while (!bytes.empty() && bytes.back() == 0)
      {
        bytes.pop_back();
      }

      if (bigEndian)
      {
        std::reverse(bytes.begin(), bytes.end());
      }

      return bytes;
    }

    static BigInteger FromByteArray(std::span<const uint8_t> bytes, bool negative, bool bigEndian = true)
    {
      if (bytes.empty())
      {
        return BigInteger();
      }

      constexpr SizeT limbBytes = LimbBits / 8;
      const size_t numLimbs = (bytes.size() + limbBytes - 1) / limbBytes;
      std::vector<DataT> limbs(numLimbs, 0);

      auto getByte = [&](size_t idx) -> uint8_t {
        if (bigEndian)
        {
          return bytes[bytes.size() - 1 - idx];
        }
        else
        {
          return bytes[idx];
        }
      };

      for (size_t i = 0; i < numLimbs; ++i)
      {
        DataT limbVal = 0;
        for (size_t b = 0; b < limbBytes; ++b)
        {
          size_t byteIdx = i * limbBytes + b;
          if (byteIdx < bytes.size())
          {
            limbVal |= (static_cast<DataT>(getByte(byteIdx)) << (b * 8));
          }
        }
        limbs[i] = limbVal;
      }

      return BigInteger(std::move(limbs), negative);
    }

    DataT operator[](const SizeT i) const
    {
      return theInteger[i];
    }

    // Properties
    SizeT size() const
    {
      return (SizeT)theInteger.size();
    }

    static BaseT Base()
    {
      // Resolves to Base2_64 by default (BIGMATH_LIMB_64=1); Base2_32 under -DBIGMATH_LIMB_64=0.
      return CurrentBase;
    }

    bool IsNegative() const
    {
      return isNegative;
    }

    bool Zero() const
    {
      return IsZero(theInteger);
    }

  public:
    // Trims Leading Zeros
    SizeT Trim()
    {
      return TrimZerosToOne(theInteger);
    }

    BigInteger &SetSign(bool sign)
    {
      this->isNegative = sign && !Zero();
      return *this;
    }

    // Negation, returns a negated copy; *this is unchanged
    BigInteger operator-() const
    {
      BigInteger r(*this);
      if (!r.Zero())
        r.isNegative = !r.isNegative;
      return r;
    }

  public:
    // Compares this with `with'
    // Returns
    // 0 if equal
    // 1 if this>with
    // -1 if this<with
    Int CompareTo(BigInteger const &with) const
    {
      // Case 1: Positive , Negative
      if (!isNegative && with.isNegative)
        return 1;
      // Case 2: Negative, Positive
      else if (isNegative && !with.isNegative)
        return -1;

      Int cmp = Compare(theInteger, with.theInteger);

      // Now, Both are Same Sign
      Int neg = 1;
      if (isNegative && with.isNegative)
        neg = -1;

      return cmp * neg;
    }
  };
}

#endif
