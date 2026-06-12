/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER
#define BIGINTEGER

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

    // Constructor, desctructor, and assignment operator
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
