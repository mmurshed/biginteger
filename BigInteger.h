/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_H
#define BIGINTEGER_H

#include <vector>
using namespace std;

#include "BigIntegerUtil.h"

namespace BigMath
{
  class BigInteger
  {
  private:
    // The Integer array to hold the number
    vector<DataT> theInteger;
    // True if the number is negative
    bool isNegative;

  public:
    BigInteger(SizeT size = 0, bool negative = false) : theInteger(size), isNegative(negative) {}
    
    BigInteger(vector<DataT> aInt, bool negative) : theInteger(aInt), isNegative(negative) {}

    // Filled with specified data
    BigInteger(SizeT size, bool negative, DataT fill) : theInteger(size), isNegative(negative)
    {
      for(SizeT i = 0; i < size; i++)
        theInteger[i] = fill;
    }
    
    // Copy constructor
    BigInteger(BigInteger const& copy) : theInteger(copy.theInteger), isNegative(copy.isNegative) {}

    // The Destructor
    ~BigInteger() {}

      // Assignment Operator
    BigInteger& operator=(BigInteger const& arg)
    {
      if(this != &arg)
      {
        theInteger = arg.theInteger;
        isNegative = arg.isNegative;
      }
      return *this;
    }

  public:
    // true if 'this' is zero
    bool IsZero() const
    {
      return BigIntegerUtil::IsZero(theInteger);
    }

    SizeT size() const 
    {
      return theInteger.size();
    }

    static ULong Base()
    {
      return BigIntegerUtil::Base2_32;
    }

    bool IsNegative() const
    {
      return isNegative;
    }
 
    // Negation, returns -*this
    BigInteger& operator-()
    {
      isNegative = !isNegative;
      return *this;
    }

    vector<DataT> const& GetInteger() const
    {
      return theInteger;
    }

    DataT operator[] (const SizeT i) const
    {
      return theInteger[i];
    }

public:
    // Trims Leading Zeros
    void TrimZeros()
    {
      BigIntegerUtil::TrimZeros(theInteger);
    }

    BigInteger& SetSign(bool sign)
    {
      this->isNegative = sign;
      return *this;
    }

public:
    // Compares this with `with' irrespective of sign
    // Returns
    // 0 if equal
    // +value if this > with
    // -value if this < with
    Int UnsignedCompareTo(BigInteger const& with)const
    {
      // Case with zero
      bool isZero = this->IsZero();
      bool otherZero = with.IsZero();

      if(isZero && otherZero)
        return 0;
      else if(!isZero && otherZero)
        return 1;
      else if(isZero && !otherZero)
        return -1;

      // Different in size
      Long diff = size();
      diff -= with.size();
      if(diff != 0)
        return diff;

      // Both ints have same number of digits
      Int cmp = 0;
      for(SizeT i = size() - 1; i >= 0; i--)
      {
        diff = theInteger[i];
        diff -= with.theInteger[i];
        if(diff != 0)
          return diff;
      }

      return cmp;
    }

    // Compares this with `with'
    // Returns
    // 0 if equal
    // 1 if this>with
    // -1 if this<with
    Int CompareTo(BigInteger const& with)const
    {
      // Case 1: Positive , Negative
      if(!isNegative && with.isNegative)
        return 1;
      // Case 2: Negative, Positive
      else if(isNegative && !with.isNegative)
        return -1;

      Int cmp = UnsignedCompareTo(with);
      
      // Now, Both are Same Sign
      Int neg = 1;
      if(isNegative && with.isNegative)
        neg = -1;

      return cmp * neg;
    }
   };

  // Operators
  bool operator==(BigInteger const& a, BigInteger const& b)
  {
    return a.CompareTo(b) == 0;
  }

  bool operator!=(BigInteger const& a, BigInteger const& b)
  {
    return a.CompareTo(b) != 0;
  }

  bool operator>=(BigInteger const& a, BigInteger const& b)
  {
    return a.CompareTo(b) >= 0;
  }

  bool operator<=(BigInteger const& a, BigInteger const& b)
  {
    return a.CompareTo(b) <= 0;
  }

  bool operator>(BigInteger const& a, BigInteger const& b)
  {
    return a.CompareTo(b)>0;
  }

  bool operator<(BigInteger const& a, BigInteger const& b)
  {
    return a.CompareTo(b)<0;
  }
}

#endif