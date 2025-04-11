/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER
#define BIGINTEGER

#include <vector>
using namespace std;

#include "algorithms/BigIntegerUtil.h"
#include "BigIntegerComparator.h"

namespace BigMath
{
  class BigInteger
  {
    // Data
  private:
    // The Integer array to hold the number
    vector<DataT> theInteger;
    // True if the number is negative
    bool isNegative;

    // Constructor, desctructor, and assignment operator
  public:  
    BigInteger(SizeT size = 0, bool negative = false) : theInteger(size), isNegative(negative) {}
    
    BigInteger(vector<DataT> const& aInt, bool negative = false) : theInteger(aInt), isNegative(negative)
    {}

    // Filled with specified data
    BigInteger(SizeT size, bool negative, DataT fill) : theInteger(size), isNegative(negative)
    {
      BigIntegerUtil::SetBit(theInteger, 0, size - 1, fill);
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

  // Accessors
  public:
    vector<DataT> const& GetInteger() const
    {
      return theInteger;
    }

    DataT operator[] (const SizeT i) const
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
      return BigIntegerUtil::Base100;
    }

    bool IsNegative() const
    {
      return isNegative;
    }

    bool IsZero() const
    {
      return BigIntegerUtil::IsZero(theInteger);
    }
 
public:
    // Trims Leading Zeros
    SizeT TrimZeros()
    {
      return BigIntegerUtil::TrimZeros(theInteger);
    }

    BigInteger& SetSign(bool sign)
    {
      this->isNegative = sign;
      return *this;
    }

    // Negation, returns -*this
    BigInteger& operator-()
    {
      isNegative = !isNegative;
      return *this;
    }

    long double ToDouble() const {
        const long double base = Base(); // 2^16
        long double result = 0.0L;
        
        // Process digits from the most significant to the least significant.
        // The digits are stored in little-endian order, so we iterate in reverse.
        for (size_t i = theInteger.size(); i > 0; i--) {
            result = result * base + static_cast<long double>(theInteger[i - 1]);
        }
        
        // Apply the sign.
        if (isNegative)
            result = -result;
        
        return result;
    }    
public:
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

      Int cmp = BigIntegerComparator::CompareTo(theInteger, with.theInteger);
      
      // Now, Both are Same Sign
      Int neg = 1;
      if(isNegative && with.isNegative)
        neg = -1;

      return cmp * neg;
    }
   };

  BigInteger const& ONE = BigInteger(1, false, 1);
  BigInteger const& ZERO = BigInteger();
}

#endif
