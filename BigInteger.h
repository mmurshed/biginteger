/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_H
#define BIGINTEGER_H

#include <vector>
#include <tuple>
using namespace std;

namespace BigMath
{
  typedef long long Long;
  typedef int Int;
  typedef unsigned long long ULong;
  typedef unsigned int UInt;

  // The Size Type
  typedef UInt SizeT;
  // The Data Type
  typedef UInt DataT;

  // The Base Used
  // const DataT BASE = 2147483647l;
  const DataT BASE = 100000000lu;
  // An invalid data
  const DataT INVALIDDATA = 4294967295UL;
  // Number of digits in `BASE'
  const SizeT BASE_DIGIT = 8;

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
    
    // Character array constructor
    BigInteger(char const* n)
    {
      BigInteger& bigInt = Parse(n);
      theInteger = bigInt.theInteger;
      isNegative = bigInt.isNegative;
    }

    // Copy constructor
    BigInteger(BigInteger const& copy) : theInteger(copy.theInteger), isNegative(copy.isNegative) {}

    // The Destructor
    ~BigInteger() {}

  public:
    // true if 'this' is zero
    bool IsZero() const
    {
      return size() == 0 || (size() == 1 && theInteger[0] == 0);
    }

    SizeT size() const 
    {
      return theInteger.size();
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
    static BigInteger& Parse(char const* num)
    {
      Int l = strlen(num);
      if(num == NULL || l == 0)
        return *new BigInteger();

      Int start = 0;
      bool isNegative = false;

      // Take care of the sign
      if (num[start] == '-')
      {
        isNegative = true;
        start++;
      }
      else if (num[start] == '+')
      {
        start++;
      }

      // Trim leading zeros
      while (num[start] == '0')
      {
        start++;
      }

      // If the resulting string is empty return 0
      l -= start;
      if (l <= 0)
        return *new BigInteger();
      
      BigInteger& bigInt = *new BigInteger(l / BASE_DIGIT + 1);
      bigInt.isNegative = isNegative;

      // Convert the string to int
      SizeT j = 0;
      --l;
      while(l >= start)
      {
        DataT digit = 0;
        DataT TEN = 1;
        
        for (SizeT i = 0; i < BASE_DIGIT && l >= start; i++, l--)
        {
          // Current digit
          int d = num[l] - '0';
          if(d < 0 || d > 9)
          {
            // Error
          }

          digit += TEN * d;
          TEN *= 10;
        }
 
        bigInt.theInteger[j++] = digit;
      }

      if (bigInt.IsZero())
      {
        bigInt.isNegative = false;
      }

      bigInt.TrimZeros();

      return bigInt;
    }

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
    // Trims Leading Zeros
    void TrimZeros()
    {
      while(size() > 0 && theInteger[size() - 1] == 0)
      {
        theInteger.pop_back();
      }
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
      Long diff = (Long)size() - (Long)with.size();
      if(diff != 0)
        return diff;

      // Both ints have same number of digits
      Int cmp = 0;
      for(SizeT i = size() - 1; i >= 0; i--)
      {
        diff = (Long)theInteger[i] - (Long)with.theInteger[i];
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

    // Converts `this' to a string representation
    string& ToString() const
    {
      int len = size() * BASE_DIGIT + 2;
      char* num = new char[len];

      SizeT j = len - 1;
      
      num[j--] = 0;

      for(SizeT i = 0; i < size(); i++)
      {
        DataT n = theInteger[i];
        SizeT l = 0;
        while(l++ < BASE_DIGIT)
        {
          num[j--] = (n % 10) + '0';
          n /= 10;
        }
      }

      // Trim zeros
      while(num[j+1] == '0')
        j++;

      if(isNegative)
        num[j--] = '-';

      string& converted = *new string(num + j + 1);

      return converted;
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