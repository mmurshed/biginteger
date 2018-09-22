/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include <vector>
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
    // Default Constructor
    BigInteger() : theInteger(), isNegative(false) {}    

    BigInteger(SizeT bytes) : theInteger(bytes), isNegative(false)
    {}

    // Filled with specified data
    BigInteger(SizeT bytes, DataT fill) : theInteger(bytes), isNegative(false)
    {
      for(SizeT i = 0; i < bytes; i++)
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

private:
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

    // true if 'this' is zero
    bool IsZero() const
    {
      return size() == 0 || (size() == 1 && theInteger[0] == 0);
    }

    SizeT size() const 
    {
      return theInteger.size();
    }

private:
    // Implentation of addition by paper-pencil method
    static BigInteger& AddUnsigned(BigInteger const& a, BigInteger const& b)
    {
      SizeT size = max(a.size(),  b.size()) + 1;
      BigInteger& result = *new BigInteger(size);

      Long carry = 0;
      Long sum = 0;

      for(SizeT i = 0; i < size; i++)
      {
        sum = carry;
        
        if(i < a.size())
          sum += a.theInteger[i];
        
        if(i < b.size())
          sum += b.theInteger[i];

        result.theInteger[i] = sum % BASE;
        carry = sum / BASE;
      }

      result.theInteger[size-1] = carry;

      result.TrimZeros();

      return result;
    }

    // Implentation of subtraction by paper-pencil method
    // Assumption: a > b
    static BigInteger& SubtractUnsigned(BigInteger const& big, BigInteger const& small)
    {
      BigInteger& result = *new BigInteger(big.size() + 1);

      Long carry = 0;
      Long diff = 0;

      for(SizeT i = 0; i < big.size(); i++)
      {
        diff = big.theInteger[i] - carry;
        if(i < small.size())
        {
          diff -= small.theInteger[i];
        }

        carry = 0;
        if(diff < 0)
        {
          diff += BASE;
          carry = 1;
        }

        result.theInteger[i] = diff;
      }

      result.TrimZeros();

      return result;
    }

      // Implentation of multiplication by paper-pencil method
      // Classical algorithm
    static BigInteger& MultiplyUnsigned(BigInteger const& a, BigInteger const& b)
    {
      if(a.IsZero() || b.IsZero())
        return *new BigInteger(); // 0 times anything is zero

      BigInteger& result = *new BigInteger(a.size() + b.size() + 1, 0);

      ULong carry = 0;
      ULong multiply;
      
      for(SizeT i = 0; i < b.size(); i++)
      {
        carry = 0;
        for(SizeT j = 0 ; j < a.size(); j++)
        {
          multiply  = (ULong)a.theInteger[j];
          multiply *= (ULong)b.theInteger[i];
          multiply += (ULong)result.theInteger[i + j];
          multiply += carry;
          
          result.theInteger[i + j] = multiply % BASE;
          carry = multiply / BASE;
        }
        result.theInteger[i + a.size()] = carry;
      }

      result.TrimZeros();

      return result;
    }

    // Karatsuba Algorithm
    // static BigInteger& FastMultiplyUnsigned(BigInteger const& a, BigInteger const& b)
    // {
    //   BigInteger aa(a);
    //   BigInteger bb(b);
    //   SizeT size = max(a.size(), b.size());
    //   while(aa.size() < size)
    //     aa.theInteger.push_back(0);
    //   while(bb.size() < size)
    //     bb.theInteger.push_back(0);
    // }

    // static BigInteger& FastMultiplyUnsigned(BigInteger const& a, BigInteger const& b)
    // {
    // }

public:
    static BigInteger& Add(BigInteger const& a, BigInteger const& b)
    {
      // Check for zero
      bool aZero = a.IsZero();
      bool bZero = b.IsZero();
      if(aZero && bZero)
        return *new BigInteger(); // 0 + 0
      if(aZero)
        return *new BigInteger(b); // 0 + b
      if(bZero)
        return *new BigInteger(a); // a + 0

      // Check if actually subtraction is needed
      if(a.isNegative && !b.isNegative)
        return SubtractUnsigned(b, a); // a is negative, b is not. return b - a
      else if(!a.isNegative && b.isNegative)
        return SubtractUnsigned(a, b); // b is negative and a is not. return a - b

      // Add
      BigInteger& result = AddUnsigned(a, b);

      // Flip the sign when adding two negative numbers
      if(a.isNegative && b.isNegative)
        result.isNegative = true;

      return result;

    }

    // Straight pen-pencil implementation for subtraction
    static BigInteger& Subtract(BigInteger const& a, BigInteger const& b)
    {
      // Check for zero
      bool aZero = a.IsZero();
      bool bZero = b.IsZero();
      if(aZero && bZero)
        return *new BigInteger(); // 0 - 0
      if(aZero)
        return -(*new BigInteger(b)); // 0 - b
      if(bZero)
        return *new BigInteger(a); // a - 0

      // Check if actually addition is needed
      if(a.isNegative && !b.isNegative)
        return AddUnsigned(a, b).SetSign(true); // a is negative, b is not. return -(a + b)
      else if(!a.isNegative && b.isNegative)
        return AddUnsigned(a, b); // b is negative and a is not. return a + b

      Int cmp = a.UnsignedCompareTo(b);
      if(cmp == 0)
        return *new BigInteger(); // a == b, a-b == 0
      else if(cmp < 0)
        return SubtractUnsigned(b, a).SetSign(true); // -(b - a)
      else if (cmp > 0)
        return SubtractUnsigned(a, b); // a - b
      
      return *new BigInteger(); // Zero
    }
 
    static BigInteger& Multiply(BigInteger const& a, BigInteger const& b)
    {
      if(a.IsZero() || b.IsZero())
        return *new BigInteger(); // 0 times anything is zero

      BigInteger& result = MultiplyUnsigned(a, b);
      if(a.isNegative != b.isNegative)
        result.isNegative = true;
      
      return result;
    }
 
    // Negation, returns -*this
    BigInteger& operator-()
    {
      isNegative = !isNegative;
      return *this;
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

  // Adds Two BigInteger
  BigInteger& operator+(BigInteger const& a, BigInteger const& b)
  {
    return BigInteger::Add(a, b);
  }

  // Subtructs Two BigInteger
  BigInteger& operator-(BigInteger const& a, BigInteger const& b)
  {
    return BigInteger::Subtract(a, b);
  }

  // Multiplies Two BigInteger
  BigInteger& operator*(BigInteger const& a, BigInteger const& b)
  {
    return BigInteger::Multiply(a, b);
  }
}


