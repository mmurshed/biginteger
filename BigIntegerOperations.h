/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_OPERATIONS
#define BIGINTEGER_OPERATIONS

#include <utility>
#include <vector>
using namespace std;

#include "BigInteger.h"
#include "BigIntegerComparator.h"
#include "algorithms/classic/ClassicAddition.h"
#include "algorithms/classic/ClassicSubtraction.h"
#include "algorithms/classic/ClassicMultiplication.h"
#include "algorithms/classic/ClassicDivision.h"
#include "algorithms/classic/CommonAlgorithms.h"
#include "algorithms/karatsuba/KaratsubaMultiplication.h"

namespace BigMath
{
  class BigIntegerOperations
  {
    // Karatsuba performs better for over 128 digits for the result
    static const SizeT KARATSUBA_THRESHOLD = 128;
    static const SizeT TOOM_COOK_THRESHOLD = 512;

    private:
    // Implentation of addition by paper-pencil method
    static BigInteger AddUnsigned(BigInteger const& a, BigInteger const& b)
    {
      return BigInteger(
        ClassicAddition::Add(
          a.GetInteger(),
          b.GetInteger(),
          BigInteger::Base())
      );
    }

    // Implentation of subtraction by paper-pencil method
    // Assumption: a > b
    static BigInteger SubtractUnsigned(BigInteger const& a, BigInteger const& b)
    {
      return BigInteger(
        ClassicSubtraction::Subtract(
          a.GetInteger(),
          b.GetInteger(),
          BigInteger::Base())
      );
    }

    static BigInteger MultiplyUnsigned(BigInteger const& a, BigInteger const& b)
    {
        return BigInteger(
          MultiplyUnsigned(
            a.GetInteger(),
            b.GetInteger(),
            BigInteger::Base())
        );
    }

    public:
    static vector<DataT> MultiplyUnsigned(vector<DataT> const& a, vector<DataT> const& b, ULong base)
    {
      SizeT size = a.size() + b.size();

      if(size <= KARATSUBA_THRESHOLD)
      {
        return ClassicMultiplication::Multiply(a, b, base);
      }

      return KaratsubaMultiplication::Multiply(a, b,base);
    }

public:
    static BigInteger Add(BigInteger const& a, BigInteger const& b)
    {
      // Check for zero
      bool aZero = a.IsZero();
      bool bZero = b.IsZero();
      if(aZero && bZero)
        return BigInteger(); // 0 + 0
      if(aZero)
        return BigInteger(b); // 0 + b
      if(bZero)
        return BigInteger(a); // a + 0

      // Check if actually subtraction is needed
      bool aNeg = a.IsNegative();
      bool bNeg = b.IsNegative();
      if(aNeg && !bNeg)
        return SubtractUnsigned(b, a); // a is negative, b is not. return b - a
      else if(!aNeg && bNeg)
        return SubtractUnsigned(a, b); // b is negative and a is not. return a - b

      // Add
      BigInteger result = AddUnsigned(a, b);

      // Flip the sign when adding two negative numbers
      if(aNeg && bNeg)
        result.SetSign(true);

      return result;
    }

    // Straight pen-pencil implementation for subtraction
    static BigInteger Subtract(BigInteger const& a, BigInteger const& b)
    {
      // Check for zero
      bool aZero = a.IsZero();
      bool bZero = b.IsZero();
      if(aZero && bZero)
        return BigInteger(); // 0 - 0
      if(aZero)
      {
        BigInteger result(b);
        result.SetSign(!b.IsNegative()); // 0 - b
        return result;
      }
      if(bZero)
        return BigInteger(a); // a - 0

      // Check if actually addition is needed
      bool aNeg = a.IsNegative();
      bool bNeg = b.IsNegative();
      if(aNeg && !bNeg)
        return AddUnsigned(a, b).SetSign(true); // a is negative, b is not. return -(a + b)
      else if(!aNeg && bNeg)
        return AddUnsigned(a, b); // b is negative and a is not. return a + b

      Int cmp = BigIntegerComparator::CompareTo(a.GetInteger(), b.GetInteger());
      if(cmp < 0) {
        BigInteger result = SubtractUnsigned(b, a);
        result.SetSign(true); // -(b - a)
        return result;
      }
      else if (cmp > 0){
        BigInteger result = SubtractUnsigned(a, b);
        return result; // a - b
      }
      
      return BigInteger(); // Zero when a == b
    }
 
    static BigInteger Multiply(BigInteger const& a, BigInteger const& b)
    {
      if(a.IsZero() || b.IsZero())
        return BigInteger(); // 0 times anything is zero

      BigInteger result = MultiplyUnsigned(a, b);
      if(a.IsNegative() != b.IsNegative())
        result.SetSign(true);
      
      return result;
    }

    static BigInteger Multiply(BigInteger const& a, ULong b)
    {
      if(a.IsZero() || b == 0)
        return BigInteger(); // 0 times anything is zero

      vector<DataT> m = ClassicMultiplication::Multiply(a.GetInteger(), b, BigInteger::Base());
      BigInteger result = BigInteger(m, false);
      if(a.IsNegative() != b < 0)
        result.SetSign(true);
      
      return result;
    }

    static pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const& a, BigInteger const& b)
    {
      if(a.IsZero() || b.IsZero())
      {
        BigInteger q = BigInteger();
        return make_pair(q,q); // case of 0
      }

      Int cmp = BigIntegerComparator::CompareTo(a.GetInteger(), b.GetInteger());
      if(cmp == 0)
      {
        vector<DataT> one(1);
        one[0] = 1;
        BigInteger q = BigInteger(one, a.IsNegative() || b.IsNegative());
        BigInteger r = BigInteger();
        return make_pair(q,r); // case of 0 // case of a/a
      }
      else if (cmp < 0)
      {
        BigInteger q = BigInteger();
        BigInteger r = BigInteger(a.GetInteger(), a.IsNegative() || b.IsNegative() );
        return make_pair(q,r); // case of a < b
      }

      // Now: a > b
      pair< vector<DataT>, vector<DataT> > result = ClassicDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());

      BigInteger q = BigInteger(result.first, false);
      BigInteger r = BigInteger(result.second, false);

      if(a.IsNegative() != b.IsNegative())
      {
        q.SetSign(true);
        r.SetSign(true);
      }
      
      return make_pair(q,r);
    }

    static BigInteger Divide(BigInteger const& a, ULong b)
    {
      if(a.IsZero() || b == 0)
      {
        return BigInteger();// case of 0
      }

      // Assume a > b
      vector<DataT> q = ClassicDivision::Divide(
        a.GetInteger(),
        b,
        BigInteger::Base());

      BigInteger result = BigInteger(q, false);

      if(a.IsNegative() != b < 0)
      {
        result.SetSign(true);
      }
      
      return result;
    }    
   };

  // Adds Two BigInteger
  BigInteger operator+(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerOperations::Add(a, b);
  }

  // Subtructs Two BigInteger
  BigInteger operator-(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerOperations::Subtract(a, b);
  }

  // Multiplies Two BigInteger
  BigInteger operator*(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerOperations::Multiply(a, b);
  }

  BigInteger operator*(BigInteger const& a, ULong b)
  {
    return BigIntegerOperations::Multiply(a, b);
  }

  // Divides Two BigInteger
  BigInteger operator/(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerOperations::DivideAndRemainder(a, b).first;
  }

  BigInteger operator/(BigInteger const& a, ULong const& b)
  {
    return BigIntegerOperations::Divide(a, b);
  }

  BigInteger operator%(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerOperations::DivideAndRemainder(a, b).second;
  }

  // Comparison operators
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

  BigInteger operator<<(BigInteger const& a, DataT b)
  {
    vector<DataT> result = CommonAlgorithms::ShiftLeft(a.GetInteger(), b);
    return BigInteger(result, a.IsNegative());
  }  
}

#endif