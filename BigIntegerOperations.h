/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_OPERATIONS
#define BIGINTEGER_OPERATIONS

#include <vector>
using namespace std;

#include "BigInteger.h"
#include "BigIntegerComparator.h"
#include "algorithms/classic/ClassicAddition.h"
#include "algorithms/classic/ClassicSubtraction.h"
#include "algorithms/classic/ClassicMultiplication.h"
#include "algorithms/karatsuba/KaratsubaMultiplication.h"
#include "algorithms/toomcook/ToomCookMultiplication.h"

namespace BigMath
{
  class BigIntegerOperations
  {
    // Karatsuba performs better for over 128 digits for the result
    static const SizeT MULTIPLICATION_SWITCH = 128;

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
      SizeT size = a.size() + b.size();

      if(size <= MULTIPLICATION_SWITCH)
      {
        return BigInteger(
          ClassicMultiplication::Multiply(
            a.GetInteger(),
            b.GetInteger(),
            BigInteger::Base())
        );
      }

      return BigInteger(
        KaratsubaMultiplication::Multiply(
          a.GetInteger(),
          b.GetInteger(),
          BigInteger::Base())
      );

      // return BigInteger(
      //   ToomCookAlgorithm::Multiply(
      //     a.GetInteger(), 
      //     b.GetInteger(), 
      //     BigInteger::Base())
      // );      
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
      if(cmp < 0)
        return SubtractUnsigned(b, a).SetSign(true); // -(b - a)
      else if (cmp > 0)
        return SubtractUnsigned(a, b); // a - b
      
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

    static vector<BigInteger> DivideAndRemainder(BigInteger const& a, BigInteger const& b)
    {
      vector<BigInteger> results(2);
      if(a.IsZero() || b.IsZero())
      {
        results[0] = BigInteger();
        results[1] = BigInteger();
        return results; // case of 0
      }

      Int cmp = BigIntegerComparator::CompareTo(a.GetInteger(), b.GetInteger());
      if(cmp == 0)
      {
        vector<DataT> one(1);
        one[0] = 1;
        results[0] = BigInteger(one, a.IsNegative() || b.IsNegative());
        results[1] = BigInteger();
        return results; // case of a/a
      }
      else if (cmp < 0)
      {
        results[0] = BigInteger();
        results[1] = BigInteger(a.GetInteger(), a.IsNegative() || b.IsNegative() );
        return results; // case of a < b
      }

      // Now: a > b
      vector< vector<DataT> > result = ClassicalAlgorithms::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());

      results[0] = BigInteger(result[0], false);
      results[1] = BigInteger(result[1], false);

      if(a.IsNegative() != b.IsNegative())
      {
        results[0].SetSign(true);
        results[1].SetSign(true);
      }
      
      return results;
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

  BigInteger operator/(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerOperations::DivideAndRemainder(a, b)[0];
  }

  BigInteger operator%(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerOperations::DivideAndRemainder(a, b)[1];
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
}

#endif