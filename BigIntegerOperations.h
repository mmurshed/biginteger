/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGEROPERATIONS_H
#define BIGINTEGEROPERATIONS_H

#include <vector>
#include <tuple>
using namespace std;

#include "BigIntegerUtil.h"
#include "BigInteger.h"
#include "ClassicalAlgorithms.h"

namespace BigMath
{
  class BigIntegerOperations
  {
    public:
    // Implentation of addition by paper-pencil method
    static BigInteger& AddUnsigned(BigInteger const& a, BigInteger const& b)
    {
      vector<DataT>& result = ClassicalAlgorithms::AddUnsigned(a.GetInteger(), b.GetInteger(), BigInteger::Base());
      return *new BigInteger(result, false);
    }

    // Implentation of subtraction by paper-pencil method
    // Assumption: a > b
    static BigInteger& SubtractUnsigned(BigInteger const& a, BigInteger const& b)
    {
      vector<DataT>& result = ClassicalAlgorithms::SubtractUnsigned(a.GetInteger(), b.GetInteger(), BigInteger::Base());
      return *new BigInteger(result, false);
    }

    static BigInteger& MultiplyUnsigned(BigInteger const& a, BigInteger const& b)
    {
      vector<DataT>& result = ClassicalAlgorithms::MultiplyUnsigned(a.GetInteger(), b.GetInteger(), BigInteger::Base());
      return *new BigInteger(result, false);
    }

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
      bool aNeg = a.IsNegative();
      bool bNeg = b.IsNegative();
      if(aNeg && !bNeg)
        return SubtractUnsigned(b, a); // a is negative, b is not. return b - a
      else if(!aNeg && bNeg)
        return SubtractUnsigned(a, b); // b is negative and a is not. return a - b

      // Add
      BigInteger& result = AddUnsigned(a, b);

      // Flip the sign when adding two negative numbers
      if(aNeg && bNeg)
        result.SetSign(true);

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
      bool aNeg = a.IsNegative();
      bool bNeg = b.IsNegative();
      if(aNeg && !bNeg)
        return AddUnsigned(a, b).SetSign(true); // a is negative, b is not. return -(a + b)
      else if(!aNeg && bNeg)
        return AddUnsigned(a, b); // b is negative and a is not. return a + b

      Int cmp = a.UnsignedCompareTo(b);
      if(cmp < 0)
        return SubtractUnsigned(b, a).SetSign(true); // -(b - a)
      else if (cmp > 0)
        return SubtractUnsigned(a, b); // a - b
      
      return *new BigInteger(); // Zero when a == b
    }
 
    static BigInteger& Multiply(BigInteger const& a, BigInteger const& b)
    {
      if(a.IsZero() || b.IsZero())
        return *new BigInteger(); // 0 times anything is zero

      BigInteger& result = MultiplyUnsigned(a, b);
      if(a.IsNegative() != b.IsNegative())
        result.SetSign(true);
      
      return result;
    }

    static vector<BigInteger>& DivideAndRemainder(BigInteger const& a, BigInteger const& b)
    {
      vector<BigInteger>& results = *new vector<BigInteger>(2);
      if(a.IsZero() || b.IsZero())
      {
        results[0] = *new BigInteger();
        results[1] = *new BigInteger();
        return results; // case of 0
      }

      Int cmp = a.UnsignedCompareTo(b);
      if(cmp == 0)
      {
        vector<DataT>& one = *new vector<DataT>(1);
        one[0] = 1;
        results[0] = *new BigInteger(one, a.IsNegative() || b.IsNegative());
        results[1] = *new BigInteger();
        return results; // case of a/a
      }
      else if (cmp < 0)
      {
        results[0] = *new BigInteger();
        results[1] = *new BigInteger(a.GetInteger(), a.IsNegative() || b.IsNegative() );
        return results; // case of a < b
      }

      // Now: a > b
      vector< vector<DataT> >& result = ClassicalAlgorithms::DivideAndRemainderUnsigned(a.GetInteger(), b.GetInteger(), BigInteger::Base());
      results[0] = *new BigInteger(result[0], false);
      results[1] = *new BigInteger(result[1], false);

      if(a.IsNegative() != b.IsNegative())
      {
        results[0].SetSign(true);
        results[1].SetSign(true);
      }
      
      return results;
    }
 
   };

     // Adds Two BigInteger
  BigInteger& operator+(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerOperations::Add(a, b);
  }

  // Subtructs Two BigInteger
  BigInteger& operator-(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerOperations::Subtract(a, b);
  }

  // Multiplies Two BigInteger
  BigInteger& operator*(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerOperations::Multiply(a, b);
  }

  BigInteger& operator/(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerOperations::DivideAndRemainder(a, b)[0];
  }

  BigInteger& operator%(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerOperations::DivideAndRemainder(a, b)[1];
  }
}

#endif