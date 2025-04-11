/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_SUBTRACTION
#define BIGINTEGER_SUBTRACTION

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../BigIntegerBuilder.h"
#include "../BigIntegerComparator.h"
#include "../algorithms/classic/ClassicAddition.h"
#include "../algorithms/classic/ClassicSubtraction.h"

namespace BigMath
{
  class BigIntegerSubtraction
  {
    public:
    // Implentation of subtraction by paper-pencil method
    static BigInteger SubtractCompared(BigInteger const& a, BigInteger const& b)
    {
      int cmp = BigIntegerComparator::CompareTo(a.GetInteger(), b.GetInteger());
      if(cmp > 0) {
        return BigInteger(
          ClassicSubtraction::Subtract(
          a.GetInteger(),
          b.GetInteger(),
          BigInteger::Base())
        );
      }
      else if (cmp < 0) {
        BigInteger result = BigInteger(
          ClassicSubtraction::Subtract(
            b.GetInteger(),
            a.GetInteger(),
            BigInteger::Base())
        );  
        result.SetSign(true);
        return result;
      }
      return BigInteger(); // Zero when a == b
    }

  public:
    // Straight pen-pencil implementation for subtraction
    static BigInteger Subtract(BigInteger const& a, BigInteger const& b)
    {
      // Check for zero
      bool aZero = a.IsZero();
      bool bZero = b.IsZero();
      
      // 0 - 0
      if(aZero && bZero)
        return BigInteger();

      // 0 - b
      if(aZero)
      {
        BigInteger result(b);
        result.SetSign(!b.IsNegative());
        return result;
      }

      // a - 0
      if(bZero)
        return BigInteger(a);

      // Check if actually addition is needed
      bool aNeg = a.IsNegative();
      bool bNeg = b.IsNegative();
      
      // either a or b is negative
      if(aNeg != bNeg)
      {
        BigInteger result = BigInteger(
          ClassicAddition::Add(
            a.GetInteger(),
            b.GetInteger(),
            BigInteger::Base())
        );
        // if a is negative and b is not, result is -(a + b)
        if(!aNeg && bNeg)
          result.SetSign(true);
        return result;
      }

      return SubtractCompared(a, b); // Zero when a == b
    } 
   };

  // Subtructs Two BigInteger
  BigInteger operator-(BigInteger const& a, BigInteger const& b)
  {
    return BigIntegerSubtraction::Subtract(a, b);
  }

  BigInteger operator-(BigInteger const& a, DataT b)
  {
    return BigIntegerSubtraction::Subtract(a, BigIntegerBuilder::From(b));
  }

}

#endif