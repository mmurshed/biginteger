/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_ADDITION
#define BIGINTEGER_ADDITION

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../algorithms/addition/ClassicAddition.h"
#include "BigIntegerSubtraction.h"

namespace BigMath
{
  class BigIntegerAddition
  {
  public:
    static BigInteger Add(BigInteger const &a, BigInteger const &b)
    {
      // Check for zero
      bool aZero = a.IsZero();
      bool bZero = b.IsZero();
      // 0 + 0
      if (aZero && bZero)
        return BigInteger();
      // 0 + b
      if (aZero)
        return BigInteger(b);
      // a + 0
      if (bZero)
        return BigInteger(a);

      // Check if actually subtraction is needed
      bool aNeg = a.IsNegative();
      bool bNeg = b.IsNegative();

      // a is negative, b is not. return b - a
      if (aNeg != bNeg)
      {
        return BigIntegerSubtraction::SubtractCompared(b, a);
      }
      // b is negative and a is not. return a - b
      else if (!aNeg && bNeg)
      {
        return BigIntegerSubtraction::SubtractCompared(a, b);
      }

      // Both are positive or both are negative
      BigInteger result = BigInteger(
          ClassicAddition::Add(
              a.GetInteger(),
              b.GetInteger(),
              BigInteger::Base()));

      // Flip the sign when adding two negative numbers
      if (aNeg && bNeg)
        result.SetSign(true);

      return result;
    }
  };

  // Adds Two BigInteger
  BigInteger operator+(BigInteger const &a, BigInteger const &b)
  {
    return BigIntegerAddition::Add(a, b);
  }
}

#endif