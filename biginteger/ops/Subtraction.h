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
#include "../common/Builder.h"
#include "../common/Comparator.h"
#include "../algorithms/Addition.h"
#include "../algorithms/Subtraction.h"

namespace BigMath
{
  // Implentation of subtraction by paper-pencil method
  BigInteger SubtractCompared(BigInteger const &a, BigInteger const &b)
  {
    int cmp = Compare(a.GetInteger(), b.GetInteger());
    if (cmp > 0)
    {
      return BigInteger(
          Subtract(
              a.GetInteger(),
              b.GetInteger(),
              BigInteger::Base()),
          false);
    }
    else if (cmp < 0)
    {
      BigInteger result = BigInteger(
          Subtract(
              b.GetInteger(),
              a.GetInteger(),
              BigInteger::Base()),
          true);
      return result;
    }
    return BigInteger(); // Zero when a == b
  }

  // Straight pen-pencil implementation for subtraction
  BigInteger Subtract(BigInteger const &a, BigInteger const &b)
  {
    // Check for zero
    bool aZero = a.Zero();
    bool bZero = b.Zero();

    // 0 - 0
    if (aZero && bZero)
      return BigInteger();

    // 0 - b
    if (aZero)
    {
      BigInteger result(b);
      result.SetSign(!b.IsNegative());
      return result;
    }

    // a - 0
    if (bZero)
      return BigInteger(a);

    // Check if actually addition is needed
    bool aNeg = a.IsNegative();
    bool bNeg = b.IsNegative();

    // either a or b is negative
    if (aNeg != bNeg)
    {
      BigInteger result = BigInteger(
          Add(
              a.GetInteger(),
              b.GetInteger(),
              BigInteger::Base()),
          aNeg && !bNeg); // if a is negative and b is not, result is -(a + b), negative
      return result;
    }

    return SubtractCompared(a, b); // Zero when a == b
  }

  // Subtructs Two BigInteger
  BigInteger operator-(BigInteger const &a, BigInteger const &b)
  {
    return Subtract(a, b);
  }

  BigInteger operator-(BigInteger const &a, DataT b)
  {
    return Subtract(a, BigIntegerBuilder::From(b));
  }

}

#endif