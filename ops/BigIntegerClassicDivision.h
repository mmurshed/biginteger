/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_CLASSIC_DIVISION
#define BIGINTEGER_CLASSIC_DIVISION

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../algorithms/classic/ClassicDivision.h"

namespace BigMath
{
  class BigIntegerClassicDivision
  {
  public:
    static pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, BigInteger const &b)
    {
      if (a.IsZero() || b.IsZero())
      {
        BigInteger q = BigInteger();
        return {q, q}; // case of 0
      }

      Int cmp = BigIntegerComparator::Compare(a.GetInteger(), b.GetInteger());
      if (cmp == 0)
      {
        vector<DataT> one(1, 1); // size 1, value 1
        BigInteger q = BigInteger(one, a.IsNegative() || b.IsNegative());
        BigInteger r = BigInteger();
        return {q, r}; // case of 0 // case of a/a
      }
      else if (cmp < 0)
      {
        BigInteger q = BigInteger();
        BigInteger r = BigInteger(a.GetInteger(), a.IsNegative() || b.IsNegative());
        return {q, r}; // case of a < b
      }

      // Now: a > b
      pair<vector<DataT>, vector<DataT>> result = ClassicDivision::DivideAndRemainder(a.GetInteger(), b.GetInteger(), BigInteger::Base());

      // Convert to BigInteger
      BigInteger q = BigInteger(result.first, a.IsNegative() || b.IsNegative());
      BigInteger r = BigInteger(result.second, a.IsNegative() || b.IsNegative());

      return {q, r};
    }
  };
}

#endif