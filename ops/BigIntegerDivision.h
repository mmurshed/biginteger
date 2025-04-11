/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_DIVISION
#define BIGINTEGER_DIVISION

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../algorithms/classic/ClassicDivision.h"
#include "../algorithms/newtonraphson/NewtonRaphsonDivision.h"

namespace BigMath
{
  class BigIntegerDivision
  {
  public:
    static pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &a, BigInteger const &b)
    {
      if (a.IsZero() || b.IsZero())
      {
        BigInteger q = BigInteger();
        return make_pair(q, q); // case of 0
      }

      Int cmp = BigIntegerComparator::CompareTo(a.GetInteger(), b.GetInteger());
      if (cmp == 0)
      {
        vector<DataT> one(1, 1); // size 1, value 1
        BigInteger q = BigInteger(one, a.IsNegative() || b.IsNegative());
        BigInteger r = BigInteger();
        return make_pair(q, r); // case of 0 // case of a/a
      }
      else if (cmp < 0)
      {
        BigInteger q = BigInteger();
        BigInteger r = BigInteger(a.GetInteger(), a.IsNegative() || b.IsNegative());
        return make_pair(q, r); // case of a < b
      }

      // Now: a > b
      pair<BigInteger, BigInteger> result = NewtonRaphsonDivision::DivideAndRemainder(a, b);

      if (a.IsNegative() != b.IsNegative())
      {
        result.first.SetSign(true);
        result.second.SetSign(true);
      }

      return result;
    }
  };

  // Divides Two BigInteger
  BigInteger operator/(BigInteger const &a, BigInteger const &b)
  {
    return BigIntegerDivision::DivideAndRemainder(a, b).first;
  }

  BigInteger operator%(BigInteger const &a, BigInteger const &b)
  {
    return BigIntegerDivision::DivideAndRemainder(a, b).second;
  }

}

#endif