/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_BUILDER
#define BIGINTEGER_BUILDER

#include <vector>
#include <string>
using namespace std;

#include "algorithms/BigIntegerUtil.h"
#include "BigIntegerParser.h"
#include "BigInteger.h"

namespace BigMath
{
  class BigIntegerBuilder
  {
  public:
    static BigInteger From(ULong n)
    {
      return From(to_string(n));
    }

    static BigInteger From(Long n)
    {
      return From(to_string(n));
    }

    static BigInteger From(string const &num)
    {
      return From(num.c_str());
    }

    static BigInteger From(char const *num)
    {
      return BigIntegerParser::Parse(num);
    }

    // Convert from long double.
    static BigInteger From(long double value)
    {
      auto [theInteger, isNegative] = VectorFrom(value);
      return BigInteger(theInteger, isNegative);
    }

    static pair<vector<DataT>, bool> VectorFrom(long double value)
    {
      // If the value is zero, simply store one limb with zero.
      if (value == 0.0L)
      {
        return {vector<DataT>{0}, false};
      }

      bool isNegative = false;
      vector<DataT> theInteger;
      // Set the sign flag and work with the absolute value.
      if (value < 0.0L)
      {
        isNegative = true;
        value = -value;
      }

      // Define the base for one limb. For example, if each limb is 16 bits:
      const long double base = BigInteger::Base();

      // While there is still a nonzero part of the value, extract a limb.
      while (value > 0.0L)
      {
        // Get the remainder of value divided by base.
        long double rem = fmodl(value, base);
        // Convert the remainder to a limb.
        DataT limb = static_cast<DataT>(rem);
        theInteger.push_back(limb);
        // Divide value by base (using floor to remove the fractional part).
        value = floorl(value / base);
      }

      // Remove any (unlikely) extra zero limbs (from the most-significant side)
      // to keep the representation as compact as possible.
      BigIntegerUtil::TrimZeros(theInteger);

      return {theInteger, isNegative};
    }
  };
}

#endif
