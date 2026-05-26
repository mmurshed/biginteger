/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_BUILDER
#define BIGINTEGER_BUILDER

#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "Util.h"
#include "Parser.h"
#include "../BigInteger.h"

namespace BigMath
{
  // Factory class for building BigInteger from primitive types.
  // All methods are static; instances aren't meant to be constructed.
  class BigIntegerBuilder
  {
  public:
    static std::vector<DataT> VectorFrom(ULong n)
    {
      return From(std::to_string(n)).GetInteger();
    }

    static BigInteger From(ULong n)              { return From(std::to_string(n)); }
    static BigInteger From(Long n)               { return From(std::to_string(n)); }
    static BigInteger From(std::string const &n) { return From(n.c_str()); }
    static BigInteger From(char const *num)      { return Parse(num); }

    // Convert from long double.
    static BigInteger From(long double value)
    {
      auto [theInteger, isNegative] = VectorFrom(value);
      return BigInteger(theInteger, isNegative);
    }

    // Decompose long double into (little-endian limb vector, sign).
    // Limb base comes from BigInteger::Base() to stay aligned with the rest of
    // the library (currently 2^32).
    static std::pair<std::vector<DataT>, bool> VectorFrom(long double value)
    {
      if (value == 0.0L)
        return {std::vector<DataT>{0}, false};

      bool isNegative = false;
      std::vector<DataT> theInteger;
      if (value < 0.0L)
      {
        isNegative = true;
        value = -value;
      }

      // For Base2_64 the BaseT sentinel is 0 (since 2^64 doesn't fit in BaseT);
      // use ldexpl to recover the numeric value. For other bases the BaseT
      // value is the numeric base.
      const long double base = (BigInteger::Base() == Base2_64)
                                   ? ldexpl(1.0L, 64)
                                   : static_cast<long double>(BigInteger::Base());

      while (value > 0.0L)
      {
        long double rem = fmodl(value, base);
        DataT limb = static_cast<DataT>(rem);
        theInteger.push_back(limb);
        value = floorl(value / base);
      }

      TrimZeros(theInteger);
      return {theInteger, isNegative};
    }
  };
}

#endif
