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

    static BigInteger From(string const& num)
    {
      return From(num.c_str());
    }

    static BigInteger From(char const* num)
    {
      return BigIntegerParser::Parse(num);
    }
   };
}

#endif
