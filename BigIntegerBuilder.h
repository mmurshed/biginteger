/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIG_INTEGER_BUILDER_H
#define BIG_INTEGER_BUILDER_H

#include <vector>
#include <string>
using namespace std;

#include "BigIntegerUtil.h"
#include "BigIntegerParser.h"
#include "BigInteger.h"
#include "ClassicalAlgorithms.h"

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
