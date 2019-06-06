/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIG_INTEGER_PARSER_H
#define BIG_INTEGER_PARSER_H

#include <vector>
#include <string>
using namespace std;

#include "BigIntegerUtil.h"
#include "BigInteger.h"
#include "ClassicalAlgorithms.h"

namespace BigMath
{
  class BigIntegerParser
  {
  public:
    static vector<DataT> Convert(ULong n)
    {
      string num = to_string(n);
      return ParseUnsigned(num.c_str(), 0, num.length());
    }

    static BigInteger Parse(char const* num)
    {
      Int len = (Int)strlen(num);
      return ParseSigned(num, 0, len);
    }

    static BigInteger ParseSigned(char const* num, Int start, Int len)
    {
      if(num == NULL || len == 0 || start < 0 || start >= len)
        return BigInteger();

      bool isNegative = false;

      // Take care of the sign
      if (num[start] == '-')
      {
        isNegative = true;
        start++;
      }
      else if (num[start] == '+')
      {
        start++;
      }

      // Trim leading zeros
      while (num[start] == '0')
      {
        start++;
      }

      // If the resulting string is empty return 0
      if (start >= len)
        return BigInteger();

      vector<DataT> bigInt = ParseUnsigned(num, start, len);

      return BigInteger(bigInt, isNegative);
    }

    static vector<DataT> ParseUnsigned(char const* num, Int start, Int len)
    {
      vector<DataT> bigIntB1 = ParseUnsignedBase10n(num, start, len, BigIntegerUtil::Base100M_Zeroes);
      vector<DataT> bigIntB2 = ClassicalAlgorithms::ConvertBase(bigIntB1, BigIntegerUtil::Base100M, BigInteger::Base());
      return bigIntB2;
    }

    // Group by n, meaning 10^n base
    static vector<DataT> ParseUnsignedBase10n(char const* num, Int start, Int len, SizeT n)
    {
      len--;
      vector<DataT> bigInt(len / n + 1);

      // Convert the string to int
      SizeT j = 0;
      len -= start;
      
      while(len >= start)
      {
        DataT digit = 0;
        DataT TEN = 1;
        
        for (SizeT i = 0; i < n && len >= start; i++, len--)
        {
          // Current digit
          int d = num[len] - '0';
          if(d < 0 || d > 9)
          {
            // Error
          }

          digit += TEN * d;
          TEN *= 10;
        }
 
        bigInt[j++] = digit;
      }
      
      BigIntegerUtil::TrimZeros(bigInt);

      return bigInt;
    }

    // Converts the integer to a string representation
    static string ToString(BigInteger const& bigInt)
    {
      return ToString(bigInt.GetInteger(), bigInt.IsNegative());
    }

    static string ToString(vector<DataT> const& bigInt, bool isNeg = false)
    {
      return ToString(bigInt, 0, bigInt.size() - 1, isNeg);
    }

    static string ToString(vector<DataT> const& bigInt, SizeT start, SizeT end, bool isNeg = false)
    {
      if(BigIntegerUtil::IsZero(bigInt))
      {
        return *new string("0");
      }

      vector<DataT> bigIntB2 = ClassicalAlgorithms::ConvertBase(
        bigInt, start, end,
        BigInteger::Base(),
        BigIntegerUtil::Base100M);
      
      Int len = (Int)bigIntB2.size() * BigIntegerUtil::Base100M_Zeroes + 2;
      char* num = new char[len];

      Int j = UnsignedBase10nToString(
        bigIntB2, 0, bigIntB2.size() - 1,
        BigIntegerUtil::Base100M_Zeroes,
        num, len);

      if(isNeg)
        num[j--] = '-';

      string converted(num + j + 1);

      delete [] num;

      return converted;
    }

    // Convert unsigned integer to base 10^n string
    static Int UnsignedBase10nToString(vector<DataT> const& a, SizeT start, SizeT end, SizeT baseDigit, char *num, SizeT len)
    {
      Int j = len - 1;
      
      num[j--] = 0;

      // Trim zeros
      while(end >= start && a[end] == 0)
        end--;

      for(Int i = start; i <= end; i++)
      {
        DataT n = a[i];
        SizeT l = 0;
        while(l++ < baseDigit)
        {
          num[j--] = (n % 10) + '0';
          n /= 10;
        }
      }

      // Trim zeros
      while(num[j+1] == '0')
        j++;
      return j;
    }    

   };
}

#endif
