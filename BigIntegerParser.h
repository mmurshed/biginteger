/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_PARSER
#define BIGINTEGER_PARSER

#include <vector>
#include <string>
using namespace std;

#include "BigInteger.h"
#include "algorithms/BigIntegerUtil.h"
#include "algorithms/classic/CommonAlgorithms.h"

namespace BigMath
{
  class BigIntegerParser
  {
  public:
    static vector<DataT> Convert(ULong n)
    {
      string num = to_string(n);
      return Parse(num.c_str()).GetInteger();
    }

    static Int FindRange(char const* num, Int& start, Int& end, bool& isNegative)
    {
      Int char_processed = start;
      isNegative = false;

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
      while(num[start] == '0')
      {
        start++;
      }

      // Find the end of the string
      end = start;
      while(num[end] >= '0' && num[end] <= '9' && num[end] != 0)
      {
        end++;
      }
      char_processed = end - char_processed;
      end--;

      return char_processed;
    }

    static BigInteger Parse(char const* num, Int* char_processed = 0)
    {
      return Parse(num, 0, char_processed);
    }

    static BigInteger Parse(char const* num, Int start, Int* char_processed = 0)
    {
      if(num == NULL || start < 0)
        return BigInteger();
      
      Int end;
      bool isNegative;
      Int processed = FindRange(num, start, end, isNegative);
      if(char_processed != 0)
        *char_processed = processed;

      // If the resulting string is empty return 0
      if (start > end)
        return BigInteger();

      vector<DataT> bigInt = ParseUnsigned(num, start, end);

      return BigInteger(bigInt, isNegative);
    }

    static vector<DataT> ParseUnsigned(char const* num, Int start, Int end)
    {
      vector<DataT> bigIntB1 = ParseUnsignedBase10n(
        num,
        start,
        end,
        BigIntegerUtil::Base100M_Zeroes);
      vector<DataT> bigIntB2 = CommonAlgorithms::ConvertBase(
        bigIntB1,
        BigIntegerUtil::Base100M,
        BigInteger::Base());
      return bigIntB2;
    }

    /*
    “little‐endian” means:
	•	bigInt[0] = coefficient of B^0  (least significant digit),
	•	bigInt[1] = coefficient of B^1,
	•	…
	•	bigInt[n - 1] = coefficient of B^{n-1} (most significant digit).å
    */
    // Group by n, meaning 10^n base
    static vector<DataT> ParseUnsignedBase10n(char const* num, Int start, Int end, SizeT digits)
    {
      vector<DataT> bigInt(end / digits + 1);

      SizeT j = 0;

      bool error = false;
      
      while(end >= start)
      {
        DataT digit = 0;
        DataT TEN = 1;
        
        for (SizeT i = 0; i < digits && end >= start; i++, end--)
        {
          // Current digit
          int d = num[end] - '0';
          if(d < 0 || d > 9)
          {
            error = true;
            break;
          }

          digit += TEN * d;
          TEN *= 10;
        }
 
        bigInt[j++] = digit;
        if(error) break;
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
      if(BigIntegerUtil::IsZero(bigInt, start, end))
      {
        return *new string("0");
      }

      vector<DataT> bigIntB2 = CommonAlgorithms::ConvertBase(
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
