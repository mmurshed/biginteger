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
    static BigInteger& Parse(char const* num)
    {
      Int len = strlen(num);
      if(num == NULL || len == 0)
        return *new BigInteger();

      Int start = 0;
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
      if (len <= start)
        return *new BigInteger();

      vector<DataT>& bigIntB1 = ParseUnsignedBase10n(num, start, len, BigIntegerUtil::Base10n_Digit);
      vector<DataT>& bigIntB2 = ClassicalAlgorithms::ConvertBase(bigIntB1, BigIntegerUtil::Base10n, BigInteger::Base());

      return *new BigInteger(bigIntB2, isNegative);
    }

    // Group by n, meaning 10^n base
    static vector<DataT>& ParseUnsignedBase10n(char const* num, Int start, Int len, SizeT n)
    {
      len--;
      vector<DataT>& bigInt = *new vector<DataT>(len / n + 1);

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
    static string& ToString(BigInteger const& bigInt)
    {
      if(bigInt.IsZero())
      {
        return *new string("0");
      }

      vector<DataT>& bigIntB2 = ClassicalAlgorithms::ConvertBase(bigInt.GetInteger(), BigInteger::Base(), BigIntegerUtil::Base10n);
      
      Int len = bigIntB2.size() * BigIntegerUtil::Base10n_Digit + 2;
      char* num = new char[len];

      Int j = UnsignedBase10nToString(bigIntB2, BigIntegerUtil::Base10n_Digit, num, len);

      if(bigInt.IsNegative())
        num[j--] = '-';

      string& converted = *new string(num + j + 1);

      delete [] num;

      return converted;
    }    

    // Convert unsigned integer to base 10^n string
    static Int UnsignedBase10nToString(vector<DataT> const& bigInt, SizeT baseDigit, char *num, SizeT len)
    {
      SizeT j = len - 1;
      
      num[j--] = 0;

      for(Int i = 0; i < bigInt.size(); i++)
      {
        DataT n = bigInt[i];
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