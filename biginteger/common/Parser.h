/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_PARSER
#define BIGINTEGER_PARSER

#include <vector>
#include <string>
#include <cassert>
using namespace std;

#include "../BigInteger.h"
#include "Util.h"
#include "../algorithms/BaseConversion.h"

namespace BigMath
{
  Int FindRange(char const *num, Int &start, Int &end, bool &isNegative)
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
    while (num[start] == '0')
    {
      start++;
    }

    // Find the end of the string
    end = start;
    while (num[end] >= '0' && num[end] <= '9' && num[end] != 0)
    {
      end++;
    }
    char_processed = end - char_processed;
    end--;

    return char_processed;
  }

  /*
  “little‐endian” means:
•	bigInt[0] = coefficient of B^0  (least significant digit),
•	bigInt[1] = coefficient of B^1,
•	…
•	bigInt[n - 1] = coefficient of B^{n-1} (most significant digit).å
  */
  // Group by n, meaning 10^n base
  vector<DataT> ParseUnsignedBase10n(const char *num, Int start, Int end, SizeT digits)
  {
    // Calculate the total number of digits (assumes num is a null-terminated string)
    SizeT len = end - start + 1;
    SizeT numBlocks = len / digits;
    SizeT remainder = len % digits;

    // Reserve space for full blocks plus one block for the remainder if any.
    vector<DataT> bigInt(numBlocks + (remainder ? 1 : 0));

    // Precompute powers of 10 for a block of the given digit length.
    // For instance, for digits == 8, power[0] = 1, power[1] = 10, ... power[7] = 10^7.
    DataT power[20]; // Ensure the array is large enough.
    power[0] = 1;
    for (SizeT i = 1; i < digits; i++)
    {
      power[i] = power[i - 1] * 10;
    }

    SizeT j = 0;
    bool error = false;

    // Process the remainder (if the number of digits isn't a multiple of `digits`)
    if (remainder > 0)
    {
      DataT digit = 0;
      // The remainder block is the most-significant digits.
      // They reside at positions: start to start + remainder - 1.
      for (SizeT i = 0; i < remainder; i++)
      {
        int d = num[start + i] - '0';
        assert(d >= 0 && d <= 9); // For performance, you might remove these in production.
        digit = digit * 10 + d;
      }
      // Place this block at the end of the vector (most significant block in little-endian)
      bigInt[numBlocks + (remainder ? 1 : 0) - 1] = digit;
    }

    // Process full blocks.
    // Each full block has exactly `digits` digits.
    // We process from the least-significant block upward.
    // The least-significant digit is at the end of the string.
    for (SizeT block = 0; block < numBlocks; block++)
    {
      DataT digit = 0;
      // Compute the block value by reading the last `digits` digits.
      // Compute the starting index for this block.
      SizeT blockEnd = end - block * digits;    // index of last digit in the block
      SizeT blockStart = blockEnd - digits + 1; // index of first digit in the block
      // Process the block in little-endian order:
      // i.e. the least-significant digit of the block comes first.
      for (SizeT i = 0; i < digits; i++)
      {
        int d = num[blockEnd - i] - '0';
        assert(d >= 0 && d <= 9);
        digit += d * power[i]; // Use precomputed multiplier for position i.
      }
      bigInt[block] = digit;
    }

    // The bigInt vector is built as little-endian:
    //   block 0 holds the least-significant digits.
    //   The remainder block (if exists) is the most-significant block.

    TrimZeros(bigInt);
    return bigInt;
  }

  vector<DataT> ParseUnsigned(char const *num, Int start, Int end)
  {
    vector<DataT> bigIntB1 = ParseUnsignedBase10n(
        num,
        start,
        end,
        Base100_Zeroes);

    if (Base100 == BigInteger::Base())
      return bigIntB1;

    vector<DataT> bigIntB2 = ConvertBase(
        bigIntB1,
        Base100,
        BigInteger::Base());
    return bigIntB2;
  }

  BigInteger Parse(char const *num, Int start, Int *char_processed = 0)
  {
    if (num == NULL || start < 0)
      return BigInteger();

    Int end;
    bool isNegative;
    Int processed = FindRange(num, start, end, isNegative);
    if (char_processed != 0)
      *char_processed = processed;

    // If the resulting string is empty return 0
    if (start > end)
      return BigInteger();

    vector<DataT> bigInt = ParseUnsigned(num, start, end);

    return BigInteger(bigInt, isNegative);
  }

  BigInteger Parse(char const *num, Int *char_processed = 0)
  {
    return Parse(num, 0, char_processed);
  }

  // Convert unsigned integer to base 10^n string
  Int UnsignedBase10nToString(vector<DataT> const &a, SizeT start, SizeT end, SizeT baseDigit, char *num, SizeT len)
  {
    Int j = len - 1;

    num[j--] = 0;

    // Trim zeros
    while (end >= start && a[end] == 0)
      end--;

    if (end < start)
    {
      // Means entire number is zero
      // so produce "0" in the string
      if (len > 1)
      {
        num[0] = '0';
        num[1] = '\0';
        return 0; // string starts at num[0]
      }
      return -1; // or handle the error if len=0
    }

    for (Int i = start; i <= end; i++)
    {
      DataT n = a[i];
      SizeT l = 0;
      while (l++ < baseDigit)
      {
        num[j--] = (n % 10) + '0';
        n /= 10;
      }
    }

    // Trim zeros
    while (num[j + 1] == '0')
      j++;
    return j;
  }

  vector<DataT> Convert(ULong n)
  {
    string num = to_string(n);
    return Parse(num.c_str()).GetInteger();
  }

  string ToString(vector<DataT> const &bigInt, SizeT start, SizeT end, bool isNeg = false)
  {
    if (IsZero(bigInt, start, end))
    {
      return string("0");
    }

    vector<DataT> bigIntB2 = ConvertBase(
        bigInt, start, end,
        BigInteger::Base(),
        Base100);

    Int len = (Int)bigIntB2.size() * Base100_Zeroes + 2;
    char *num = new char[len];

    Int j = UnsignedBase10nToString(
        bigIntB2, 0, bigIntB2.size() - 1,
        Base100_Zeroes,
        num, len);

    if (isNeg)
      num[j--] = '-';

    string converted(num + j + 1);

    delete[] num;

    return converted;
  }

  string ToString(vector<DataT> const &bigInt, bool isNeg = false)
  {
    return ToString(bigInt, 0, bigInt.size() - 1, isNeg);
  }

  // Converts the integer to a string representation
  string ToString(BigInteger const &bigInt)
  {
    return ToString(bigInt.GetInteger(), bigInt.IsNegative());
  }
}

#endif
