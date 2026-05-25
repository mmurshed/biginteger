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

#include "Util.h"
#include "../algorithms/Addition.h"
#include "../algorithms/multiplication/ClassicMultiplication.h"
#include "../algorithms/division/ClassicDivision.h"

namespace BigMath
{
  // Base 10^18 grouping for fast decimal <-> base-2^32 conversion.
  // 10^18 < 2^60 fits in uint64; products with a 32-bit limb need __int128.
  static constexpr ULong Base10_18 = 1000000000000000000ULL;
  static constexpr SizeT Base10_18_Zeroes = 18;

  inline Int FindRange(char const *num, Int &start, Int &end, bool &isNegative)
  {
    Int char_processed = start;
    isNegative = false;

    if (num[start] == '-')
    {
      isNegative = true;
      start++;
    }
    else if (num[start] == '+')
    {
      start++;
    }

    while (num[start] == '0')
      start++;

    end = start;
    while (num[end] >= '0' && num[end] <= '9')
      end++;

    char_processed = end - char_processed;
    end--;

    return char_processed;
  }

  // Parse decimal digit range [start..end] into base-2^32 little-endian limbs.
  // Reads 18 digits at a time into a uint64 chunk, then folds via
  // ClassicMultiplication::MultiplyTo + AddTo on the running result.
  inline vector<DataT> ParseUnsigned(char const *num, Int start, Int end)
  {
    vector<DataT> r;
    if (start > end)
      return r;

    SizeT len = (SizeT)(end - start + 1);
    r.reserve(len / 9 + 2);
    r.push_back(0);

    SizeT remainder = len % Base10_18_Zeroes;
    Int pos = start;

    if (remainder > 0)
    {
      ULong chunk = 0;
      for (SizeT i = 0; i < remainder; ++i)
        chunk = chunk * 10 + (ULong)(num[pos++] - '0');
      AddTo(r, chunk, Base2_32);
    }

    while (pos <= end)
    {
      ULong chunk = 0;
      for (SizeT i = 0; i < Base10_18_Zeroes; ++i)
        chunk = chunk * 10 + (ULong)(num[pos++] - '0');
      ClassicMultiplication::MultiplyTo(r, Base10_18, Base2_32);
      AddTo(r, chunk, Base2_32);
    }

    return r;
  }

  inline BigInteger Parse(char const *num, Int start, Int *char_processed = 0)
  {
    if (num == NULL || start < 0)
      return BigInteger();

    Int end;
    bool isNegative;
    Int processed = FindRange(num, start, end, isNegative);
    if (char_processed != 0)
      *char_processed = processed;

    if (start > end)
      return BigInteger();

    vector<DataT> bigInt = ParseUnsigned(num, start, end);

    return BigInteger(bigInt, isNegative);
  }

  inline BigInteger Parse(char const *num, Int *char_processed = 0)
  {
    return Parse(num, 0, char_processed);
  }

  // Direct conversion of a 64-bit value into base-2^32 limbs.
  inline vector<DataT> Convert(ULong n)
  {
    if (n == 0)
      return vector<DataT>{0};
    vector<DataT> r;
    while (n)
    {
      r.push_back((DataT)(n & 0xFFFFFFFFULL));
      n >>= 32;
    }
    return r;
  }

  inline string ToString(vector<DataT> const &bigInt, SizeT start, SizeT end, bool isNeg = false)
  {
    if (IsZero(bigInt, start, end))
      return string("0");

    vector<DataT> r(bigInt.begin() + start, bigInt.begin() + end + 1);
    while (r.size() > 1 && r.back() == 0)
      r.pop_back();

    vector<ULong> chunks;
    chunks.reserve(r.size() + 1);
    while (!(r.size() == 1 && r[0] == 0))
      chunks.push_back((ULong)ClassicDivision::DivModTo(r, Base10_18, Base2_32));

    string s;
    s.reserve(chunks.size() * Base10_18_Zeroes + 2);
    if (isNeg)
      s.push_back('-');

    char buf[20];
    char *bufEnd = buf + 20;

    // Most-significant chunk: no zero padding.
    auto it = chunks.rbegin();
    {
      ULong v = *it++;
      char *p = bufEnd;
      do
      {
        *--p = (char)('0' + (v % 10));
        v /= 10;
      } while (v);
      s.append(p, bufEnd - p);
    }

    // Remaining chunks: zero-pad to 18 digits.
    for (; it != chunks.rend(); ++it)
    {
      ULong v = *it;
      for (Int i = (Int)Base10_18_Zeroes - 1; i >= 0; --i)
      {
        buf[i] = (char)('0' + (v % 10));
        v /= 10;
      }
      s.append(buf, Base10_18_Zeroes);
    }

    return s;
  }

  inline string ToString(vector<DataT> const &bigInt, bool isNeg = false)
  {
    return ToString(bigInt, 0, bigInt.size() - 1, isNeg);
  }

  inline string ToString(BigInteger const &bigInt)
  {
    return ToString(bigInt.GetInteger(), bigInt.IsNegative());
  }
}

#endif
