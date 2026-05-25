/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_PARSER
#define BIGINTEGER_PARSER

#include <vector>
#include <string>
#include <unordered_map>
using namespace std;

#include "Util.h"
#include "../algorithms/Addition.h"
#include "../algorithms/Multiplication.h"
#include "../algorithms/Squaring.h"
#include "../algorithms/multiplication/ClassicMultiplication.h"
#include "../algorithms/division/ClassicDivision.h"
#include "../algorithms/division/NewtonDivision.h"

namespace BigMath
{
  // Base 10^18 grouping for fast decimal <-> base-2^32 conversion.
  // 10^18 < 2^60 fits in uint64; products with a 32-bit limb need ULong128.
  static constexpr ULong Base10_18 = 1000000000000000000ULL;
  static constexpr SizeT Base10_18_Zeroes = 18;
  static constexpr SizeT DecimalDcThreshold = 8192;

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

  inline ULong ParseChunk(char const *num, Int start, SizeT len)
  {
    ULong chunk = 0;
    for (SizeT i = 0; i < len; ++i)
      chunk = chunk * 10 + (ULong)(num[start + (Int)i] - '0');
    return chunk;
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

  inline vector<DataT> ParseUnsignedLinear(char const *num, Int start, Int end)
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
      ULong chunk = ParseChunk(num, pos, remainder);
      pos += (Int)remainder;
      AddTo(r, chunk, Base2_32);
    }

    while (pos <= end)
    {
      ULong chunk = ParseChunk(num, pos, Base10_18_Zeroes);
      pos += (Int)Base10_18_Zeroes;
      ClassicMultiplication::MultiplyTo(r, Base10_18, Base2_32);
      AddTo(r, chunk, Base2_32);
    }

    return r;
  }

  inline vector<DataT> Pow10(SizeT digits)
  {
    static thread_local unordered_map<SizeT, vector<DataT>> cache;

    auto it = cache.find(digits);
    if (it != cache.end())
      return it->second;

    vector<DataT> value;
    if (digits == 0)
    {
      value = vector<DataT>{1};
    }
    else if (digits <= Base10_18_Zeroes)
    {
      ULong p = 1;
      for (SizeT i = 0; i < digits; ++i)
        p *= 10;
      value = Convert(p);
    }
    else if (digits % 2 == 0)
    {
      // Even split: Pow10(d) = Pow10(d/2)². One Square (single FFT in NTT path)
      // instead of Multiply(p, p) (two FFTs).
      vector<DataT> p = Pow10(digits / 2);
      value = Square(p, Base2_32);
    }
    else
    {
      SizeT lo = digits / 2;
      SizeT hi = digits - lo;
      value = Multiply(Pow10(hi), Pow10(lo), Base2_32);
    }

    return cache.emplace(digits, value).first->second;
  }

  inline vector<DataT> ParseUnsignedDivideConquer(char const *num, Int start, Int end)
  {
    SizeT len = (SizeT)(end - start + 1);
    if (len <= DecimalDcThreshold)
      return ParseUnsignedLinear(num, start, end);

    SizeT lowDigits = len / 2;
    Int split = end - (Int)lowDigits;

    vector<DataT> high = ParseUnsignedDivideConquer(num, start, split);
    vector<DataT> low = ParseUnsignedDivideConquer(num, split + 1, end);
    vector<DataT> scale = Pow10(lowDigits);
    vector<DataT> scaledHigh = Multiply(high, scale, Base2_32);
    vector<DataT> result = Add(scaledHigh, low, Base2_32);
    TrimZerosToOne(result);
    return result;
  }

  // Parse decimal digit range [start..end] into base-2^32 little-endian limbs.
  // Small ranges use a simple 10^18 fold; large ranges use divide-and-conquer
  // so conversion benefits from Karatsuba/NTT multiplication.
  inline vector<DataT> ParseUnsigned(char const *num, Int start, Int end)
  {
    if (start > end)
      return vector<DataT>();
    return ParseUnsignedDivideConquer(num, start, end);
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

  // Linear divmod-10^18 formatter. Appends to `out`. If padTo > 0, pads to exactly
  // padTo decimal digits with leading zeros (used by the D&C lower halves).
  inline void ToStringLinearAppend(vector<DataT> r, SizeT padTo, string &out)
  {
    if (IsZero(r))
    {
      if (padTo == 0)
        out.push_back('0');
      else
        out.append(padTo, '0');
      return;
    }

    while (r.size() > 1 && r.back() == 0)
      r.pop_back();

    vector<ULong> chunks;
    chunks.reserve(r.size() + 1);
    while (!(r.size() == 1 && r[0] == 0))
      chunks.push_back((ULong)ClassicDivision::DivModTo(r, Base10_18, Base2_32));

    // Natural digit count: full 18 per chunk except the top one.
    int topDigits = 0;
    {
      ULong v = chunks.back();
      do { ++topDigits; v /= 10; } while (v);
    }
    SizeT natural = (chunks.size() - 1) * Base10_18_Zeroes + (SizeT)topDigits;

    if (padTo > natural)
      out.append(padTo - natural, '0');

    char buf[20];
    char *bufEnd = buf + 20;

    auto it = chunks.rbegin();
    {
      ULong v = *it++;
      char *p = bufEnd;
      do { *--p = (char)('0' + (v % 10)); v /= 10; } while (v);
      out.append(p, bufEnd - p);
    }
    for (; it != chunks.rend(); ++it)
    {
      ULong v = *it;
      for (Int i = (Int)Base10_18_Zeroes - 1; i >= 0; --i)
      {
        buf[i] = (char)('0' + (v % 10));
        v /= 10;
      }
      out.append(buf, Base10_18_Zeroes);
    }
  }

  // Threshold (in decimal digits) below which the linear formatter wins.
  // D&C pays a Newton-reciprocal setup per level; only amortizes for big inputs.
#ifndef BIGMATH_TOSTR_DC_THRESHOLD
#define BIGMATH_TOSTR_DC_THRESHOLD 2048
#endif
  static constexpr SizeT ToStringDcThreshold = BIGMATH_TOSTR_DC_THRESHOLD;

  // Chain entry: 10^digits + its precomputed Newton-reciprocal divider.
  // Built top-down so each level splits its parent in half — balanced T(N)=2T(N/2)+M(N).
  struct DecimalDcEntry
  {
    SizeT digits;
    vector<DataT> value;
    std::shared_ptr<NewtonDivision::Divider> divider;
  };

  // Chain[0] = N/2-digit divisor (top split). Chain[1] = N/4. ... down until level is
  // small enough to hand off to the linear base-case formatter.
  inline vector<DecimalDcEntry> BuildDecimalDcChain(SizeT topDigits)
  {
    vector<DecimalDcEntry> chain;
    for (SizeT d = topDigits; d >= ToStringDcThreshold / 2; d /= 2)
    {
      DecimalDcEntry e;
      e.digits = d;
      e.value = Pow10(d);
      e.divider = std::make_shared<NewtonDivision::Divider>(e.value, Base2_32);
      chain.push_back(std::move(e));
    }
    return chain;
  }

  // Recursive splitter. Writes `n` as decimals into `out`. padTo > 0 → pad with leading
  // zeros so the join is well-defined. Top-level call uses padTo = 0.
  inline void ToStringDivConquer(
      vector<DataT> n,
      vector<DecimalDcEntry> const &chain,
      SizeT level,
      SizeT padTo,
      string &out)
  {
    if (IsZero(n))
    {
      if (padTo == 0)
        out.push_back('0');
      else
        out.append(padTo, '0');
      return;
    }

    if (level >= chain.size())
    {
      ToStringLinearAppend(std::move(n), padTo, out);
      return;
    }

    // n smaller than this level's divisor → top half empty; skip down.
    if (Compare(n, chain[level].value) < 0)
    {
      ToStringDivConquer(std::move(n), chain, level + 1, padTo, out);
      return;
    }

    auto qr = chain[level].divider->DivideAndRemainder(n);
    SizeT half = chain[level].digits;
    SizeT topPad = (padTo > half) ? padTo - half : 0;
    ToStringDivConquer(std::move(qr.first), chain, level + 1, topPad, out);
    ToStringDivConquer(std::move(qr.second), chain, level + 1, half, out);
  }

  inline string ToString(vector<DataT> const &bigInt, SizeT start, SizeT end, bool isNeg = false)
  {
    if (IsZero(bigInt, start, end))
      return string("0");

    vector<DataT> r(bigInt.begin() + start, bigInt.begin() + end + 1);
    while (r.size() > 1 && r.back() == 0)
      r.pop_back();

    // Approximate decimal-digit count: each limb (32 bits) ≈ 9.633 digits.
    // Overestimate is fine — used only to pre-reserve and pick D&C level.
    SizeT approxDigits = (SizeT)((double)r.size() * 9.633) + 1;

    string s;
    s.reserve(approxDigits + (isNeg ? 2 : 1));
    if (isNeg)
      s.push_back('-');

    if (approxDigits < ToStringDcThreshold)
    {
      ToStringLinearAppend(std::move(r), 0, s);
      return s;
    }

    auto chain = BuildDecimalDcChain(approxDigits / 2);
    ToStringDivConquer(std::move(r), chain, 0, 0, s);
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
