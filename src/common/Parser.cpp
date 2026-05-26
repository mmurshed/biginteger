/**
 * BigMath: Decimal parsing and formatting implementation.
 *
 * Architectural overview: STRING_CONVERSION.md.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/common/Parser.h"
#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/Squaring.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/algorithms/division/ClassicDivision.h"
#include "biginteger/algorithms/division/NewtonDivision.h"

#include <memory>
#include <cmath>
#include <unordered_map>
#include <utility>

namespace BigMath
{
  Int FindRange(char const *num, Int &start, Int &end, bool &isNegative)
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

  ULong ParseChunk(char const *num, Int start, SizeT len)
  {
    ULong chunk = 0;
    for (SizeT i = 0; i < len; ++i)
      chunk = chunk * 10 + (ULong)(num[start + (Int)i] - '0');
    return chunk;
  }

  std::vector<DataT> Convert(ULong n)
  {
    if (n == 0)
      return std::vector<DataT>{0};
#if BIGMATH_LIMB_64
    // 64-bit limbs fit any ULong in one limb.
    return std::vector<DataT>{(DataT)n};
#else
    std::vector<DataT> r;
    while (n)
    {
      r.push_back((DataT)(n & 0xFFFFFFFFULL));
      n >>= 32;
    }
    return r;
#endif
  }

  std::vector<DataT> ParseUnsignedLinear(char const *num, Int start, Int end)
  {
    std::vector<DataT> r;
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
      AddTo(r, chunk, CurrentBase);
    }

    while (pos <= end)
    {
      ULong chunk = ParseChunk(num, pos, Base10_18_Zeroes);
      pos += (Int)Base10_18_Zeroes;
      ClassicMultiplication::MultiplyTo(r, Base10_18, CurrentBase);
      AddTo(r, chunk, CurrentBase);
    }

    return r;
  }

  // Memoized recursive doubling. thread_local cache lives in this TU only — every
  // consumer that includes Parser.h shares the same cache per-thread (a real
  // benefit of the .cpp split over the prior header-only design).
  std::vector<DataT> Pow10(SizeT digits)
  {
    static thread_local std::unordered_map<SizeT, std::vector<DataT>> cache;

    auto it = cache.find(digits);
    if (it != cache.end())
      return it->second;

    std::vector<DataT> value;
    if (digits == 0)
    {
      value = std::vector<DataT>{1};
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
      // Pow10(d) = Pow10(d/2)². One Square (single FFT in NTT) instead of
      // Multiply(p, p) (two FFTs).
      std::vector<DataT> p = Pow10(digits / 2);
      value = Square(p, CurrentBase);
    }
    else
    {
      SizeT lo = digits / 2;
      SizeT hi = digits - lo;
      value = Multiply(Pow10(hi), Pow10(lo), CurrentBase);
    }

    return cache.emplace(digits, value).first->second;
  }

  std::vector<DataT> ParseUnsignedDivideConquer(char const *num, Int start, Int end)
  {
    SizeT len = (SizeT)(end - start + 1);
    if (len <= DecimalDcThreshold)
      return ParseUnsignedLinear(num, start, end);

    SizeT lowDigits = len / 2;
    Int split = end - (Int)lowDigits;

    std::vector<DataT> high = ParseUnsignedDivideConquer(num, start, split);
    std::vector<DataT> low = ParseUnsignedDivideConquer(num, split + 1, end);
    std::vector<DataT> scale = Pow10(lowDigits);
    std::vector<DataT> scaledHigh = Multiply(high, scale, CurrentBase);
    std::vector<DataT> result = Add(scaledHigh, low, CurrentBase);
    TrimZerosToOne(result);
    return result;
  }

  // Parse decimal digit range [start..end] into base-2^32 limbs.
  std::vector<DataT> ParseUnsigned(char const *num, Int start, Int end)
  {
    if (start > end)
      return std::vector<DataT>();
    return ParseUnsignedDivideConquer(num, start, end);
  }

  BigInteger Parse(char const *num, Int start, Int *char_processed)
  {
    if (num == nullptr || start < 0)
      return BigInteger();

    Int end;
    bool isNegative;
    Int processed = FindRange(num, start, end, isNegative);
    if (char_processed != nullptr)
      *char_processed = processed;

    if (start > end)
      return BigInteger();

    std::vector<DataT> bigInt = ParseUnsigned(num, start, end);
    return BigInteger(bigInt, isNegative);
  }

  BigInteger Parse(char const *num, Int *char_processed)
  {
    return Parse(num, 0, char_processed);
  }

  namespace
  {
    constexpr char DigitPairs[] =
        "00010203040506070809"
        "10111213141516171819"
        "20212223242526272829"
        "30313233343536373839"
        "40414243444546474849"
        "50515253545556575859"
        "60616263646566676869"
        "70717273747576777879"
        "80818283848586878889"
        "90919293949596979899";

    void AppendUnsignedDecimal(ULong v, std::string &out)
    {
      char buf[20];
      char *p = buf + sizeof(buf);
      while (v >= 100)
      {
        ULong q = v / 100;
        unsigned rem = (unsigned)(v - q * 100);
        p -= 2;
        p[0] = DigitPairs[2 * rem];
        p[1] = DigitPairs[2 * rem + 1];
        v = q;
      }
      if (v >= 10)
      {
        p -= 2;
        p[0] = DigitPairs[2 * v];
        p[1] = DigitPairs[2 * v + 1];
      }
      else
      {
        *--p = (char)('0' + v);
      }
      out.append(p, (buf + sizeof(buf)) - p);
    }

    void AppendPaddedUnsignedDecimal(ULong v, SizeT width, std::string &out)
    {
      char buf[20];
      SizeT pos = width;
      while (pos >= 2)
      {
        ULong q = v / 100;
        unsigned rem = (unsigned)(v - q * 100);
        pos -= 2;
        buf[pos] = DigitPairs[2 * rem];
        buf[pos + 1] = DigitPairs[2 * rem + 1];
        v = q;
      }
      if (pos == 1)
        buf[0] = (char)('0' + v);
      out.append(buf, width);
    }
  }

  // Linear divmod-10^19 formatter. Appends to `out`. padTo > 0 pads to exactly
  // padTo decimal digits with leading zeros (used by D&C lower halves).
  void ToStringLinearAppend(std::vector<DataT> r, SizeT padTo, std::string &out)
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

    std::vector<ULong> chunks;
    chunks.reserve(r.size() + 1);
    while (!(r.size() == 1 && r[0] == 0))
      chunks.push_back((ULong)ClassicDivision::DivModTo(r, Base10_19, CurrentBase));

    int topDigits = 0;
    {
      ULong v = chunks.back();
      do { ++topDigits; v /= 10; } while (v);
    }
    SizeT natural = (chunks.size() - 1) * Base10_19_Zeroes + (SizeT)topDigits;

    if (padTo > natural)
      out.append(padTo - natural, '0');

    auto it = chunks.rbegin();
    {
      ULong v = *it++;
      AppendUnsignedDecimal(v, out);
    }
    for (; it != chunks.rend(); ++it)
      AppendPaddedUnsignedDecimal(*it, Base10_19_Zeroes, out);
  }

  // Chain entry: 10^digits + its precomputed Newton-reciprocal divider.
  // Built top-down so each level splits its parent in half — balanced T(N) = 2T(N/2) + M(N).
  namespace
  {
    SizeT EstimateDecimalDigits(std::vector<DataT> const &r)
    {
#if BIGMATH_LIMB_64
      constexpr SizeT limbBits = 64;
      SizeT highBits = limbBits - (SizeT)__builtin_clzll((ULong)r.back());
#else
      constexpr SizeT limbBits = 32;
      SizeT highBits = limbBits - (SizeT)__builtin_clz((unsigned)r.back());
#endif
      long double bits = (long double)(r.size() - 1) * (long double)limbBits + (long double)highBits;
      return (SizeT)std::floor(bits * 0.3010299956639811952137388947244930267682L) + 1;
    }

    struct DecimalDcEntry
    {
      SizeT digits;
      std::vector<DataT> value;
      std::shared_ptr<NewtonDivision::Divider> divider;
    };

    std::vector<DecimalDcEntry> BuildDecimalDcChain(SizeT topDigits)
    {
      std::vector<DecimalDcEntry> chain;
      for (SizeT d = topDigits; d >= ToStringDcThreshold / 2; d /= 2)
      {
        DecimalDcEntry e;
        e.digits = d;
        e.value = Pow10(d);
        e.divider = std::make_shared<NewtonDivision::Divider>(e.value, CurrentBase);
        chain.push_back(std::move(e));
      }
      return chain;
    }

    std::vector<DecimalDcEntry> const &GetDecimalDcChain(SizeT topDigits)
    {
      static thread_local std::unordered_map<SizeT, std::vector<DecimalDcEntry>> cache;
      auto it = cache.find(topDigits);
      if (it != cache.end())
        return it->second;

      auto inserted = cache.emplace(topDigits, BuildDecimalDcChain(topDigits));
      return inserted.first->second;
    }

    void ToStringDivConquer(
        std::vector<DataT> n,
        std::vector<DecimalDcEntry> const &chain,
        SizeT level,
        SizeT padTo,
        std::string &out)
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

      // n smaller than this level's divisor → top half empty; descend.
      if (Compare(n, chain[level].value) < 0)
      {
        ToStringDivConquer(std::move(n), chain, level + 1, padTo, out);
        return;
      }

      std::vector<DataT> q;
      std::vector<DataT> r;
      chain[level].divider->DivideAndRemainderInto(n, q, r);
      SizeT half = chain[level].digits;
      SizeT topPad = (padTo > half) ? padTo - half : 0;
      ToStringDivConquer(std::move(q), chain, level + 1, topPad, out);
      ToStringDivConquer(std::move(r), chain, level + 1, half, out);
    }
  } // namespace

  std::string ToString(std::vector<DataT> const &bigInt, SizeT start, SizeT end, bool isNeg)
  {
    if (IsZero(bigInt, start, end))
      return std::string("0");

    std::vector<DataT> r(bigInt.begin() + start, bigInt.begin() + end + 1);
    while (r.size() > 1 && r.back() == 0)
      r.pop_back();

    SizeT approxDigits = EstimateDecimalDigits(r);

    std::string s;
    s.reserve(approxDigits + (isNeg ? 2 : 1));
    if (isNeg)
      s.push_back('-');

    if (approxDigits < ToStringDcThreshold)
    {
      ToStringLinearAppend(std::move(r), 0, s);
      return s;
    }

    auto const &chain = GetDecimalDcChain(approxDigits / 2);
    ToStringDivConquer(std::move(r), chain, 0, 0, s);
    return s;
  }

  std::string ToString(std::vector<DataT> const &bigInt, bool isNeg)
  {
    return ToString(bigInt, 0, bigInt.size() - 1, isNeg);
  }

  std::string ToString(BigInteger const &bigInt)
  {
    return ToString(bigInt.GetInteger(), bigInt.IsNegative());
  }
}
