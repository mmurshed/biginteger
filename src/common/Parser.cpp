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
#include "biginteger/algorithms/Division.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/Squaring.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/algorithms/division/ClassicDivision.h"
#include "biginteger/algorithms/division/NewtonDivision.h"

#include "biginteger/common/DecimalConversionScope.h"
#include "biginteger/common/Parallel.h"

#include <cstring>
#include <memory>
#include <bit>
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

  static std::vector<DataT> CombineDecimalHalves(std::vector<DataT> const &high,
                                                 std::vector<DataT> const &low,
                                                 SizeT lowDigits)
  {
    detail::DecimalConversionScope conversion;
    std::vector<DataT> scale = Pow10(lowDigits);
    std::vector<DataT> scaledHigh = Multiply(high, scale, CurrentBase);
    std::vector<DataT> result = Add(scaledHigh, low, CurrentBase);
    TrimZerosToOne(result);
    return result;
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
    return CombineDecimalHalves(high, low, lowDigits);
  }

// Parallel parse fan-out: depth of the deterministic split replay before the
// one ParallelDo dispatch (2^splits leaf ranges), and the digit count below
// which the fan-out isn't worth the dispatch.
#ifndef BIGMATH_PARSE_PARALLEL_SPLITS
#define BIGMATH_PARSE_PARALLEL_SPLITS 3
#endif
#ifndef BIGMATH_PARSE_PARALLEL_THRESHOLD
#define BIGMATH_PARSE_PARALLEL_THRESHOLD 100000
#endif

#if BIGMATH_USE_THREADS
  namespace
  {
    // The split points are a pure function of (start, end), so the collect
    // and combine phases replay the exact ParseUnsignedDivideConquer tree:
    // leaves parse in parallel, the combines run on the caller afterwards
    // (their Multiply uses the threaded NTT; workers' nested ParallelDo runs
    // inline — Parallel.cpp tl_chunkDepth guard).
    void CollectParseRanges(Int start, Int end, int splits,
                            std::vector<std::pair<Int, Int>> &ranges)
    {
      SizeT len = (SizeT)(end - start + 1);
      if (splits == 0 || len <= DecimalDcThreshold)
      {
        ranges.push_back({start, end});
        return;
      }
      SizeT lowDigits = len / 2;
      Int split = end - (Int)lowDigits;
      CollectParseRanges(start, split, splits - 1, ranges);
      CollectParseRanges(split + 1, end, splits - 1, ranges);
    }

    std::vector<DataT> CombineParsed(Int start, Int end, int splits,
                                     std::vector<std::vector<DataT>> &leaves,
                                     SizeT &idx)
    {
      SizeT len = (SizeT)(end - start + 1);
      if (splits == 0 || len <= DecimalDcThreshold)
        return std::move(leaves[idx++]);

      SizeT lowDigits = len / 2;
      Int split = end - (Int)lowDigits;
      std::vector<DataT> high = CombineParsed(start, split, splits - 1, leaves, idx);
      std::vector<DataT> low = CombineParsed(split + 1, end, splits - 1, leaves, idx);
      return CombineDecimalHalves(high, low, lowDigits);
    }
  }
#endif

  // Parse decimal digit range [start..end] into base-2^32 limbs.
  std::vector<DataT> ParseUnsigned(char const *num, Int start, Int end)
  {
    if (start > end)
      return std::vector<DataT>();

#if BIGMATH_USE_THREADS
    SizeT len = (SizeT)(end - start + 1);
    if (len >= (SizeT)BIGMATH_PARSE_PARALLEL_THRESHOLD && ParallelNumThreads() > 1)
    {
      std::vector<std::pair<Int, Int>> ranges;
      ranges.reserve((SizeT)1 << BIGMATH_PARSE_PARALLEL_SPLITS);
      CollectParseRanges(start, end, BIGMATH_PARSE_PARALLEL_SPLITS, ranges);
      if (ranges.size() > 1)
      {
        std::vector<std::vector<DataT>> leaves(ranges.size());
        ParallelDo((Int)ranges.size(), [&](Int s, Int e) {
          for (Int i = s; i < e; ++i)
            leaves[(SizeT)i] =
                ParseUnsignedDivideConquer(num, ranges[(SizeT)i].first, ranges[(SizeT)i].second);
        });
        SizeT idx = 0;
        return CombineParsed(start, end, BIGMATH_PARSE_PARALLEL_SPLITS, leaves, idx);
      }
    }
#endif
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

// Parallel ToString fan-out (T1): depth of serial descent before the one
// ParallelDo dispatch (2^splits subtrees), and the digit count below which
// the fan-out isn't worth the extra copy + dispatch.
#ifndef BIGMATH_TOSTR_PARALLEL_SPLITS
#define BIGMATH_TOSTR_PARALLEL_SPLITS 3
#endif
#ifndef BIGMATH_TOSTR_PARALLEL_THRESHOLD
#define BIGMATH_TOSTR_PARALLEL_THRESHOLD 100000
#endif

// Smallest chain level (power-of-ten divisor, in decimal digits) that divides
// through its cached Newton reciprocal. Smaller levels go through the regular
// division dispatcher (Knuth D / Burnikel-Ziegler), which needs no reciprocal
// and costs about half as much there: dividing a (2d − 1)-digit value by 10^d
// on M1 Max, single-threaded, the dispatcher took 0.43-0.59x the cached-Newton
// time up to d = 6144 (319 limbs) and 0.75-0.79x at 11000-12000; the two meet
// near 13000 digits (675 limbs); above that Newton is ahead by 1.0-1.7x up to
// 22000 digits and 1.3-2.2x from 24576. With a warm thread pool Newton's
// products thread and the crossover halves (0.83-0.88x at 5500-6000 digits,
// 1.04x at 6500, 1.16-1.69x at 7500-11000).
// The threaded value applies to divisions on the calling thread; subtrees
// converted inside the parallel fan-out run serially on a worker and use the
// single-thread value. Reciprocals are built only for levels that use them.
#ifndef BIGMATH_TOSTR_NEWTON_MIN_DIGITS
#define BIGMATH_TOSTR_NEWTON_MIN_DIGITS 13000
#endif
#ifndef BIGMATH_TOSTR_NEWTON_MIN_DIGITS_THREADED
#define BIGMATH_TOSTR_NEWTON_MIN_DIGITS_THREADED 6500
#endif

  // Chain entry: 10^digits, plus its precomputed Newton-reciprocal divider
  // when the level is large enough to use one (see
  // BIGMATH_TOSTR_NEWTON_MIN_DIGITS). Built top-down so each level splits its
  // parent in half — balanced T(N) = 2T(N/2) + M(N).
  namespace
  {
    constexpr int ToStringParallelSplits = BIGMATH_TOSTR_PARALLEL_SPLITS;
    constexpr SizeT ToStringParallelThreshold = BIGMATH_TOSTR_PARALLEL_THRESHOLD;
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
      std::shared_ptr<NewtonDivision::Divider> divider; // null below the Newton minimum
    };

    static_assert(BIGMATH_TOSTR_NEWTON_MIN_DIGITS_THREADED <= BIGMATH_TOSTR_NEWTON_MIN_DIGITS,
                  "the threaded Newton minimum must not exceed the single-thread one");

    // Newton minimum for divisions made on the calling thread, for a chain
    // whose top level is `topDigits`. Below the threaded minimum no level can
    // use Newton in any build, so the thread count is not asked for: the
    // first ParallelNumThreads() call constructs the pool and starts its
    // workers, which would cost a small first conversion tens of µs.
    SizeT ToStringNewtonMinDigits(SizeT topDigits)
    {
      if (topDigits < (SizeT)BIGMATH_TOSTR_NEWTON_MIN_DIGITS_THREADED)
        return (SizeT)BIGMATH_TOSTR_NEWTON_MIN_DIGITS;
      return ParallelNumThreads() > 1 ? (SizeT)BIGMATH_TOSTR_NEWTON_MIN_DIGITS_THREADED
                                      : (SizeT)BIGMATH_TOSTR_NEWTON_MIN_DIGITS;
    }

    // q, r = n divmod chain[level].value: through the cached reciprocal when
    // the level has one and is at least newtonMinDigits, else the dispatcher.
    void SplitAtLevel(std::vector<DataT> const &n, DecimalDcEntry const &e, SizeT newtonMinDigits,
                      std::vector<DataT> &q, std::vector<DataT> &r)
    {
      if (e.divider && e.digits >= newtonMinDigits)
      {
        e.divider->DivideAndRemainderInto(n, q, r);
        return;
      }
      auto qr = DivideAndRemainder(n, e.value, CurrentBase);
      q = std::move(qr.first);
      r = std::move(qr.second);
    }

    std::vector<DecimalDcEntry> BuildDecimalDcChain(SizeT topDigits)
    {
      // The lowest minimum any caller uses (the parallel fan-out's workers
      // use the higher single-thread one).
      const SizeT buildMin = ToStringNewtonMinDigits(topDigits);
      std::vector<DecimalDcEntry> chain;
      for (SizeT d = topDigits; d >= ToStringDcThreshold / 2; d /= 2)
      {
        DecimalDcEntry e;
        e.digits = d;
        e.value = Pow10(d);
        // NOTE: building the dividers in a ParallelDo over levels was tried
        // (2026-06-12) and REGRESSES cold ToString 10-55%: the level-1/2
        // reciprocal builds lose their internal threaded NTT when run
        // serial-inline on a worker, which costs more than the fan-out wins.
        if (d >= buildMin)
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
        SizeT newtonMinDigits,
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
        ToStringDivConquer(std::move(n), chain, level + 1, padTo, newtonMinDigits, out);
        return;
      }

      std::vector<DataT> q;
      std::vector<DataT> r;
      SplitAtLevel(n, chain[level], newtonMinDigits, q, r);
      SizeT half = chain[level].digits;
      SizeT topPad = (padTo > half) ? padTo - half : 0;
      ToStringDivConquer(std::move(q), chain, level + 1, topPad, newtonMinDigits, out);
      ToStringDivConquer(std::move(r), chain, level + 1, half, newtonMinDigits, out);
    }

#if BIGMATH_USE_THREADS
    // T1 (tostring_chain_plan.md): the D&C recursion is embarrassingly
    // parallel below the first few divmods — every subtree's output is a
    // fixed-width field (value < 10^padTo by construction), so the substrings
    // land at disjoint, precomputable offsets. Descend serially through the
    // top `splits` levels (those divmods are the big O(M(L)) ones and thread
    // internally via NTT), collect the subtrees, then convert them with ONE
    // ParallelDo dispatch — the pool is non-reentrant by design and nested
    // dispatches run inline (Parallel.cpp tl_chunkDepth guard).
    struct ToStringSubtree
    {
      std::vector<DataT> n;
      SizeT level;
      SizeT padTo;
    };

    void CollectToStringSubtrees(
        std::vector<DataT> n,
        std::vector<DecimalDcEntry> const &chain,
        SizeT level,
        SizeT padTo,
        int splits,
        std::vector<ToStringSubtree> &items)
    {
      if (splits == 0 || level >= chain.size() || IsZero(n))
      {
        items.push_back({std::move(n), level, padTo});
        return;
      }
      // n smaller than this level's divisor → top half empty; descend without
      // consuming a split (no fan-out happened).
      if (Compare(n, chain[level].value) < 0)
      {
        CollectToStringSubtrees(std::move(n), chain, level + 1, padTo, splits, items);
        return;
      }
      std::vector<DataT> q;
      std::vector<DataT> r;
      SplitAtLevel(n, chain[level], ToStringNewtonMinDigits(chain[0].digits), q, r);
      SizeT half = chain[level].digits;
      SizeT topPad = (padTo > half) ? padTo - half : 0;
      CollectToStringSubtrees(std::move(q), chain, level + 1, topPad, splits - 1, items);
      CollectToStringSubtrees(std::move(r), chain, level + 1, half, splits - 1, items);
    }

    void ToStringDivConquerParallel(
        std::vector<DataT> n,
        std::vector<DecimalDcEntry> const &chain,
        std::string &out)
    {
      std::vector<ToStringSubtree> items;
      items.reserve((SizeT)1 << ToStringParallelSplits);
      CollectToStringSubtrees(std::move(n), chain, 0, 0, ToStringParallelSplits, items);

      // Head subtree carries padTo = 0 (variable width — no leading zeros) so
      // its length is unknown until converted; do it serially. Everything
      // after it is fixed-width.
      ToStringDivConquer(std::move(items[0].n), chain, items[0].level, items[0].padTo,
                         ToStringNewtonMinDigits(chain[0].digits), out);

      SizeT m = (SizeT)items.size();
      if (m <= 1)
        return;

      std::vector<SizeT> offset(m);
      SizeT pos = out.size();
      for (SizeT i = 1; i < m; ++i)
      {
        offset[i] = pos;
        pos += items[i].padTo;
      }
      out.resize(pos);

      ParallelDo((Int)(m - 1), [&](Int start, Int end) {
        detail::DecimalConversionScope conversion;
        for (Int t = start; t < end; ++t)
        {
          SizeT i = (SizeT)t + 1;
          std::string part;
          part.reserve(items[i].padTo);
          // Serial on the worker (nested dispatches run inline), so the
          // single-thread Newton minimum applies.
          ToStringDivConquer(std::move(items[i].n), chain, items[i].level, items[i].padTo,
                             (SizeT)BIGMATH_TOSTR_NEWTON_MIN_DIGITS, part);
          if (part.size() != items[i].padTo)
            std::abort(); // fixed-width invariant broken — corrupt output otherwise
          std::memcpy(&out[offset[i]], part.data(), part.size());
        }
      });
    }
#endif
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

    // Round the chain top up to a 1/16-octave grid so nearby digit counts
    // share one cached chain. The cache key used to be the EXACT digit count:
    // workloads converting many values of slightly different lengths rebuilt
    // the whole divider chain (the dominant cold cost) on every call. The
    // rounded top makes the first split at worst ~56/44 instead of exactly
    // half — benign, unlike the 73/27 power-of-2-tower split documented in
    // STRING_CONVERSION.md.
    SizeT grid = std::max<SizeT>(SizeT(1) << (std::bit_width(approxDigits) - 5), 1);
    SizeT rounded = ((approxDigits + grid - 1) / grid) * grid;

    detail::DecimalConversionScope conversion;
    auto const &chain = GetDecimalDcChain(rounded / 2);
#if BIGMATH_USE_THREADS
    if (approxDigits >= ToStringParallelThreshold && ParallelNumThreads() > 1 &&
        chain.size() >= (SizeT)ToStringParallelSplits)
    {
      ToStringDivConquerParallel(std::move(r), chain, s);
      return s;
    }
#endif
    ToStringDivConquer(std::move(r), chain, 0, 0, ToStringNewtonMinDigits(chain[0].digits), s);
    return s;
  }

  std::string ToString(std::vector<DataT> const &bigInt, bool isNeg)
  {
    return ToString(bigInt, 0, bigInt.size() - 1, isNeg);
  }

  std::string ToString(BigInteger const &bigInt)
  {
    return ToString(bigInt.AsVector().get(), bigInt.IsNegative());
  }
}
