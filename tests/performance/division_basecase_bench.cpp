// Basecase division benchmark: direct FastDivision vs GMP across the scalar
// band that feeds both production FastDivision and BZ recursion basecases.
//
// CSV columns:
//   bz_basecase,divisor_limbs,ratio,dividend_limbs,pattern,algorithm,mode,reps,bigmath_ms,gmp_ms,ratio_bm_gmp,check
//
// Optional args:
//   --quick       reduced matrix for smoke runs
//   --reps N      fixed repetitions instead of the adaptive policy
//   --seed N      random seed
//
// bz_basecase is the BZ recursion floor this binary was built with:
// BIGMATH_DC_DIV_THRESHOLD for 64-bit limbs, BIGMATH_BZ_RECURSION_THRESHOLD
// for the legacy 32-in-64 layout. The division_basecase_bench_bz<N> CMake
// targets link a library built with the same value, so the bench and the
// library run one definition of the BZ code.
//
// Both sides allocate their outputs on every timed call: BigMath returns new
// vectors, and the GMP side initializes and clears fresh mpz outputs.

#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/division/BurnikelZieglerDivision.h"
#include "biginteger/algorithms/division/FastDivision.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/common/Util.h"

#include <gmp.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace BigMath;

namespace
{
  volatile SizeT sink = 0;

#if BIGMATH_LIMB_64
  constexpr SizeT kBzBasecase = BIGMATH_DC_DIV_THRESHOLD;
#else
  constexpr SizeT kBzBasecase = BIGMATH_BZ_RECURSION_THRESHOLD;
#endif

  enum class Pattern
  {
    Random,
    AllMax,
    TopEqual,
    SparseLow,
    Exact
  };

  const char *PatternName(Pattern p)
  {
    switch (p)
    {
      case Pattern::Random: return "random";
      case Pattern::AllMax: return "all_max";
      case Pattern::TopEqual: return "top_equal";
      case Pattern::SparseLow: return "sparse_low";
      case Pattern::Exact: return "exact";
    }
    return "?";
  }

  DataT LimbMax()
  {
#if BIGMATH_LIMB_64
    return (DataT)0xFFFFFFFFFFFFFFFFULL;
#else
    return (DataT)0xFFFFFFFFULL;
#endif
  }

  DataT TopBit()
  {
#if BIGMATH_LIMB_64
    return (DataT)0x8000000000000000ULL;
#else
    return (DataT)0x80000000ULL;
#endif
  }

  DataT RandomLimb(std::mt19937_64 &gen)
  {
#if BIGMATH_LIMB_64
    return (DataT)gen();
#else
    return (DataT)(gen() & 0xFFFFFFFFULL);
#endif
  }

  void Trim(std::vector<DataT> &v)
  {
    TrimZerosToOne(v);
  }

  std::vector<DataT> RandomNumber(SizeT limbs, std::mt19937_64 &gen)
  {
    std::vector<DataT> value(limbs);
    for (SizeT i = 0; i < limbs; ++i)
      value[i] = RandomLimb(gen);
    if (!value.empty())
      value.back() |= TopBit();
    Trim(value);
    return value;
  }

  std::vector<DataT> MaxNumber(SizeT limbs)
  {
    return std::vector<DataT>(limbs, LimbMax());
  }

  std::vector<DataT> SparseLowNumber(SizeT limbs, std::mt19937_64 &gen)
  {
    std::vector<DataT> value(limbs, 0);
    if (limbs == 0)
      return std::vector<DataT>{0};
    value[0] = (DataT)(RandomLimb(gen) | 1);
    if (limbs > 2)
      value[limbs / 2] = (DataT)(RandomLimb(gen) | 1);
    value.back() = (DataT)(TopBit() | (RandomLimb(gen) & (TopBit() - 1)));
    Trim(value);
    return value;
  }

  std::vector<DataT> LowShifted(std::vector<DataT> const &value, SizeT shift)
  {
    std::vector<DataT> out(shift, 0);
    out.insert(out.end(), value.begin(), value.end());
    Trim(out);
    return out;
  }

  std::pair<std::vector<DataT>, std::vector<DataT>> MakeCase(
      SizeT divisorLimbs,
      double ratio,
      Pattern pattern,
      std::mt19937_64 &gen)
  {
    SizeT dividendLimbs = std::max(divisorLimbs + 1, (SizeT)(divisorLimbs * ratio));
    std::vector<DataT> b;
    std::vector<DataT> a;

    switch (pattern)
    {
      case Pattern::Random:
        b = RandomNumber(divisorLimbs, gen);
        a = RandomNumber(dividendLimbs, gen);
        break;
      case Pattern::AllMax:
        b = MaxNumber(divisorLimbs);
        a = MaxNumber(dividendLimbs);
        break;
      case Pattern::TopEqual:
        b = RandomNumber(divisorLimbs, gen);
        a = RandomNumber(dividendLimbs, gen);
        if (dividendLimbs >= divisorLimbs + 1 && divisorLimbs >= 2)
        {
          a[dividendLimbs - 1] = b.back();
          a[dividendLimbs - 2] = b[divisorLimbs - 2];
        }
        break;
      case Pattern::SparseLow:
        b = SparseLowNumber(divisorLimbs, gen);
        a = SparseLowNumber(dividendLimbs, gen);
        break;
      case Pattern::Exact:
        b = RandomNumber(divisorLimbs, gen);
        a = LowShifted(b, dividendLimbs - divisorLimbs);
        break;
    }

    if (Compare(a, b) <= 0)
      a.push_back(1);
    Trim(a);
    Trim(b);
    return {a, b};
  }

  void ImportToMpz(mpz_t z, std::vector<DataT> const &v)
  {
#if BIGMATH_LIMB_64
    mpz_import(z, v.size(), -1, sizeof(uint64_t), 0, 0, v.data());
#else
    std::vector<uint32_t> words(v.size());
    for (SizeT i = 0; i < v.size(); ++i)
      words[i] = (uint32_t)v[i];
    mpz_import(z, words.size(), -1, sizeof(uint32_t), 0, 0, words.data());
#endif
  }

  bool MpzEqualsVector(mpz_t z, std::vector<DataT> const &v)
  {
    mpz_t tmp;
    mpz_init(tmp);
    ImportToMpz(tmp, v);
    bool ok = (mpz_cmp(z, tmp) == 0);
    mpz_clear(tmp);
    return ok;
  }

  template <class F>
  double BestMs(F &&fn, int reps)
  {
    double best = std::numeric_limits<double>::max();
    for (int i = 0; i < reps; ++i)
    {
      auto start = std::chrono::steady_clock::now();
      SizeT s = fn();
      auto end = std::chrono::steady_clock::now();
      sink += s;
      best = std::min(best, std::chrono::duration<double, std::milli>(end - start).count());
    }
    return best;
  }

  int RepsFor(SizeT divisorLimbs, double ratio)
  {
    if (divisorLimbs <= 32) return 15;
    if (divisorLimbs <= 128) return 9;
    if (divisorLimbs <= 512) return ratio >= 5.0 ? 5 : 7;
    return ratio >= 5.0 ? 3 : 5;
  }

  void RunRow(
      SizeT divisorLimbs,
      double ratio,
      Pattern pattern,
      bool useBz,
      bool quotientOnly,
      int fixedReps,
      std::mt19937_64 &gen)
  {
    auto [a, b] = MakeCase(divisorLimbs, ratio, pattern, gen);
    int reps = fixedReps > 0 ? fixedReps : RepsFor(divisorLimbs, ratio);

    mpz_t ga, gb, gq, gr;
    mpz_init(ga);
    mpz_init(gb);
    mpz_init(gq);
    mpz_init(gr);
    ImportToMpz(ga, a);
    ImportToMpz(gb, b);

    bool ok = true;
    if (quotientOnly)
    {
      auto q = useBz
          ? BurnikelZieglerDivision::Divide(a, b, BigInteger::Base())
          : FastDivision::Divide(a, b, BigInteger::Base());
      mpz_tdiv_q(gq, ga, gb);
      ok = MpzEqualsVector(gq, q);
    }
    else
    {
      auto qr = useBz
          ? BurnikelZieglerDivision::DivideAndRemainder(a, b, BigInteger::Base(), true)
          : FastDivision::DivideAndRemainder(a, b, BigInteger::Base(), true);
      mpz_tdiv_qr(gq, gr, ga, gb);
      ok = MpzEqualsVector(gq, qr.first) && MpzEqualsVector(gr, qr.second);
    }

    double bm = BestMs([&]() {
      if (quotientOnly)
      {
        auto q = useBz
            ? BurnikelZieglerDivision::Divide(a, b, BigInteger::Base())
            : FastDivision::Divide(a, b, BigInteger::Base());
        return (SizeT)q.size();
      }
      auto qr = useBz
          ? BurnikelZieglerDivision::DivideAndRemainder(a, b, BigInteger::Base(), true)
          : FastDivision::DivideAndRemainder(a, b, BigInteger::Base(), true);
      return (SizeT)(qr.first.size() + qr.second.size());
    }, reps);

    // Fresh outputs per call, like the BigMath side, so allocation is not
    // charged to one library only.
    double gp = BestMs([&]() {
      mpz_t q, r;
      mpz_init(q);
      SizeT size;
      if (quotientOnly)
      {
        mpz_tdiv_q(q, ga, gb);
        size = (SizeT)mpz_size(q);
      }
      else
      {
        mpz_init(r);
        mpz_tdiv_qr(q, r, ga, gb);
        size = (SizeT)mpz_size(q) + (SizeT)mpz_size(r);
        mpz_clear(r);
      }
      mpz_clear(q);
      return size;
    }, reps);

    std::cout << kBzBasecase << ','
              << divisorLimbs << ','
              << std::fixed << std::setprecision(2) << ratio << ','
              << a.size() << ','
              << PatternName(pattern) << ','
              << (useBz ? "bz" : "fast") << ','
              << (quotientOnly ? "q" : "qr") << ','
              << reps << ','
              << std::setprecision(6) << bm << ','
              << gp << ','
              << (gp > 0 ? bm / gp : 0.0) << ','
              << (ok ? "ok" : "FAIL") << '\n';

    mpz_clear(ga);
    mpz_clear(gb);
    mpz_clear(gq);
    mpz_clear(gr);

    if (!ok)
    {
      std::cerr << "correctness failure at divisor_limbs=" << divisorLimbs
                << " ratio=" << ratio << " pattern=" << PatternName(pattern)
                << " algorithm=" << (useBz ? "bz" : "fast")
                << " mode=" << (quotientOnly ? "q" : "qr") << '\n';
      std::abort();
    }
  }

  void PrintUsage(char const *argv0)
  {
    std::cerr << "usage: " << argv0 << " [--quick] [--reps N] [--seed N]\n";
  }
}

int main(int argc, char **argv)
{
  bool quick = false;
  int fixedReps = 0;
  uint64_t seed = 0xBA5ECA5E5EEDULL;

  for (int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];
    if (arg == "--quick")
    {
      quick = true;
    }
    else if (arg == "--reps" && i + 1 < argc)
    {
      fixedReps = std::atoi(argv[++i]);
    }
    else if (arg == "--seed" && i + 1 < argc)
    {
      seed = std::strtoull(argv[++i], nullptr, 0);
    }
    else
    {
      PrintUsage(argv[0]);
      return 1;
    }
  }

  std::vector<SizeT> divisors = quick
      ? std::vector<SizeT>{2, 8, 32, 128, 512, 1024}
      : std::vector<SizeT>{2, 3, 4, 8, 16, 32, 64, 96, 128, 192, 256, 384, 512, 768, 1024};
  std::vector<double> ratios = quick
      ? std::vector<double>{1.10, 1.50, 3.00, 10.00}
      : std::vector<double>{1.05, 1.10, 1.25, 1.50, 2.00, 3.00, 5.00, 10.00};
  std::vector<Pattern> patterns = quick
      ? std::vector<Pattern>{Pattern::Random, Pattern::TopEqual, Pattern::Exact}
      : std::vector<Pattern>{Pattern::Random, Pattern::AllMax, Pattern::TopEqual, Pattern::SparseLow, Pattern::Exact};

  std::mt19937_64 gen(seed);
  std::cout << "bz_basecase,divisor_limbs,ratio,dividend_limbs,pattern,algorithm,mode,reps,bigmath_ms,gmp_ms,ratio_bm_gmp,check\n";
  for (SizeT n : divisors)
  {
    for (double ratio : ratios)
    {
      for (Pattern pattern : patterns)
      {
        RunRow(n, ratio, pattern, false, true, fixedReps, gen);
        RunRow(n, ratio, pattern, false, false, fixedReps, gen);
        if (n >= 64)
        {
          RunRow(n, ratio, pattern, true, true, fixedReps, gen);
          RunRow(n, ratio, pattern, true, false, fixedReps, gen);
        }
      }
    }
  }

  std::cerr << "checksum=" << sink << '\n';
  return 0;
}
