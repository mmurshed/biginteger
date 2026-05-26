#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/Division.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/division/BurnikelZieglerDivision.h"
#include "biginteger/algorithms/division/FastDivision.h"
#include "biginteger/algorithms/division/NewtonDivision.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
#include "biginteger/algorithms/multiplication/NTTMultiplication.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/common/Util.h"

using namespace BigMath;
using namespace std;

namespace
{
  volatile SizeT sink = 0;

  vector<DataT> RandomNumber(SizeT limbs, mt19937_64 &gen)
  {
    uniform_int_distribution<uint64_t> digit(0, 0xFFFFFFFFULL);
    vector<DataT> value(limbs);
    for (SizeT i = 0; i < limbs; ++i)
      value[i] = (DataT)digit(gen);
    if (!value.empty())
      value.back() = (DataT)(0x80000000ULL | (digit(gen) & 0x7FFFFFFFULL));
    TrimZeros(value);
    return value.empty() ? vector<DataT>{0} : value;
  }

  template <typename Fn>
  double BestMs(Fn &&fn, int reps)
  {
    double best = numeric_limits<double>::max();
    for (int i = 0; i < reps; ++i)
    {
      auto start = chrono::high_resolution_clock::now();
      SizeT checksum = fn();
      auto end = chrono::high_resolution_clock::now();
      sink += checksum;
      best = min(best, chrono::duration<double, milli>(end - start).count());
    }
    return best;
  }

  int RepsForLimbs(SizeT limbs)
  {
    if (limbs <= 128)
      return 7;
    if (limbs <= 4096)
      return 5;
    return 3;
  }

  string Winner(initializer_list<pair<string, double>> values)
  {
    auto best = values.begin();
    for (auto it = values.begin(); it != values.end(); ++it)
      if (it->second < best->second)
        best = it;
    return best->first;
  }

  void PrintHeader(string const &title)
  {
    cout << "\n== " << title << " ==\n";
  }

  void TuneMultiplication(bool full)
  {
    PrintHeader("Multiplication Dispatch");
    cout << "limbs_each,total_limbs,classic_ms,karatsuba_ms,ntt_ms,winner\n";

    vector<SizeT> sizes = {
        8, 16, 24, 32, 40, 48, 64, 96, 128, 192, 256, 384, 512,
        1024, 2048, 4096, 8192, 12000, 16000};
    if (full)
    {
      sizes.push_back(20000);
      sizes.push_back(24576);
      sizes.push_back(32768);
    }

    mt19937_64 gen(0xB16B00B5);
    SizeT suggestedClassicTotal = CLASSIC_MULTIPLICATION_THRESHOLD;
    SizeT suggestedNttTotal = 0;
    SizeT previousTotal = 0;
    bool foundClassicCrossover = false;

    for (SizeT limbs : sizes)
    {
      vector<DataT> a = RandomNumber(limbs, gen);
      vector<DataT> b = RandomNumber(limbs, gen);
      int reps = RepsForLimbs(limbs);

      double classicMs = numeric_limits<double>::quiet_NaN();
      if (limbs <= 256)
      {
        classicMs = BestMs([&]() {
          auto r = ClassicMultiplication::Multiply(a, b, BigInteger::Base());
          return (SizeT)r.size();
        }, reps);
      }

      double karaMs = BestMs([&]() {
        auto r = KaratsubaMultiplication::Multiply(a, b, BigInteger::Base());
        return (SizeT)r.size();
      }, reps);

      double nttMs = BestMs([&]() {
        auto r = NTTMultiplication::Multiply(a, b, BigInteger::Base());
        return (SizeT)r.size();
      }, reps);

      string winner = limbs <= 256
                          ? Winner({{"classic", classicMs}, {"karatsuba", karaMs}, {"ntt", nttMs}})
                          : Winner({{"karatsuba", karaMs}, {"ntt", nttMs}});

      if (!foundClassicCrossover && limbs <= 256 && karaMs < classicMs)
      {
        suggestedClassicTotal = previousTotal;
        foundClassicCrossover = true;
      }
      if (suggestedNttTotal == 0 && nttMs < karaMs)
        suggestedNttTotal = limbs * 2;

      cout << limbs << ',' << limbs * 2 << ','
           << fixed << setprecision(4) << classicMs << ','
           << karaMs << ',' << nttMs << ','
           << winner << '\n';
      previousTotal = limbs * 2;
    }

    cout << "suggested CLASSIC_MULTIPLICATION_THRESHOLD ~= "
         << suggestedClassicTotal
         << '\n';
    cout << "suggested NTT_MULTIPLICATION_THRESHOLD ~= "
         << (suggestedNttTotal ? suggestedNttTotal : NTT_MULTIPLICATION_THRESHOLD)
         << '\n';
    cout << "compile override example: -DBIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD="
         << suggestedClassicTotal
         << " -DBIGMATH_NTT_MULTIPLICATION_THRESHOLD="
         << (suggestedNttTotal ? suggestedNttTotal : NTT_MULTIPLICATION_THRESHOLD)
         << '\n';

    PrintHeader("Multiplication Min-Limb Gate");
    cout << "small_limbs,large_limbs,classic_ms,karatsuba_ms,ntt_ms,winner\n";

    vector<SizeT> smallSizes = {4, 8, 16, 24, 32, 40, 48, 64};
    SizeT largeLimbs = full ? 4096 : 1024;
    SizeT suggestedMinLimb = 0;

    for (SizeT small : smallSizes)
    {
      vector<DataT> a = RandomNumber(largeLimbs, gen);
      vector<DataT> b = RandomNumber(small, gen);
      int reps = RepsForLimbs(small);

      double classicMs = BestMs([&]() {
        auto r = ClassicMultiplication::Multiply(a, b, BigInteger::Base());
        return (SizeT)r.size();
      }, reps);
      double karaMs = BestMs([&]() {
        auto r = KaratsubaMultiplication::Multiply(a, b, BigInteger::Base());
        return (SizeT)r.size();
      }, reps);
      double nttMs = BestMs([&]() {
        auto r = NTTMultiplication::Multiply(a, b, BigInteger::Base());
        return (SizeT)r.size();
      }, reps);

      string winner = Winner({{"classic", classicMs}, {"karatsuba", karaMs}, {"ntt", nttMs}});
      if (winner == "classic")
        suggestedMinLimb = small;

      cout << small << ',' << largeLimbs << ','
           << fixed << setprecision(4) << classicMs << ','
           << karaMs << ',' << nttMs << ','
           << winner << '\n';
    }

    cout << "suggested CLASSIC_MIN_LIMB_THRESHOLD ~= "
         << suggestedMinLimb
         << '\n';
    cout << "compile override example: -DBIGMATH_CLASSIC_MIN_LIMB_THRESHOLD="
         << suggestedMinLimb
         << '\n';
  }

  void TuneDivision(bool full)
  {
    PrintHeader("Division Dispatch");
    cout << "dividend_limbs,divisor_limbs,ratio,fast_ms,bz_ms,newton_ms,winner\n";

    vector<pair<SizeT, double>> cases = {
        {512, 2.0}, {768, 2.0}, {1024, 2.0},
        {2048, 1.5}, {2048, 2.0}, {2048, 3.0},
        {4096, 1.5}, {4096, 2.0}, {4096, 4.0},
        {8192, 1.25}, {8192, 1.5}, {8192, 2.0}, {8192, 4.0},
        {16384, 1.5}, {16384, 2.0},
        {24576, 1.5}};
    if (full)
    {
      cases.push_back({16384, 1.5});
      cases.push_back({24576, 2.0});
      cases.push_back({32768, 1.5});
    }

    mt19937_64 gen(0xD1715100);
    SizeT suggestedBzDivisor = 0;
    SizeT suggestedNewtonMedium = 0;
    SizeT suggestedNewtonLarge = 0;
    SizeT suggestedNewtonSkewNumerator = 0;
    SizeT suggestedNewtonSkewDenominator = 1;

    for (auto const &[divisorLimbs, ratio] : cases)
    {
      SizeT dividendLimbs = max((SizeT)(divisorLimbs + 1), (SizeT)(divisorLimbs * ratio));
      vector<DataT> a = RandomNumber(dividendLimbs, gen);
      vector<DataT> b = RandomNumber(divisorLimbs, gen);
      int reps = RepsForLimbs(divisorLimbs);

      double fastMs = BestMs([&]() {
        auto qr = FastDivision::DivideAndRemainder(a, b, BigInteger::Base(), true);
        return (SizeT)(qr.first.size() + qr.second.size());
      }, reps);

      double bzMs = BestMs([&]() {
        auto qr = BurnikelZieglerDivision::DivideAndRemainder(a, b, BigInteger::Base(), true);
        return (SizeT)(qr.first.size() + qr.second.size());
      }, reps);

      double newtonMs = numeric_limits<double>::quiet_NaN();
      if (divisorLimbs >= 8192)
      {
        newtonMs = BestMs([&]() {
          auto qr = NewtonDivision::DivideAndRemainder(a, b, BigInteger::Base(), true);
          return (SizeT)(qr.first.size() + qr.second.size());
        }, max(1, reps / 2));
      }

      string winner = divisorLimbs >= 8192
                          ? Winner({{"fast", fastMs}, {"bz", bzMs}, {"newton", newtonMs}})
                          : Winner({{"fast", fastMs}, {"bz", bzMs}});

      if (suggestedBzDivisor == 0 && bzMs < fastMs)
        suggestedBzDivisor = divisorLimbs;
      if (suggestedNewtonMedium == 0 && divisorLimbs >= 8192 && ratio >= 2.0 && newtonMs < min(fastMs, bzMs))
        suggestedNewtonMedium = divisorLimbs;
      if (suggestedNewtonLarge == 0 && divisorLimbs >= 8192 && ratio < 2.0 && newtonMs < min(fastMs, bzMs))
        suggestedNewtonLarge = divisorLimbs;
      if (suggestedNewtonSkewNumerator == 0 && divisorLimbs >= 8192 && newtonMs < min(fastMs, bzMs))
      {
        suggestedNewtonSkewNumerator = (SizeT)(ratio * 4.0 + 0.5);
        suggestedNewtonSkewDenominator = 4;
        while (suggestedNewtonSkewDenominator > 1 &&
               suggestedNewtonSkewNumerator % 2 == 0 &&
               suggestedNewtonSkewDenominator % 2 == 0)
        {
          suggestedNewtonSkewNumerator /= 2;
          suggestedNewtonSkewDenominator /= 2;
        }
      }

      cout << dividendLimbs << ',' << divisorLimbs << ','
           << fixed << setprecision(2) << ratio << ','
           << setprecision(4) << fastMs << ',' << bzMs << ',' << newtonMs << ','
           << winner << '\n';
    }

    cout << "suggested BZ divisor floor ~= "
         << (suggestedBzDivisor ? suggestedBzDivisor : BZ_DIVISOR_THRESHOLD)
         << '\n';
    cout << "suggested NEWTON_MEDIUM_B ~= "
         << (suggestedNewtonMedium ? suggestedNewtonMedium : NEWTON_MEDIUM_B)
         << '\n';
    cout << "suggested NEWTON_LARGE_B ~= "
         << (suggestedNewtonLarge ? suggestedNewtonLarge : NEWTON_LARGE_B)
         << '\n';
    cout << "suggested NEWTON skew ratio ~= "
         << (suggestedNewtonSkewNumerator ? suggestedNewtonSkewNumerator : NEWTON_SKEW_NUMERATOR)
         << '/'
         << (suggestedNewtonSkewNumerator ? suggestedNewtonSkewDenominator : NEWTON_SKEW_DENOMINATOR)
         << '\n';
    cout << "compile override example: -DBIGMATH_BZ_DIVISOR_THRESHOLD="
         << (suggestedBzDivisor ? suggestedBzDivisor : BZ_DIVISOR_THRESHOLD)
         << " -DBIGMATH_NEWTON_MEDIUM_B="
         << (suggestedNewtonMedium ? suggestedNewtonMedium : NEWTON_MEDIUM_B)
         << " -DBIGMATH_NEWTON_LARGE_B="
         << (suggestedNewtonLarge ? suggestedNewtonLarge : NEWTON_LARGE_B)
         << " -DBIGMATH_NEWTON_SKEW_NUMERATOR="
         << (suggestedNewtonSkewNumerator ? suggestedNewtonSkewNumerator : NEWTON_SKEW_NUMERATOR)
         << " -DBIGMATH_NEWTON_SKEW_DENOMINATOR="
         << (suggestedNewtonSkewNumerator ? suggestedNewtonSkewDenominator : NEWTON_SKEW_DENOMINATOR)
         << '\n';
  }
}

int main(int argc, char **argv)
{
  bool full = argc > 1 && string(argv[1]) == "--full";

  cout << "BigInteger dispatch tuner\n";
  cout << unitbuf;
  cout << "mode=" << (full ? "full" : "quick") << '\n';
  cout << "base=" << BigInteger::Base() << '\n';

  TuneMultiplication(full);
  TuneDivision(full);

  cerr << "checksum=" << sink << '\n';
  return 0;
}
