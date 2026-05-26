// Shape-focused division benchmark for dispatcher tuning.

#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/Division.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/division/BurnikelZieglerDivision.h"
#include "biginteger/algorithms/division/FastDivision.h"
#include "biginteger/algorithms/division/NewtonDivision.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/common/Util.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

using namespace BigMath;

namespace
{
  volatile SizeT sink = 0;

  std::vector<DataT> RandomNumber(SizeT limbs, std::mt19937_64 &gen)
  {
    std::vector<DataT> value(limbs);
    for (SizeT i = 0; i < limbs; ++i)
    {
#if BIGMATH_LIMB_64
      value[i] = (DataT)gen();
#else
      value[i] = (DataT)(gen() & 0xFFFFFFFFULL);
#endif
    }
    if (!value.empty())
    {
#if BIGMATH_LIMB_64
      value.back() |= 0x8000000000000000ULL;
#else
      value.back() |= 0x80000000ULL;
#endif
    }
    TrimZerosToOne(value);
    return value;
  }

  template <class F>
  double BestMs(F &&fn, int reps)
  {
    double best = std::numeric_limits<double>::max();
    for (int i = 0; i < reps; ++i)
    {
      auto start = std::chrono::steady_clock::now();
      SizeT check = fn();
      auto end = std::chrono::steady_clock::now();
      sink += check;
      best = std::min(best, std::chrono::duration<double, std::milli>(end - start).count());
    }
    return best;
  }

  int RepsFor(SizeT divisorLimbs)
  {
    if (divisorLimbs <= 256) return 7;
    if (divisorLimbs <= 2048) return 5;
    if (divisorLimbs <= 8192) return 3;
    return 1;
  }

  std::string Winner(double fast, double bz, double newton, double dispatch)
  {
    std::string name = "fast";
    double best = fast;
    if (bz < best)
    {
      best = bz;
      name = "bz";
    }
    if (newton < best)
    {
      best = newton;
      name = "newton";
    }
    if (dispatch < best)
      name = "dispatch";
    return name;
  }

  void CheckResult(
      std::vector<DataT> const &a,
      std::vector<DataT> const &b,
      std::pair<std::vector<DataT>, std::vector<DataT>> const &qr)
  {
    std::vector<DataT> qb = Multiply(qr.first, b, BigInteger::Base());
    std::vector<DataT> recomposed = Add(qb, qr.second, BigInteger::Base());
    TrimZerosToOne(recomposed);
    if (Compare(recomposed, a) != 0 || Compare(qr.second, b) >= 0)
    {
      std::cerr << "division result mismatch\n";
      std::abort();
    }
  }

  struct Case
  {
    const char *shape;
    SizeT dividendLimbs;
    SizeT divisorLimbs;
  };
}

int main(int argc, char **argv)
{
  std::vector<SizeT> divisors = {512, 1024, 2048, 4096, 8192};
  if (argc > 1)
  {
    divisors.clear();
    for (int i = 1; i < argc; ++i)
      divisors.push_back((SizeT)std::strtoull(argv[i], nullptr, 10));
  }

  std::cout << "base=" << BigInteger::Base() << "\n";
  std::cout << "shape,dividend_limbs,divisor_limbs,fast_ms,bz_ms,newton_ms,dispatch_ms,winner\n";

  std::mt19937_64 gen(0xD1715100B16B00B5ULL);
  for (SizeT n : divisors)
  {
    std::vector<Case> cases = {
        {"2n/n", 2 * n, n},
        {"3n/2n", 3 * n, 2 * n},
        {"3n/n", 3 * n, n},
        {"5n/2n", 5 * n, 2 * n},
        {"5n/n", 5 * n, n},
        {"10n/n", 10 * n, n},
    };

    for (auto const &c : cases)
    {
      std::vector<DataT> a = RandomNumber(c.dividendLimbs, gen);
      std::vector<DataT> b = RandomNumber(c.divisorLimbs, gen);
      if (Compare(a, b) <= 0)
        a.push_back(1);

      auto reference = FastDivision::DivideAndRemainder(a, b, BigInteger::Base(), true);
      CheckResult(a, b, reference);

      int reps = RepsFor(c.divisorLimbs);

      double fast = BestMs([&]() {
        auto qr = FastDivision::DivideAndRemainder(a, b, BigInteger::Base(), true);
        return (SizeT)(qr.first.size() + qr.second.size());
      }, reps);

      double bz = BestMs([&]() {
        auto qr = BurnikelZieglerDivision::DivideAndRemainder(a, b, BigInteger::Base(), true);
        CheckResult(a, b, qr);
        return (SizeT)(qr.first.size() + qr.second.size());
      }, reps);

      double newton = BestMs([&]() {
        auto qr = NewtonDivision::DivideAndRemainder(a, b, BigInteger::Base(), true);
        CheckResult(a, b, qr);
        return (SizeT)(qr.first.size() + qr.second.size());
      }, std::max(1, reps / 2));

      double dispatch = BestMs([&]() {
        auto qr = DivideAndRemainder(a, b, BigInteger::Base(), true);
        CheckResult(a, b, qr);
        return (SizeT)(qr.first.size() + qr.second.size());
      }, reps);

      std::cout << c.shape << ','
                << c.dividendLimbs << ','
                << c.divisorLimbs << ','
                << std::fixed << std::setprecision(4)
                << fast << ',' << bz << ',' << newton << ',' << dispatch << ','
                << Winner(fast, bz, newton, dispatch) << '\n';
    }
  }

  std::cerr << "checksum=" << sink << '\n';
  return 0;
}
