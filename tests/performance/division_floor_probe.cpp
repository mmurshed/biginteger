// Focused floor probe for division dispatch tuning.
// Usage: division_floor_probe <divisor:ratio> [...]
// e.g. division_floor_probe 1024:2.5 1536:2.5 2560:1.6
// Equal interleaved reps for fast/bz/newton/dispatch (the shape bench gives
// Newton reps/2, which overstates it on small shapes).

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
      value[i] = (DataT)gen();
    if (!value.empty())
      value.back() |= 0x8000000000000000ULL;
    TrimZerosToOne(value);
    return value;
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
}

int main(int argc, char **argv)
{
  std::vector<std::pair<SizeT, double>> cases;
  for (int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];
    auto colon = arg.find(':');
    if (colon == std::string::npos)
      continue;
    cases.push_back({(SizeT)std::strtoull(arg.substr(0, colon).c_str(), nullptr, 10),
                     std::strtod(arg.substr(colon + 1).c_str(), nullptr)});
  }
  if (cases.empty())
  {
    std::cerr << "usage: division_floor_probe <divisor:ratio> [...]\n";
    return 1;
  }

  const int reps = 7;
  std::cout << "divisor_limbs,ratio,dividend_limbs,fast_ms,bz_ms,newton_ms,dispatch_ms,winner\n";

  std::mt19937_64 gen(0xD1715100B16B00B5ULL);
  for (auto const &[n, ratio] : cases)
  {
    SizeT na = std::max((SizeT)(n + 1), (SizeT)((double)n * ratio));
    std::vector<DataT> a = RandomNumber(na, gen);
    std::vector<DataT> b = RandomNumber(n, gen);
    if (Compare(a, b) <= 0)
      a.push_back(1);

    double fast = std::numeric_limits<double>::max();
    double bz = fast, newton = fast, dispatch = fast;

    // Warm-up + correctness once per algorithm.
    CheckResult(a, b, FastDivision::DivideAndRemainder(a, b, BigInteger::Base(), true));
    CheckResult(a, b, BurnikelZieglerDivision::DivideAndRemainder(a, b, BigInteger::Base(), true));
    CheckResult(a, b, NewtonDivision::DivideAndRemainder(a, b, BigInteger::Base(), true));
    CheckResult(a, b, DivideAndRemainder(a, b, BigInteger::Base(), true));

    // Interleaved equal reps: one round = one timed run of each algorithm.
    for (int r = 0; r < reps; ++r)
    {
      auto time1 = [&](auto &&fn, double &best) {
        auto start = std::chrono::steady_clock::now();
        auto qr = fn();
        auto end = std::chrono::steady_clock::now();
        sink += (SizeT)(qr.first.size() + qr.second.size());
        best = std::min(best, std::chrono::duration<double, std::milli>(end - start).count());
      };
      time1([&] { return FastDivision::DivideAndRemainder(a, b, BigInteger::Base(), true); }, fast);
      time1([&] { return BurnikelZieglerDivision::DivideAndRemainder(a, b, BigInteger::Base(), true); }, bz);
      time1([&] { return NewtonDivision::DivideAndRemainder(a, b, BigInteger::Base(), true); }, newton);
      time1([&] { return DivideAndRemainder(a, b, BigInteger::Base(), true); }, dispatch);
    }

    const char *winner = "fast";
    double best = fast;
    if (bz < best) { best = bz; winner = "bz"; }
    if (newton < best) { best = newton; winner = "newton"; }
    if (dispatch < best) winner = "dispatch";

    std::cout << n << ',' << std::fixed << std::setprecision(2) << ratio << ',' << na << ','
              << std::setprecision(4)
              << fast << ',' << bz << ',' << newton << ',' << dispatch << ','
              << winner << '\n';
  }

  std::cerr << "checksum=" << sink << '\n';
  return 0;
}
