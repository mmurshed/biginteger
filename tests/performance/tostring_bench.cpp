// Focused ToString benchmark for decimal conversion experiments.
// Reports cold (first call at a size — includes the divider-chain build)
// and warm (best-of-3 with the chain cached) separately; single-iteration
// chain rows vary ±20-60% with process state, so the split is mandatory
// (PR #112 lesson, tostring_chain_plan.md methodology).

#include "biginteger/BigInteger.h"
#include "biginteger/common/Parser.h"

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

static std::string GenerateDigits(int digits, std::uint64_t seed)
{
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<int> dist(0, 9);

  std::string s;
  s.reserve((std::size_t)digits);
  s.push_back((char)('1' + (rng() % 9)));
  for (int i = 1; i < digits; ++i)
    s.push_back((char)('0' + dist(rng)));
  return s;
}

template <class F>
static double OnceMs(F &&fn)
{
  auto start = std::chrono::steady_clock::now();
  fn();
  auto end = std::chrono::steady_clock::now();
  return std::chrono::duration<double, std::milli>(end - start).count();
}

template <class F>
static double BestMs(F &&fn, int iters)
{
  double best = std::numeric_limits<double>::max();
  for (int i = 0; i < iters; ++i)
    best = std::min(best, OnceMs(fn));
  return best;
}

static int Iterations(int digits)
{
  if (digits <= 10000) return 20;
  if (digits <= 100000) return 7;
  if (digits <= 1000000) return 5;
  return 3;
}

int main(int argc, char **argv)
{
  std::vector<int> sizes;
  for (int i = 1; i < argc; ++i)
    sizes.push_back(std::atoi(argv[i]));
  if (sizes.empty())
    sizes = {100000, 500000, 1000000, 2000000};

  std::cout << "digits,tostring_cold_ms,tostring_warm_ms,parse_cold_ms,parse_warm_ms\n";
  for (int digits : sizes)
  {
    std::string input = GenerateDigits(digits, 0xB16B00B5ULL ^ (std::uint64_t)digits);

    BigInteger value;
    double parseCold = OnceMs([&]() { value = Parse(input.c_str()); });

    std::string out;
    double toStringCold = OnceMs([&]() { out = ToString(value); });
    if (out != input)
    {
      std::cerr << "round trip failed at " << digits << " digits\n";
      return 1;
    }

    int iters = Iterations(digits);
    double parseWarm = BestMs([&]() {
      BigInteger parsed = Parse(input.c_str());
      if (parsed.size() == 0) std::abort();
    }, iters);
    double toStringWarm = BestMs([&]() {
      std::string s = ToString(value);
      if (s.size() != input.size()) std::abort();
    }, iters);

    std::cout << digits << ','
              << std::fixed << std::setprecision(3)
              << toStringCold << ',' << toStringWarm << ','
              << parseCold << ',' << parseWarm << '\n'
              << std::flush;
  }
  return 0;
}
