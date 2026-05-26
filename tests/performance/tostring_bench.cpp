// Focused ToString benchmark for decimal conversion experiments.

#include "biginteger/BigInteger.h"
#include "biginteger/common/Parser.h"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
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
static double TimeMs(F &&fn, int iters)
{
  auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < iters; ++i)
    fn();
  auto end = std::chrono::steady_clock::now();
  return std::chrono::duration<double, std::milli>(end - start).count() / iters;
}

static int Iterations(int digits)
{
  if (digits <= 1000) return 80;
  if (digits <= 10000) return 20;
  if (digits <= 50000) return 6;
  if (digits <= 100000) return 3;
  return 1;
}

int main(int argc, char **argv)
{
  std::vector<int> sizes;
  if (argc > 1)
  {
    for (int i = 1; i < argc; ++i)
      sizes.push_back(std::atoi(argv[i]));
  }
  else
  {
    sizes = {1000, 10000, 50000, 100000, 200000};
  }

  std::cout << "digits,iters,parse_ms,tostring_ms\n";
  for (int digits : sizes)
  {
    std::string input = GenerateDigits(digits, 0xB16B00B5ULL ^ (std::uint64_t)digits);
    BigInteger value = Parse(input.c_str());
    std::string warm = ToString(value);
    if (warm != input)
    {
      std::cerr << "round trip failed at " << digits << " digits\n";
      return 1;
    }

    int iters = Iterations(digits);
    double parseMs = TimeMs([&]() {
      BigInteger parsed = Parse(input.c_str());
      if (parsed.size() == 0) std::abort();
    }, iters);

    double toStringMs = TimeMs([&]() {
      std::string out = ToString(value);
      if (out.size() != input.size()) std::abort();
    }, iters);

    std::cout << digits << ',' << iters << ','
              << std::fixed << std::setprecision(4)
              << parseMs << ',' << toStringMs << '\n';
  }
  return 0;
}
