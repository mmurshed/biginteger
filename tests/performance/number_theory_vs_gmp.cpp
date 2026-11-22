#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <string_view>
#include <vector>

#include <gmp.h>

#include "biginteger/ops/Operations.h"

namespace
{
  using namespace BigMath;
  using Clock = std::chrono::steady_clock;

  struct Mpz
  {
    mpz_t value;

    Mpz() { mpz_init(value); }
    ~Mpz() { mpz_clear(value); }
    Mpz(const Mpz &) = delete;
    Mpz &operator=(const Mpz &) = delete;
  };

  void Import(Mpz &target, const BigInteger &source)
  {
    const auto bytes = source.ToByteArray();
    mpz_import(target.value, bytes.size(), 1, 1, 0, 0, bytes.data());
    if (source.IsNegative())
      mpz_neg(target.value, target.value);
  }

  BigInteger RandomInteger(unsigned bits, std::mt19937_64 &rng, bool negative = false)
  {
    std::vector<uint8_t> bytes((bits + 7) / 8);
    for (auto &byte : bytes)
      byte = static_cast<uint8_t>(rng());
    const unsigned highBit = (bits - 1) % 8;
    bytes.front() = (bytes.front() & ((1u << highBit) - 1)) | (1u << highBit);
    return BigInteger::FromByteArray(bytes, negative);
  }

  template <class BigMathCall, class GmpCall>
  void Measure(const char *operation, const char *shape, unsigned bits, unsigned cases,
               unsigned calls, BigMathCall runBigMath, GmpCall runGmp)
  {
    constexpr size_t rounds = 7;
    std::array<double, rounds> bigmath{}, gmp{};
    const auto time = [calls](auto &run) {
      const auto start = Clock::now();
      for (unsigned i = 0; i < calls; ++i)
        run();
      return std::chrono::duration<double, std::micro>(Clock::now() - start).count() / calls;
    };
    for (size_t i = 0; i < rounds; ++i)
    {
      if (i % 2 == 0)
      {
        bigmath[i] = time(runBigMath);
        gmp[i] = time(runGmp);
      }
      else
      {
        gmp[i] = time(runGmp);
        bigmath[i] = time(runBigMath);
      }
    }
    std::sort(bigmath.begin(), bigmath.end());
    std::sort(gmp.begin(), gmp.end());
    std::printf("%s,%s,%u,%u,%u,%.6f,%.6f\n",
                operation, shape, bits, cases, calls, bigmath[rounds / 2], gmp[rounds / 2]);
  }

  void BenchmarkPowMod(unsigned bits, bool odd, std::mt19937_64 &rng)
  {
    BigInteger base, exponent, modulus, result;
    Mpz gBase, gExponent, gModulus, expected, actual;
    const auto runBigMath = [&] { result = PowMod(base, exponent, modulus); };
    const auto runGmp = [&] { mpz_powm(expected.value, gBase.value, gExponent.value, gModulus.value); };
    constexpr unsigned cases = 24;
    for (unsigned i = 0; i < cases; ++i)
    {
      base = RandomInteger(bits + (i % 3) * 9, rng, i % 2 == 0);
      exponent = RandomInteger(i == cases - 1 ? bits : 1 + (i * 131) % bits, rng);
      auto bytes = RandomInteger(bits, rng).ToByteArray();
      bytes.back() = (bytes.back() & 0xfe) | static_cast<unsigned>(odd);
      modulus = BigInteger::FromByteArray(bytes, false);
      Import(gBase, base);
      Import(gExponent, exponent);
      Import(gModulus, modulus);
      runBigMath();
      runGmp();
      Import(actual, result);
      if (mpz_cmp(actual.value, expected.value) != 0)
        throw std::runtime_error("PowMod result differs from GMP");
    }
    Measure("powmod", odd ? "odd" : "even", bits, cases, std::max(4u, 16384u / bits),
            runBigMath, runGmp);
  }

  void BenchmarkGcd(unsigned bits, std::mt19937_64 &rng)
  {
    BigInteger a, b, divisor;
    ExtendedGcdResult extended;
    Mpz gA, gB, expected, x, y, actual;
    const auto runBigMath = [&] { divisor = Gcd(a, b); };
    const auto runGmp = [&] { mpz_gcd(expected.value, gA.value, gB.value); };
    const auto runExtended = [&] { extended = ExtendedGcd(a, b); };
    const auto runGmpExtended = [&] { mpz_gcdext(expected.value, x.value, y.value, gA.value, gB.value); };
    constexpr unsigned cases = 24;
    for (unsigned i = 0; i < cases; ++i)
    {
      const bool plantFactor = i % 2 == 0;
      const unsigned factorBits = bits / 4;
      const unsigned inputBits = plantFactor ? bits - factorBits : bits;
      a = RandomInteger(inputBits, rng, i % 4 < 2);
      b = RandomInteger(inputBits, rng, i % 3 == 0);
      if (plantFactor)
      {
        const auto factor = RandomInteger(factorBits, rng);
        a = a * factor;
        b = b * factor;
      }
      Import(gA, a);
      Import(gB, b);
      runBigMath();
      runGmp();
      Import(actual, divisor);
      if (mpz_cmp(actual.value, expected.value) != 0)
        throw std::runtime_error("Gcd result differs from GMP");
      runExtended();
      if (extended.gcd != divisor || a * extended.x + b * extended.y != divisor)
        throw std::runtime_error("ExtendedGcd returned an invalid divisor or coefficients");
    }
    const unsigned calls = std::max(4u, 16384u / bits);
    Measure("gcd", "balanced", bits, cases, calls, runBigMath, runGmp);
    Measure("extended_gcd", "balanced", bits, cases, calls, runExtended, runGmpExtended);
  }
}

int main(int argc, char **argv)
{
  const std::string_view operation = argc == 2 ? argv[1] : "all";
  if (argc > 2 || (operation != "all" && operation != "powmod" && operation != "gcd"))
  {
    std::fputs("Usage: number_theory_vs_gmp [powmod|gcd]\n", stderr);
    return 1;
  }

  std::mt19937_64 rng(0xB16B00B5);
  std::puts("operation,shape,bits,reference_checks,calls_per_round,bigmath_us,gmp_us");
  if (operation == "all" || operation == "powmod")
    for (unsigned bits : {256, 512, 1024, 2048, 4096})
      for (bool odd : {true, false})
        BenchmarkPowMod(bits, odd, rng);
  rng.seed(0x6CD);
  if (operation == "all" || operation == "gcd")
    for (unsigned bits : {256, 512, 1024, 2048, 4096})
      BenchmarkGcd(bits, rng);
}
