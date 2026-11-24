#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <gmp.h>

#include "biginteger/ops/Operations.h"
#include "../../src/ops/Factorial.h"
#include "../../src/ops/PrimeSearch.h"

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

  BigInteger Export(Mpz const &input)
  {
    std::vector<uint8_t> bytes((mpz_sizeinbase(input.value, 2) + 7) / 8);
    size_t count;
    mpz_export(bytes.data(), &count, 1, 1, 0, 0, input.value);
    return BigInteger::FromByteArray(std::span(bytes).first(count), false);
  }

  template <class BigMathCall, class GmpCall>
  std::pair<double, double> MeasureTimes(unsigned calls, BigMathCall runBigMath, GmpCall runGmp)
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
    return {bigmath[rounds / 2], gmp[rounds / 2]};
  }

  template <class BigMathCall, class GmpCall>
  void Measure(const char *operation, const char *shape, unsigned bits, unsigned cases,
               unsigned calls, BigMathCall runBigMath, GmpCall runGmp)
  {
    const auto [bigmath, gmp] = MeasureTimes(calls, runBigMath, runGmp);
    std::printf("%s,%s,%u,%u,%u,%.6f,%.6f\n",
                operation, shape, bits, cases, calls, bigmath, gmp);
  }

  void BenchmarkFactorial()
  {
    std::puts("n,output_bits,calls_per_round,bigmath_us,gmp_us,bigmath_control_us,product_tree_us");
    for (UInt n : {0u, 20u, 21u, 100u, 128u, 256u, 1000u, 10000u, 100000u, 1000000u})
    {
      BigInteger result, productTree;
      Mpz expected;
      const auto runFactorial = [&] { result = Factorial(n); };
      const auto runGmp = [&] { mpz_fac_ui(expected.value, n); };
      const auto runTree = [&] { productTree = detail::FactorialProductTree(n); };
      runFactorial();
      runGmp();
      runTree();
      if (result != Export(expected) || result != productTree)
        throw std::runtime_error("factorial differs from reference");
      const unsigned calls = std::max(1u, 50000u / std::max(20u, n));
      const auto [bigmath, gmp] = MeasureTimes(calls, runFactorial, runGmp);
      const auto [control, tree] = MeasureTimes(calls, runFactorial, runTree);
      std::printf("%u,%zu,%u,%.6f,%.6f,%.6f,%.6f\n", n, mpz_sizeinbase(expected.value, 2),
                  calls, bigmath, gmp, control, tree);
    }
  }

  void BenchmarkPowMod(unsigned bits, bool odd, bool sparse, std::mt19937_64 &rng,
                       unsigned cases = 24, unsigned exponentBits = 0)
  {
    BigInteger base, exponent, modulus, result;
    Mpz gBase, gExponent, gModulus, expected, actual;
    const auto runBigMath = [&] { result = PowMod(base, exponent, modulus); };
    const auto runGmp = [&] { mpz_powm(expected.value, gBase.value, gExponent.value, gModulus.value); };
    const auto denseBits = exponentBits == 0 ? bits : exponentBits;
    for (unsigned i = 0; i < cases; ++i)
    {
      base = RandomInteger(bits + (i % 3) * 9, rng, i % 2 == 0);
      if (sparse)
      {
        const DataT value = 65537;
        exponent = BigInteger::FromLimbs({&value, 1}, false);
      }
      else
        exponent = RandomInteger(i == cases - 1 ? denseBits : 1 + (i * 131) % denseBits, rng);
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
    std::string shape = odd ? "odd" : "even";
    if (sparse)
      shape += "_65537";
    else if (exponentBits != 0)
      shape += "_" + std::to_string(exponentBits);
    const unsigned calls = sparse ? std::max(16u, 131072u / bits)
                                  : std::max(bits <= 4096 ? 4u : 1u, 16384u / bits);
    Measure("powmod", shape.c_str(), bits, cases, calls, runBigMath, runGmp);
  }

  void BenchmarkGcd(unsigned bits, std::string_view shape, std::mt19937_64 &rng,
                    unsigned cases = 24)
  {
    BigInteger a, b, divisor;
    ExtendedGcdResult extended;
    Mpz gA, gB, expected, x, y, actual;
    const auto runBigMath = [&] { divisor = Gcd(a, b); };
    const auto runGmp = [&] { mpz_gcd(expected.value, gA.value, gB.value); };
    const auto runExtended = [&] { extended = ExtendedGcd(a, b); };
    const auto runGmpExtended = [&] { mpz_gcdext(expected.value, x.value, y.value, gA.value, gB.value); };
    const bool fibonacci = shape == "fibonacci";
    const bool skewed = shape == "skewed";
    const bool largeFactor = shape == "common_factor";
    for (unsigned i = 0; i < cases; ++i)
    {
      if (fibonacci)
      {
        // Consecutive Fibonacci numbers of approximately the requested width.
        const auto index = static_cast<unsigned long>(bits * 1.44042009) + 3;
        mpz_fib2_ui(gA.value, gB.value, index);
        a = Export(gA);
        b = Export(gB);
      }
      else
      {
        const bool plantFactor = largeFactor || (!skewed && i % 2 == 0);
        const unsigned factorBits = largeFactor ? 3 * bits / 4 : bits / 4;
        const unsigned inputBits = plantFactor ? bits - factorBits : bits;
        a = RandomInteger(inputBits, rng, i % 4 < 2);
        b = RandomInteger(skewed ? bits / 8 : inputBits, rng, i % 3 == 0);
        if (plantFactor)
        {
          const auto factor = RandomInteger(factorBits, rng);
          a = a * factor;
          b = b * factor;
        }
        Import(gA, a);
        Import(gB, b);
      }
      runBigMath();
      runGmp();
      Import(actual, divisor);
      if (mpz_cmp(actual.value, expected.value) != 0)
        throw std::runtime_error("Gcd result differs from GMP");
      runExtended();
      Import(x, extended.x);
      Import(y, extended.y);
      mpz_mul(actual.value, gA.value, x.value);
      mpz_addmul(actual.value, gB.value, y.value);
      if (extended.gcd != divisor || mpz_cmp(actual.value, expected.value) != 0)
        throw std::runtime_error("ExtendedGcd returned an invalid divisor or coefficients");
    }
    const unsigned calls = std::max(4u, 16384u / bits);
    Measure("gcd", shape.data(), bits, cases, calls, runBigMath, runGmp);
    Measure("extended_gcd", shape.data(), bits, cases, calls, runExtended, runGmpExtended);
  }

  void BenchmarkSqrtRem(unsigned bits, std::mt19937_64 &rng)
  {
    BigInteger n;
    std::pair<BigInteger, BigInteger> result;
    Mpz input, root, remainder, actualRoot, actualRemainder;
    const auto runBigMath = [&] { result = SqrtRem(n); };
    const auto runGmp = [&] { mpz_sqrtrem(root.value, remainder.value, input.value); };
    const unsigned cases = bits <= 65536 ? 24 : 8;
    constexpr DataT unit = 1;
    const auto one = BigInteger::FromLimbs({&unit, 1}, false);
    for (unsigned i = 0; i < cases; ++i)
    {
      if (i % 4 == 3)
        n = RandomInteger(bits, rng);
      else
      {
        const auto floor = RandomInteger(bits / 2, rng);
        n = floor * floor;
        if (i % 4 == 0)
          n = n - one;
        if (i % 4 == 2)
          n = n + one;
      }
      Import(input, n);
      runBigMath();
      runGmp();
      Import(actualRoot, result.first);
      Import(actualRemainder, result.second);
      if (mpz_cmp(actualRoot.value, root.value) != 0 ||
          mpz_cmp(actualRemainder.value, remainder.value) != 0)
        throw std::runtime_error("SqrtRem root or remainder differs from GMP");
    }
    Measure("sqrtrem", "random", bits, cases, std::max(4u, 131072u / bits), runBigMath, runGmp);
  }

  void BenchmarkRootRem()
  {
    std::mt19937_64 rng(0x72E36);
    std::puts("input_bits,degree,reference_checks,calls_per_round,bigmath_us,gmp_us");
    for (unsigned bits : {64u, 1024u, 16384u, 262144u})
      for (UInt degree : {3u, 5u, 16u, 127u, 1024u})
      {
        BigInteger n;
        std::pair<BigInteger, BigInteger> result;
        Mpz input, root, remainder, actualRoot, actualRemainder, base, power;
        const auto runBigMath = [&] { result = RootRem(n, degree); };
        const auto runGmp = [&] { mpz_rootrem(root.value, remainder.value, input.value, degree); };
        const auto one = RandomInteger(1, rng);
        constexpr unsigned cases = 12;
        for (unsigned i = 0; i < cases; ++i)
        {
          if (i % 4 == 3)
            n = RandomInteger(bits, rng);
          else
          {
            Import(base, RandomInteger(std::max(2u, bits / degree), rng));
            mpz_pow_ui(power.value, base.value, degree);
            n = Export(power);
            if (i % 4 == 0)
              n = n - one;
            if (i % 4 == 2)
              n = n + one;
          }
          if ((degree & 1) && i % 3 == 0)
            n = -n;
          Import(input, n);
          runBigMath();
          runGmp();
          Import(actualRoot, result.first);
          Import(actualRemainder, result.second);
          if (mpz_cmp(actualRoot.value, root.value) != 0 ||
              mpz_cmp(actualRemainder.value, remainder.value) != 0)
            throw std::runtime_error("RootRem differs from GMP");
        }
        const unsigned calls = std::clamp(65536u / bits, 2u, 256u);
        const auto [bigmath, gmp] = MeasureTimes(calls, runBigMath, runGmp);
        std::printf("%u,%u,%u,%u,%.6f,%.6f\n", bits, degree, cases, calls, bigmath, gmp);
      }
  }

  void BenchmarkPrimality(unsigned bits, std::string_view shape, std::mt19937_64 &rng)
  {
    constexpr int bigMathRounds = 32, gmpReps = 56;
    const bool prime = shape == "prime";
    const unsigned cases = bits <= 512 ? 16 : (bits <= 1024 ? 8 : 3);
    BigInteger candidate;
    Mpz input, other;
    auto checkBigMath = [&] {
      std::atomic_signal_fence(std::memory_order_seq_cst);
      if (IsProbablePrime(candidate, bigMathRounds) != prime)
        throw std::runtime_error("unexpected BigMath primality result");
    };
    auto checkGmp = [&] {
      std::atomic_signal_fence(std::memory_order_seq_cst);
      if ((mpz_probab_prime_p(input.value, gmpReps) != 0) != prime)
        throw std::runtime_error("unexpected GMP primality result");
    };
    for (unsigned i = 0; i < cases; ++i)
    {
      if (shape == "semiprime")
      {
        // Start in the upper quarter of each half-width range so the
        // product of the two factors has the requested bit length.
        const auto generateFactor = [&](Mpz &target) {
          Import(target, RandomInteger(bits / 2, rng));
          mpz_setbit(target.value, bits / 2 - 2);
          mpz_nextprime(target.value, target.value);
        };
        generateFactor(input);
        generateFactor(other);
        mpz_mul(input.value, input.value, other.value);
      }
      else
      {
        Import(input, RandomInteger(bits, rng));
        if (prime)
          mpz_nextprime(input.value, input.value);
        else
        {
          mpz_setbit(input.value, 0);
          mpz_add_ui(input.value, input.value, 2 * mpz_fdiv_ui(input.value, 3));
        }
      }
      if (mpz_sizeinbase(input.value, 2) != bits)
        throw std::runtime_error("primality fixture crossed its intended bit length");
      candidate = Export(input);
      checkBigMath();
      checkGmp();
    }
    const unsigned calls = shape == "small_factor" ? 1024
        : std::max(1u, (prime && bits > 64 ? 1024u : 8192u) / bits);
    Measure("primality", shape.data(), bits, cases, calls, checkBigMath, checkGmp);
  }

  void BenchmarkPseudoprimes()
  {
    for (const char *text : {"3215031751", "341550071728321", "3825123056546413051",
                             "318665857834031151167461"})
    {
      Mpz input;
      if (mpz_set_str(input.value, text, 10) != 0)
        throw std::runtime_error("invalid pseudoprime fixture");
      const auto candidate = Export(input);
      const auto checkBigMath = [&] {
        std::atomic_signal_fence(std::memory_order_seq_cst);
        if (IsProbablePrime(candidate))
          throw std::runtime_error("BigMath accepted a known pseudoprime");
      };
      const auto checkGmp = [&] {
        std::atomic_signal_fence(std::memory_order_seq_cst);
        if (mpz_probab_prime_p(input.value, 56) != 0)
          throw std::runtime_error("GMP accepted a known pseudoprime");
      };
      checkBigMath();
      checkGmp();
      const auto bits = static_cast<unsigned>(mpz_sizeinbase(input.value, 2));
      Measure("primality", "strong_pseudoprime", bits, 1, 128, checkBigMath, checkGmp);
    }
  }

  void BenchmarkNextProbablePrime()
  {
    std::mt19937_64 rng(0x4E455854);
    std::puts("bits,shape,guaranteed_composite_run,gap,odd_candidates,miller_rabin_candidates,bigmath_ms,gmp_ms");
    for (unsigned bits : {256, 1024, 4096})
      for (const bool longGap : {false, true})
      {
        Mpz input, reference, delta;
        BigInteger result;
        Import(input, RandomInteger(bits, rng));
        unsigned guaranteed = 0;
        if (longGap)
        {
          // If m is divisible by every integer through k, each of
          // m+2, ..., m+k is composite. Start just before that run.
          const unsigned k = bits * 2 / 3;
          Mpz multiple;
          mpz_set_ui(multiple.value, 1);
          for (unsigned i = 2; i <= k; ++i)
            mpz_lcm_ui(multiple.value, multiple.value, i);
          mpz_fdiv_r(delta.value, input.value, multiple.value);
          mpz_sub(input.value, input.value, delta.value);
          mpz_add_ui(input.value, input.value, 1);
          guaranteed = k - 1;
        }
        if (mpz_sizeinbase(input.value, 2) != bits)
          throw std::runtime_error("next-prime fixture crossed its intended bit length");

        const auto start = Export(input);
        PrimeSearch::SearchStats stats;
        const auto countedResult = PrimeSearch::FindNext(start, 32, &stats);
        const auto runBigMath = [&] { result = NextProbablePrime(start); };
        const auto runGmp = [&] { mpz_nextprime(reference.value, input.value); };
        const auto [bigmath, gmp] = MeasureTimes(1, runBigMath, runGmp);
        if (result != Export(reference) || result != countedResult)
          throw std::runtime_error("next prime differs from GMP");
        mpz_sub(delta.value, reference.value, input.value);
        const auto gap = mpz_get_ui(delta.value);
        const auto candidates = (gap + (mpz_even_p(input.value) != 0)) / 2;
        if (gap <= guaranteed || stats.screenedCandidates > candidates + 7)
          throw std::runtime_error("invalid next-prime candidate counts");
        std::printf("%u,%s,%u,%lu,%lu,%zu,%.6f,%.6f\n", bits, longGap ? "long_gap" : "ordinary",
                    guaranteed, gap, candidates, stats.screenedCandidates, bigmath / 1000, gmp / 1000);
        std::fflush(stdout);
      }
  }
}

int main(int argc, char **argv)
{
  const std::string_view operation = argc == 2 ? argv[1] : "all";
  if (operation == "rootrem")
  {
    BenchmarkRootRem();
    return 0;
  }
  if (operation == "factorial")
  {
    BenchmarkFactorial();
    return 0;
  }
  if (operation == "nextprime")
  {
    BenchmarkNextProbablePrime();
    return 0;
  }
  if (argc > 2 || (operation != "all" && operation != "powmod" &&
                   operation != "powmod-large" && operation != "gcd" && operation != "gcd-large" &&
                   operation != "sqrtrem" && operation != "prime"))
  {
    std::fputs("Usage: number_theory_vs_gmp [powmod|powmod-large|gcd|gcd-large|sqrtrem|rootrem|prime|nextprime|factorial]\n", stderr);
    return 1;
  }

  std::mt19937_64 rng(0xB16B00B5);
  std::puts("operation,shape,bits,reference_checks,calls_per_round,bigmath_us,gmp_us");
  if (operation == "all" || operation == "powmod")
  {
    for (bool sparse : {false, true})
      for (unsigned bits : {32, 64, 128, 256, 512, 1024, 2048, 4096})
        for (bool odd : {true, false})
          BenchmarkPowMod(bits, odd, sparse, rng);
  }
  if (operation == "powmod-large")
  {
    for (unsigned bits : {3072, 4096, 4097, 6144, 8192, 8193, 12288, 16384})
      BenchmarkPowMod(bits, true, false, rng, 4);
    for (unsigned bits : {24576, 32768, 32769, 49152, 65536, 65537})
      BenchmarkPowMod(bits, true, false, rng, 4, 257);
    for (unsigned bits : {4097, 8192, 8193, 16384, 32768, 65536})
      BenchmarkPowMod(bits, true, true, rng, 4);
  }
  rng.seed(0x6CD);
  if (operation == "all" || operation == "gcd")
    for (const auto shape : {"balanced", "skewed", "common_factor"})
      for (unsigned bits : {256, 512, 1024, 2048, 4096})
        BenchmarkGcd(bits, shape, rng);
  if (operation == "gcd-large")
    for (const auto shape : {"balanced", "skewed", "common_factor", "fibonacci"})
      for (unsigned bits : {2048, 4096, 8192, 16384, 32768, 65536, 262144})
        BenchmarkGcd(bits, shape, rng, std::string_view(shape) == "fibonacci" ? 1 : 8);
  rng.seed(0x5A27);
  if (operation == "all" || operation == "sqrtrem")
    for (unsigned bits : {32, 64, 65, 127, 128, 129, 192, 254, 255, 256, 257,
                          1024, 4096, 16384, 65536, 262144, 1048576})
      BenchmarkSqrtRem(bits, rng);
  rng.seed(0xAB1979);
  if (operation == "all" || operation == "prime")
  {
    for (unsigned bits : {32, 64, 128, 256, 512, 1024, 2048, 4096})
      for (const auto shape : {"prime", "small_factor", "semiprime"})
        BenchmarkPrimality(bits, shape, rng);
    BenchmarkPseudoprimes();
  }
}
