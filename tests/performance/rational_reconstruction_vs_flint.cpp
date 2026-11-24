#include <cstdio>
#include <random>
#include <stdexcept>

#include "biginteger/ops/Operations.h"
#include "flint_integer.h"
#include "paired_timing.h"

#include <flint/fmpq.h>

namespace
{
  using namespace BigMath;
  using bench::FlintInteger;

  BigInteger PowerOfTwo(unsigned bit)
  {
    const DataT top = DataT{1} << (bit % LimbBits);
    return BigInteger::FromLimbs({&top, 1}, false) << (bit / LimbBits);
  }

  BigInteger RandomMagnitude(unsigned bits, std::mt19937_64 &rng, bool negative, bool odd)
  {
    std::vector<uint8_t> bytes((bits + 7) / 8);
    for (auto &byte : bytes)
      byte = static_cast<uint8_t>(rng());
    const unsigned high = 1u << ((bits - 1) % 8);
    bytes.front() = (bytes.front() & (high - 1)) | high;
    if (odd)
      bytes.back() |= 1;
    return BigInteger::FromByteArray(bytes, negative);
  }

  void Benchmark(unsigned bits, unsigned numeratorBits, std::mt19937_64 &rng)
  {
    const auto one = PowerOfTwo(0);
    const auto modulus = PowerOfTwo(bits - 1);
    const unsigned denominatorBits = bits - numeratorBits - 3;
    const auto maxNumerator = PowerOfTwo(numeratorBits) - one;
    const auto maxDenominator = PowerOfTwo(denominatorBits) - one;

    FlintInteger m, maxN, maxD, a, n, d, common, inverse, refN, refD;
    m.Set(modulus);
    maxN.Set(maxNumerator);
    maxD.Set(maxDenominator);
    BigInteger residue;
    std::optional<RationalReconstructionResult> result;
    int found = 0;
    const auto runBigMath = [&] { result = ReconstructRational(residue, modulus, maxNumerator, maxDenominator); };
    const auto runFlint = [&] {
      found = _fmpq_reconstruct_fmpz_2(refN.value, refD.value, a.value, m.value, maxN.value, maxD.value);
    };

    constexpr unsigned checks = 8;
    for (unsigned i = 0; i < checks; ++i)
    {
      const auto numerator = RandomMagnitude(numeratorBits, rng, i % 2 != 0, false);
      BigInteger denominator;
      n.Set(numerator);
      do
      {
        denominator = RandomMagnitude(denominatorBits, rng, false, true);
        d.Set(denominator);
        fmpz_gcd(common.value, n.value, d.value);
      } while (!fmpz_is_one(common.value));

      if (!fmpz_invmod(inverse.value, d.value, m.value))
        throw std::runtime_error("fixture has no modular inverse");
      fmpz_mul(a.value, n.value, inverse.value);
      fmpz_mod(a.value, a.value, m.value);
      residue = a.Get();
      runBigMath();
      runFlint();
      if (!result || !found || result->numerator != numerator || result->denominator != denominator ||
          result->numerator != refN.Get() || result->denominator != refD.Get())
        throw std::runtime_error("rational reconstruction differs from FLINT");
    }

    const auto calls = std::clamp(16384u / bits, 1u, 32u);
    const auto [bigmath, flint] = bench::MeasureTimes(calls, runBigMath, runFlint);
    std::printf("%u,%u,%u,%u,%u,%.6f,%.6f\n", bits, numeratorBits, denominatorBits,
                checks, calls, bigmath, flint);
  }
}

int main()
{
  std::mt19937_64 rng(0x76F17);
  std::puts("modulus_bits,numerator_bound_bits,denominator_bound_bits,checks,calls,bigmath_us,flint_us");
  for (unsigned bits : {1024u, 16384u, 65536u, 262144u})
    for (unsigned numeratorBits : {bits / 2 - 2, bits / 8})
      Benchmark(bits, numeratorBits, rng);
}
