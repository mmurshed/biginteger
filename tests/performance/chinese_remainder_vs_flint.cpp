#include <cstdio>
#include <random>
#include <stdexcept>

#include "biginteger/ops/Operations.h"
#include "flint_integer.h"
#include "paired_timing.h"

namespace
{
  using namespace BigMath;
  using bench::FlintInteger;

  struct FlintPlan
  {
    fmpz_comb_t plan;
    fmpz_comb_temp_t scratch;

    explicit FlintPlan(std::span<const ulong> moduli)
    {
      fmpz_comb_init(plan, moduli.data(), static_cast<slong>(moduli.size()));
      fmpz_comb_temp_init(scratch, plan);
    }

    ~FlintPlan()
    {
      fmpz_comb_temp_clear(scratch);
      fmpz_comb_clear(plan);
    }

    FlintPlan(FlintPlan const &) = delete;
    FlintPlan &operator=(FlintPlan const &) = delete;
  };

  void Benchmark(unsigned count, unsigned wordBits, std::mt19937_64 &rng)
  {
    std::vector<ulong> flintModuli(count);
    FlintInteger prime;
    fmpz_set_ui(prime.value, ulong{1} << (wordBits - 1));
    for (auto &modulus : flintModuli)
    {
      fmpz_nextprime(prime.value, prime.value, 0);
      modulus = fmpz_get_ui(prime.value);
    }
    const std::vector<DataT> moduli(flintModuli.begin(), flintModuli.end());
    const ChineseRemainderPlan plan(moduli);
    FlintPlan flint(flintModuli);
    std::vector<uint8_t> digits(plan.Modulus().ToByteArray().size());
    BigInteger input, reconstructed;
    std::vector<DataT> residues;
    std::vector<ulong> flintResidues(count);
    FlintInteger fInput, fResult;
    const auto reduce = [&] { residues = plan.Reduce(input); };
    const auto reduceFlint = [&] { fmpz_multi_mod_ui(flintResidues.data(), fInput.value, flint.plan, flint.scratch); };
    const auto reconstruct = [&] { reconstructed = plan.Reconstruct(residues); };
    const auto reconstructFlint = [&] {
      fmpz_multi_CRT_ui(fResult.value, flintResidues.data(), flint.plan, flint.scratch, 0);
    };

    constexpr unsigned checks = 8;
    for (unsigned i = 0; i < checks; ++i)
    {
      for (auto &digit : digits)
        digit = static_cast<uint8_t>(rng());
      input = BigInteger::FromByteArray(digits, i % 2 == 0);
      // Check signed values and magnitudes above the product before
      // timing a nonnegative representative.
      const auto signedResidues = plan.Reduce(input);
      input = input % plan.Modulus();
      if (input.IsNegative())
        input = input + plan.Modulus();
      fInput.Set(input);
      reduce();
      reduceFlint();
      reconstruct();
      reconstructFlint();
      if (residues != signedResidues ||
          !std::equal(residues.begin(), residues.end(), flintResidues.begin()) ||
          reconstructed != input || !fmpz_equal(fResult.value, fInput.value))
        throw std::runtime_error("Chinese remaindering differs from FLINT");
    }

    const auto calls = std::clamp(16384u / count, 2u, 512u);
    const auto setup = bench::MeasureTimes(calls, [&] { const ChineseRemainderPlan state(moduli); },
                                                [&] { const FlintPlan state(flintModuli); });
    const auto reduction = bench::MeasureTimes(calls, reduce, reduceFlint);
    const auto reconstruction = bench::MeasureTimes(calls, reconstruct, reconstructFlint);
    std::printf("%u,%u,%u,%u,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", count, wordBits, checks, calls,
                setup.first, setup.second, reduction.first, reduction.second,
                reconstruction.first, reconstruction.second);
  }
}

int main()
{
  std::mt19937_64 rng(0xC771218);
  std::puts("moduli,word_bits,checks,calls,bigmath_setup_us,flint_setup_us,"
            "bigmath_reduce_us,flint_reduce_us,bigmath_reconstruct_us,flint_reconstruct_us");
  for (unsigned wordBits : {31u, 61u})
    for (unsigned count : {8u, 128u, 1024u, 8192u})
      Benchmark(count, wordBits, rng);
}
