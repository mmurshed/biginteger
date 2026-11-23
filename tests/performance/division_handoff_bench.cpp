// One case per process so retained workspace cannot carry across cases.
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string_view>

#include <gmp.h>

#include "benchsuite/case_measurement.h"
#include "biginteger/ops/Division.h"

using namespace BigMath;

namespace
{
  struct Fixture
  {
    mpz_t a, b, q, r, check;
    gmp_randstate_t rng;

    Fixture()
    {
      mpz_inits(a, b, q, r, check, nullptr);
      gmp_randinit_mt(rng);
      gmp_randseed_ui(rng, 0x1964);
    }

    ~Fixture()
    {
      gmp_randclear(rng);
      mpz_clears(a, b, q, r, check, nullptr);
    }

    void Generate(mpz_t out, unsigned long bits)
    {
      mpz_urandomb(out, rng, bits);
      mpz_setbit(out, bits - 1);
    }

    std::vector<DataT> Export(mpz_t value)
    {
      std::vector<DataT> limbs((mpz_sizeinbase(value, 2) + 63) / 64);
      size_t count = 0;
      mpz_export(limbs.data(), &count, -1, sizeof(DataT), 0, 0, value);
      limbs.resize(std::max<size_t>(1, count));
      return limbs;
    }

    bool Matches(std::span<const DataT> value, mpz_t expected)
    {
      mpz_import(check, value.size(), -1, sizeof(DataT), 0, 0, value.data());
      return mpz_cmp(check, expected) == 0;
    }
  };
}

int main(int argc, char **argv)
{
  if (argc != 8 || BigInteger::Base() != Base2_64)
  {
    std::fprintf(stderr, "usage: division_handoff_bench A_LIMBS B_LIMBS dispatch|bz q|qr random|exact SHIFT WARM_REPS\n");
    return 2;
  }
  const auto na = std::stoul(argv[1]), nb = std::stoul(argv[2]);
  const std::string_view path = argv[3], mode = argv[4], input = argv[5];
  const auto shift = std::stoul(argv[6]);
  const int reps = std::stoi(argv[7]);
  if (na <= nb || nb == 0 || shift >= 64 ||
      (path != "dispatch" && path != "bz") || (mode != "q" && mode != "qr") ||
      (input != "random" && input != "exact"))
    throw std::invalid_argument("invalid division benchmark case");

  const bool wantRemainder = mode == "qr";
  Fixture fixture;
  fixture.Generate(fixture.a, na * 64);
  fixture.Generate(fixture.b, nb * 64 - shift);
  mpz_tdiv_qr(fixture.q, fixture.r, fixture.a, fixture.b);
  if (input == "exact")
  {
    mpz_mul(fixture.a, fixture.q, fixture.b);
    mpz_set_ui(fixture.r, 0);
  }
  const auto av = fixture.Export(fixture.a), bv = fixture.Export(fixture.b);
  if (av.size() != na || bv.size() != nb)
    throw std::runtime_error("generated operand width changed");
  const BigInteger a(av, false), b(bv, false);

  std::pair<std::vector<DataT>, std::vector<DataT>> limbs;
  std::pair<BigInteger, BigInteger> result;
  bench::CaseMeasurement time;
#if defined(__APPLE__) || defined(__linux__)
  time.memory = true;
#endif
  time.Run(reps, [&] {
    if (path == "bz")
      limbs = BurnikelZieglerDivision::DivideAndRemainder(av, bv, Base2_64, wantRemainder);
    else if (wantRemainder)
      result = DivideAndRemainder(a, b);
    else
      result.first = a / b;
  });

  const std::span<const DataT> q = path == "bz" ? std::span<const DataT>(limbs.first) : result.first.Limbs();
  const std::span<const DataT> r = path == "bz" ? std::span<const DataT>(limbs.second) : result.second.Limbs();
  if (!fixture.Matches(q, fixture.q) || (wantRemainder && !fixture.Matches(r, fixture.r)))
    throw std::runtime_error("division benchmark disagrees with GMP");
  std::printf("%.9g,%.9g,%llu,%llu,%llu\n", time.first_ms, time.warm_ms,
              static_cast<unsigned long long>(time.setup_peak),
              static_cast<unsigned long long>(time.first_peak),
              static_cast<unsigned long long>(time.warm_peak));
}
