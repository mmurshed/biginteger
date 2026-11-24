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
    mpz_t a, b, q, r, check, exactQuotient;
    gmp_randstate_t rng;

    Fixture()
    {
      mpz_inits(a, b, q, r, check, exactQuotient, nullptr);
      gmp_randinit_mt(rng);
      gmp_randseed_ui(rng, 0x1964);
    }

    ~Fixture()
    {
      gmp_randclear(rng);
      mpz_clears(a, b, q, r, check, exactQuotient, nullptr);
    }

    void Generate(mpz_t out, unsigned long bits)
    {
      mpz_urandomb(out, rng, bits);
      mpz_setbit(out, bits - 1);
    }

    std::vector<DataT> Export(mpz_t value)
    {
      std::vector<DataT> limbs((mpz_sizeinbase(value, 2) + LimbBits - 1) / LimbBits);
      size_t count = 0;
      mpz_export(limbs.data(), &count, -1, sizeof(DataT), 0, 64 - LimbBits, value);
      limbs.resize(std::max<size_t>(1, count));
      return limbs;
    }

    bool Matches(std::span<const DataT> value, mpz_t expected)
    {
      mpz_import(check, value.size(), -1, sizeof(DataT), 0, 64 - LimbBits, value.data());
      return mpz_cmp(check, expected) == 0;
    }
  };
}

int main(int argc, char **argv)
{
  if (argc != 8 && argc != 9)
  {
    std::fprintf(stderr, "usage: division_handoff_bench A_LIMBS B_LIMBS dispatch|bz|exact|gmp-exact "
                         "q|qr random|exact SHIFT WARM_REPS [DIVISOR_V2]\n");
    return 2;
  }
  const auto na = std::stoul(argv[1]), nb = std::stoul(argv[2]);
  const std::string_view path = argv[3], mode = argv[4], input = argv[5];
  const auto shift = std::stoul(argv[6]);
  const int reps = std::stoi(argv[7]);
  const bool exactPath = path == "exact" || path == "gmp-exact";
  if (nb == 0 || na < nb || (na == nb && input != "exact") || shift >= LimbBits ||
      (path != "dispatch" && path != "bz" && !exactPath) || (mode != "q" && mode != "qr") ||
      (input != "random" && input != "exact") || (exactPath && (input != "exact" || mode != "q")))
    throw std::invalid_argument("invalid division benchmark case");

  const bool wantRemainder = mode == "qr";
  Fixture fixture;
  fixture.Generate(fixture.a, na * LimbBits);
  fixture.Generate(fixture.b, nb * LimbBits - shift);
  if (argc == 9)
  {
    const auto zeros = std::stoul(argv[8]);
    if (zeros >= nb * LimbBits - shift)
      throw std::invalid_argument("divisor has too many trailing zero bits");
    mpz_tdiv_q_2exp(fixture.b, fixture.b, zeros);
    mpz_setbit(fixture.b, 0);
    mpz_mul_2exp(fixture.b, fixture.b, zeros);
  }
  mpz_tdiv_qr(fixture.q, fixture.r, fixture.a, fixture.b);
  if (input == "exact")
  {
    if (mpz_sgn(fixture.q) == 0)
      mpz_set_ui(fixture.q, 1);
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
    if (path == "gmp-exact")
      mpz_divexact(fixture.exactQuotient, fixture.a, fixture.b);
    else if (path == "exact")
      result.first = DivideExact(a, b);
    else if (path == "bz")
      limbs = BurnikelZieglerDivision::DivideAndRemainder(av, bv, BigInteger::Base(), wantRemainder);
    else if (wantRemainder)
      result = DivideAndRemainder(a, b);
    else
      result.first = a / b;
  });

  const std::span<const DataT> q = path == "bz" ? std::span<const DataT>(limbs.first) : result.first.Limbs();
  const std::span<const DataT> r = path == "bz" ? std::span<const DataT>(limbs.second) : result.second.Limbs();
  const bool quotientMatches = path == "gmp-exact" ? mpz_cmp(fixture.exactQuotient, fixture.q) == 0
                                                    : fixture.Matches(q, fixture.q);
  if (!quotientMatches || (wantRemainder && !fixture.Matches(r, fixture.r)))
    throw std::runtime_error("division benchmark disagrees with GMP");
  std::printf("%.9g,%.9g,%llu,%llu,%llu\n", time.first_ms, time.warm_ms,
              static_cast<unsigned long long>(time.setup_peak),
              static_cast<unsigned long long>(time.first_peak),
              static_cast<unsigned long long>(time.warm_peak));
}
