// Small-operand benchmark: BigInteger operators against GMP mpz on operands
// of 1-16 limbs, where call overhead and allocation decide the cost.
//
// Usage: small_ops_vs_gmp [label]   (label goes in the first CSV column)
//
// BigMath is timed through the public operators, which return a new
// BigInteger each call. GMP is timed writing into an output that is reused
// across calls (mpz_mul(c, a, b) in a loop), its idiomatic use. Each figure
// is the minimum over 7 rounds of the mean time per call.
//
// Operands: a and b have n full limbs for +, -, *; for / and % the dividend
// has 2n limbs and the divisor n. Only the public BigInteger API and limb
// indexing are used, so the file builds against older trees too.

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include <gmp.h>

#include "biginteger/BigInteger.h"
#include "biginteger/ops/Operations.h"

using namespace BigMath;
using Clock = std::chrono::steady_clock;

namespace
{
  volatile uint64_t sink = 0;

  template <class F>
  double NsPerCall(F &&f, long calls)
  {
    double best = 1e300;
    for (int round = 0; round < 7; ++round)
    {
      auto t0 = Clock::now();
      for (long i = 0; i < calls; ++i)
        f();
      best = std::min(best, std::chrono::duration<double, std::nano>(Clock::now() - t0).count() / calls);
    }
    return best;
  }

  std::vector<DataT> RandomLimbs(SizeT n, std::mt19937_64 &gen)
  {
    std::vector<DataT> v(n);
    for (auto &x : v)
    {
#if BIGMATH_LIMB_64
      x = (DataT)gen();
#else
      x = (DataT)(gen() & 0xFFFFFFFFULL);
#endif
    }
    v.back() |= 1;
    return v;
  }

  void ToMpz(mpz_t z, std::vector<DataT> const &v)
  {
#if BIGMATH_LIMB_64
    mpz_import(z, v.size(), -1, sizeof(DataT), 0, 0, v.data());
#else
    mpz_import(z, v.size(), -1, sizeof(DataT), 0, 32, v.data());
#endif
  }
}

int main(int argc, char **argv)
{
  const std::string label = argc > 1 ? argv[1] : "bigmath";
  std::mt19937_64 gen(0x5A11ULL);
  std::printf("version,op,limbs,bigmath_ns,gmp_ns,ratio\n");

  for (SizeT n : {1, 2, 3, 4, 6, 8, 12, 16})
  {
    auto av = RandomLimbs(n, gen), bv = RandomLimbs(n, gen), dv = RandomLimbs(2 * n, gen);
    BigInteger a(av, false), b(bv, false), d(dv, false);
    mpz_t ga, gb, gd, gc;
    mpz_inits(ga, gb, gd, gc, NULL);
    ToMpz(ga, av);
    ToMpz(gb, bv);
    ToMpz(gd, dv);

    const long calls = 400000;
    struct Row
    {
      const char *op;
      double bm, gmp;
    };
    Row rows[] = {
        {"mul", NsPerCall([&] { BigInteger c = a * b; sink += c[0]; }, calls),
         NsPerCall([&] { mpz_mul(gc, ga, gb); sink += mpz_getlimbn(gc, 0); }, calls)},
        {"add", NsPerCall([&] { BigInteger c = a + b; sink += c[0]; }, calls),
         NsPerCall([&] { mpz_add(gc, ga, gb); sink += mpz_getlimbn(gc, 0); }, calls)},
        {"sub", NsPerCall([&] { BigInteger c = a - b; sink += c[0]; }, calls),
         NsPerCall([&] { mpz_sub(gc, ga, gb); sink += mpz_getlimbn(gc, 0); }, calls)},
        {"div", NsPerCall([&] { BigInteger c = d / b; sink += c[0]; }, calls),
         NsPerCall([&] { mpz_tdiv_q(gc, gd, gb); sink += mpz_getlimbn(gc, 0); }, calls)},
        {"mod", NsPerCall([&] { BigInteger c = d % b; sink += c[0]; }, calls),
         NsPerCall([&] { mpz_tdiv_r(gc, gd, gb); sink += mpz_getlimbn(gc, 0); }, calls)},
        {"copy", NsPerCall([&] { BigInteger c = a; sink += c[0]; }, calls),
         NsPerCall([&] { mpz_set(gc, ga); sink += mpz_getlimbn(gc, 0); }, calls)},
    };
    for (auto const &r : rows)
      std::printf("%s,%s,%u,%.2f,%.2f,%.3f\n", label.c_str(), r.op, (unsigned)n, r.bm, r.gmp, r.bm / r.gmp);
    mpz_clears(ga, gb, gd, gc, NULL);
  }
  std::fprintf(stderr, "checksum %llu\n", (unsigned long long)sink);
  return 0;
}
