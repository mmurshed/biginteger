#ifndef BIGMATH_BENCH_FLINT_INTEGER
#define BIGMATH_BENCH_FLINT_INTEGER

#include <gmp.h>
#include <flint/fmpz.h>

#include "biginteger/BigInteger.h"

namespace bench
{
  struct FlintInteger
  {
    fmpz_t value;

    FlintInteger() { fmpz_init(value); }
    ~FlintInteger() { fmpz_clear(value); }
    FlintInteger(FlintInteger const &) = delete;
    FlintInteger &operator=(FlintInteger const &) = delete;

    void Set(BigMath::BigInteger const &source)
    {
      const auto bytes = source.ToByteArray();
      mpz_t word;
      mpz_init(word);
      mpz_import(word, bytes.size(), 1, 1, 0, 0, bytes.data());
      if (source.IsNegative())
        mpz_neg(word, word);
      fmpz_set_mpz(value, word);
      mpz_clear(word);
    }

    BigMath::BigInteger Get() const
    {
      std::vector<uint8_t> bytes((fmpz_bits(value) + 7) / 8);
      mpz_t word;
      mpz_init(word);
      fmpz_get_mpz(word, value);
      size_t count = 0;
      mpz_export(bytes.data(), &count, 1, 1, 0, 0, word);
      mpz_clear(word);
      return BigMath::BigInteger::FromByteArray(std::span(bytes).first(count), fmpz_sgn(value) < 0);
    }
  };
}

#endif
