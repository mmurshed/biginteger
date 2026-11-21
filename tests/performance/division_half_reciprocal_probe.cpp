// Compare a guarded half-size reciprocal with the public dispatcher at near-2:1
// shapes. This experimental divisor is used only by this benchmark.
#include "biginteger/algorithms/Division.h"
#include "biginteger/algorithms/division/NewtonDivision.h"

#include <gmp.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <random>

using namespace BigMath;

namespace
{
  using Limbs = std::vector<DataT>;
  volatile DataT checksum;

  // Inputs are normalized 64-bit magnitudes. Start with a remainder smaller
  // than the divisor and extend it by at most half a divisor per quotient block.
  // Dividing the operand tops can overestimate a digit. Comparing the complete
  // back-product recovers the exact digit and remainder.
  std::pair<Limbs, Limbs> HalfReciprocalQR(const Limbs &a, const Limbs &divisor)
  {
    size_t n = divisor.size(), block = (n + 1) / 2, drop = n - (block + 2);
    NewtonDivision::Divider divider(Limbs(divisor.begin() + drop, divisor.end()), Base2_64);
    size_t pos = a.size() - n;
    Limbs remainder(a.begin() + pos, a.end()), quotient(pos + 1);
    if (Compare(remainder, divisor) >= 0)
    {
      remainder = Subtract(remainder, divisor, Base2_64);
      quotient[pos] = 1;
    }
    while (pos)
    {
      size_t take = std::min(pos, block);
      Limbs chunk(a.begin() + pos - take, a.begin() + pos);
      chunk.insert(chunk.end(), remainder.begin(), remainder.end());
      TrimZerosToOne(chunk);
      pos -= take;
      Limbs digit{0};
      if (Compare(chunk, divisor) >= 0)
      {
        digit = divider.Divide(Limbs(chunk.begin() + drop, chunk.end()));
        Limbs product = Multiply(digit, divisor, Base2_64);
        while (Compare(product, chunk) > 0)
        {
          digit = Subtract(digit, Limbs{1}, Base2_64);
          product = Subtract(product, divisor, Base2_64);
        }
        remainder = Subtract(chunk, product, Base2_64);
      }
      else
        remainder = std::move(chunk);
      if (Compare(remainder, divisor) >= 0 || digit.size() > take)
        throw std::runtime_error("invalid quotient block");
      std::copy(digit.begin(), digit.end(), quotient.begin() + pos);
    }
    TrimZerosToOne(quotient);
    TrimZerosToOne(remainder);
    return {std::move(quotient), std::move(remainder)};
  }

  template<class Fn>
  double MedianMicros(Fn fn)
  {
    std::array<double, 7> times;
    for (double &us : times)
    {
      auto start = std::chrono::steady_clock::now();
      checksum = fn();
      us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - start).count();
    }
    std::sort(times.begin(), times.end());
    return times[times.size() / 2];
  }

  void Import(mpz_t out, const Limbs &limbs)
  {
    mpz_import(out, limbs.size(), -1, sizeof(DataT), 0, 0, limbs.data());
  }
}

int main()
{
  if (LimbBits != 64)
  {
    std::cerr << "This probe requires 64-bit limbs.\n";
    return 2;
  }
  std::mt19937_64 gen(129);
  auto operand = [&](size_t n) {
    Limbs limbs(n);
    for (auto &limb : limbs) limb = gen();
    limbs.back() |= DataT{1} << 63;
    return limbs;
  };
  auto setup = [](const Limbs &d) {
    return MedianMicros([&]() { return NewtonDivision::ApproxReciprocal(d)[0]; });
  };
  mpz_t ga, gb, gq, gr, check;
  mpz_inits(ga, gb, gq, gr, check, nullptr);
  std::cout << "divisor_limbs,dividend_limbs,full_setup_us,half_setup_us,dispatch_qr_us,half_qr_us,gmp_qr_us\n";
  for (int n : {1023, 1024, 1025, 1536, 2047, 2048, 2049, 4095, 4096, 4097})
  {
    Limbs d = operand(n);
    Import(gb, d);
    double fullSetup = setup(d), halfSetup = setup(Limbs(d.end() - ((n + 1) / 2 + 2), d.end()));
    for (int na : {7 * n / 4, 2 * n, 9 * n / 4})
    {
      Limbs a = operand(na);
      Import(ga, a);
      auto qr = HalfReciprocalQR(a, d);
      if (qr != DivideAndRemainder(a, d, Base2_64, true))
        throw std::runtime_error("quotient or remainder differs from dispatcher");
      mpz_tdiv_qr(gq, gr, ga, gb);
      Import(check, qr.first);
      if (mpz_cmp(check, gq)) throw std::runtime_error("quotient differs from GMP");
      Import(check, qr.second);
      if (mpz_cmp(check, gr)) throw std::runtime_error("remainder differs from GMP");
      double half = MedianMicros([&]() { return HalfReciprocalQR(a, d).first[0]; });
      double dispatch = MedianMicros([&]() { return DivideAndRemainder(a, d, Base2_64, true).first[0]; });
      double gmp = MedianMicros([&]() { mpz_tdiv_qr(gq, gr, ga, gb); return mpz_get_ui(gq); });
      std::cout << n << ',' << na << ',' << fullSetup << ',' << halfSetup << ','
                << dispatch << ',' << half << ',' << gmp << '\n';
    }
  }
  mpz_clears(ga, gb, gq, gr, check, nullptr);
}
