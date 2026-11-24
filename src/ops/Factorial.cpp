#include "Factorial.h"

#include <array>
#include <bit>

#include "biginteger/common/BitShifts.h"
#include "biginteger/common/Builder.h"
#include "biginteger/ops/Multiplication.h"

namespace BigMath
{
  namespace
  {
    BigInteger ProductTree(std::span<const DataT> factors)
    {
      if (factors.size() == 1)
        return BigInteger::FromLimbs(factors, false);
      const auto mid = factors.size() / 2;
      return ProductTree(factors.first(mid)) * ProductTree(factors.subspan(mid));
    }

    class PackedProduct
    {
      std::vector<DataT> factors;
      DataT pending = 1;

    public:
      void Append(DataT factor)
      {
        const ULong128 combined = ULong128{pending} * factor;
        if (combined <= LimbMask)
          pending = static_cast<DataT>(combined);
        else
        {
          factors.push_back(pending);
          pending = factor;
        }
      }

      BigInteger Product()
      {
        if (factors.empty())
          return BigInteger::FromLimbs({&pending, 1}, false);
        factors.push_back(pending);
        return ProductTree(factors);
      }
    };

    std::vector<UInt> SieveOddPrimes(UInt n)
    {
      std::vector<bool> composite(n / 2 + 1);
      std::vector<UInt> primes;
      for (ULong p = 3; p <= n; p += 2)
      {
        if (composite[p / 2])
          continue;
        primes.push_back(static_cast<UInt>(p));
        if (p <= n / p)
          for (ULong multiple = p * p; multiple <= n; multiple += 2 * p)
            composite[multiple / 2] = true;
      }
      return primes;
    }

    BigInteger OddFactorial(UInt n, std::span<const UInt> oddPrimes)
    {
      if (n <= 20)
      {
        ULong oddPart = 1;
        for (UInt k = 2; k <= n; ++k)
          oddPart *= k >> std::countr_zero(k);
        return BigIntegerBuilder::From(oddPart);
      }

      const auto mid = n / 2;
      const auto halfPrimes = std::upper_bound(oddPrimes.begin(), oddPrimes.end(), mid);
      auto half = OddFactorial(mid, oddPrimes.first(halfPrimes - oddPrimes.begin()));
      PackedProduct swing;
      for (UInt p : oddPrimes)
      {
        // The exponent in n! / floor(n/2)!^2 sums floor(n/p^i) mod 2.
        if (p > mid)
          swing.Append(p);
        else
          for (UInt quotient = n / p; quotient != 0; quotient /= p)
            if (quotient & 1)
              swing.Append(p);
      }
      return (half * half) * swing.Product();
    }
  }

  BigInteger detail::FactorialProductTree(UInt n)
  {
    PackedProduct product;
    for (ULong factor = 2; factor <= n; ++factor)
      product.Append(factor);
    return product.Product();
  }

  BigInteger Factorial(UInt n)
  {
    static constexpr auto small = [] {
      std::array<ULong, 21> values{1};
      for (std::size_t i = 1; i < values.size(); ++i)
        values[i] = values[i - 1] * i;
      return values;
    }();
    if (n < small.size())
      return BigIntegerBuilder::From(small[n]);

    constexpr UInt minSwing = LimbBits == 64 ? 512 : 64;
    if (n < minSwing)
      return detail::FactorialProductTree(n);
    const auto oddPart = OddFactorial(n, SieveOddPrimes(n));
    const auto twos = n - std::popcount(n);
    auto limbs = ShiftLeftBits(oddPart.Limbs(), twos % LimbBits, LimbBits);
    limbs.insert(limbs.begin(), twos / LimbBits, 0);
    return BigInteger(std::move(limbs), false);
  }
}
