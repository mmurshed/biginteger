#include "unit_test_framework.h"

#include <limits>
#include <random>

#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Subtraction.h"
#include "biginteger/algorithms/division/FastDivision.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/algorithms/multiplication/FermatOuterMultiplication.h"
#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"

using namespace BigMath;
using namespace BigMath::FermatOuter;

namespace
{
  using Limbs = std::vector<DataT>;

  Limbs RingReference(Limbs value, const Limbs &modulus)
  {
    TrimZerosToOne(value);
    auto out = FastDivision::DivideAndRemainder(value, modulus, Base2_64).second;
    out.resize(modulus.size());
    return out;
  }

  Limbs RandomValue(size_t limbs, std::mt19937_64 &rng)
  {
    Limbs value(limbs);
    for (auto &x : value) x = rng();
    return value;
  }
}

REGISTER_TEST(FermatOuter, RingArithmeticMatchesDivision)
{
  std::mt19937_64 rng(0xe5a74);
  PointwiseScratch scratch;
  for (size_t n : {1, 2, 3, 5, 8})
  {
    Limbs modulus(n + 1), zero(n + 1), one(n + 1), negativeOne(n + 1), maximal(n + 1, ~DataT{0});
    modulus.front() = modulus.back() = one.front() = negativeOne.back() = 1;
    maximal.back() = 0;
    std::vector<Limbs> values{zero, one, negativeOne, maximal};
    for (int i = 0; i < 15; ++i)
    {
      auto value = RandomValue(n, rng);
      value.push_back(0);
      values.push_back(std::move(value));
    }

    for (auto top : {DataT{0}, DataT{1}, DataT{2}, DataT{3}, ~DataT{0}})
      for (const auto &value : values)
      {
        auto a = value;
        a.back() = top;
        const auto want = RingReference(a, modulus);
        Normalize(a);
        ASSERT_TRUE(a == want);
      }

    for (const auto &a : values)
      for (const auto &b : values)
      {
        auto x = a, y = b;
        TrimZerosToOne(x);
        TrimZerosToOne(y);
        Limbs out(n + 1);
        const auto sum = RingReference(Add(x, y, Base2_64), modulus);
        const auto difference = RingReference(Subtract(Add(x, modulus, Base2_64), y, Base2_64), modulus);
        AddMod(out, a, b);
        ASSERT_TRUE(out == sum);
        SubMod(out, a, b);
        ASSERT_TRUE(out == difference);
        out = a;
        AddMod(out, out, b);
        ASSERT_TRUE(out == sum);
        out = b;
        SubMod(out, a, out);
        ASSERT_TRUE(out == difference);
        const auto product = RingReference(ClassicMultiplication::Multiply(x, y, Base2_64), modulus);
        out = a;
        MulMod(out, out, b, scratch);
        ASSERT_TRUE(out == product);
        out = b;
        MulMod(out, a, out, scratch);
        ASSERT_TRUE(out == product);
      }

    for (const auto &a : values)
      for (size_t shift = 0; shift < 128 * n; ++shift)
      {
        Limbs power(1 + shift / 64);
        power.back() = DataT{1} << (shift % 64);
        auto x = a;
        TrimZerosToOne(x);
        const auto want = RingReference(ClassicMultiplication::Multiply(x, power, Base2_64), modulus);
        Limbs out(n + 1);
        MulPow2Mod(out, a, shift);
        ASSERT_TRUE(out == want);
        auto reduced = RingReference(power, modulus);
        MulMod(reduced, a, reduced, scratch);
        ASSERT_TRUE(out == reduced);
      }
  }
}

REGISTER_TEST(FermatOuter, PointProductsReuseStorageAcrossWidths)
{
  std::mt19937_64 rng(0x70a4f);
  constexpr size_t width = 4096;
  Limbs modulus(width + 1), out(width + 1);
  modulus.front() = modulus.back() = 1;
  PointwiseScratch scratch;

  auto a = RandomValue(width, rng), b = RandomValue(width, rng);
  auto expected = RingReference(KaratsubaMultiplication::Multiply(a, b, Base2_64), modulus);
  a.push_back(0);
  b.push_back(0);
  MulMod(out, a, b, scratch);
  ASSERT_TRUE(out == expected);

  const auto capacity = scratch.product.capacity();
  const Limbs zero(width + 1);
  MulMod(out, a, zero, scratch);
  ASSERT_TRUE(out == zero);
  ASSERT_TRUE(scratch.product.empty());

  // The ring slots stay wide even when the next product fits the scalar,
  // schoolbook or recursive kernel. Include an in-place aliased square.
  for (size_t limbs : {1, 8, 137})
  {
    auto x = RandomValue(limbs, rng);
    expected = RingReference(ClassicMultiplication::Multiply(x, x, Base2_64), modulus);
    x.resize(width + 1);
    MulMod(x, x, x, scratch);
    ASSERT_TRUE(x == expected);
    ASSERT_EQ(scratch.product.capacity(), capacity);
  }
}

REGISTER_TEST(FermatOuter, RoundTripsSingleAndParallelLayers)
{
  std::mt19937_64 rng(0x906c);
  for (const auto &p : {Parameters{1, 2, 1, 1}, Parameters{3, 8, 1, 3},
                       Parameters{6, 64, 1, 5}, Parameters{8, 256, 1, 6},
                       Parameters{10, 1024, 1, 8}, Parameters{6, 64, 4096, 8193}})
  {
    const size_t stride = p.ringLimbs + 1;
    auto plane = RandomValue(p.count * stride, rng);
    for (size_t i = 0; i < p.count; ++i)
    {
      plane[(i + 1) * stride - 1] = i % 3 == 0;
      if (i % 3 == 0)
        std::fill_n(plane.data() + i * stride, p.ringLimbs, 0); // -1
    }
    const auto original = plane;
    Transform fft(p);
    fft.Forward(plane);
    fft.Inverse(plane);
    ASSERT_TRUE(plane == original);
  }
}

REGISTER_TEST(FermatOuter, ProductsMatchIndependentMultiplication)
{
  std::mt19937_64 rng(0x5e147);
  for (unsigned exponent : {3, 4, 5, 6, 8})
    for (size_t piece : {1, 2, 3, 5, 8})
    {
      const size_t count = size_t{1} << exponent;
      const size_t half = count / 2;
      for (const auto &[na, nb] : {std::pair{half * piece, half * piece + 1},
                                  std::pair{piece, (count - 1) * piece},
                                  std::pair{(half - 1) * piece + 1, (half - 1) * piece + 1}})
      {
        const auto p = MakeParameters(na, nb, exponent, piece);
        ASSERT_TRUE(p.has_value());
        auto a = RandomValue(na, rng), b = RandomValue(nb, rng);
        std::fill(b.begin(), b.end(), ~DataT{0});
        if (exponent >= 6)
          for (size_t i = 0; i + 1 < na; i += 2) a[i] = 0;
        const auto want = KaratsubaMultiplication::Multiply(a, b, Base2_64);
        ASSERT_TRUE(MultiplyFermat(a, b, *p) == want);
        ASSERT_TRUE(MultiplyFermat(b, a, *p) == want);
        ASSERT_TRUE(MultiplyFermat(Limbs(na), b, *p).empty());

        if (const auto square = MakeParameters(na, na, exponent, piece))
        {
          std::fill(a.begin(), a.end(), ~DataT{0});
          ASSERT_TRUE(MultiplyFermat(a, a, *square) == KaratsubaMultiplication::Multiply(a, a, Base2_64));
        }
      }
    }
}

REGISTER_TEST(FermatOuter, SelectsParametersAndRejectsOutOfRangeSizes)
{
  ASSERT_TRUE(!MakeParameters(5, 5, 3, 1));
  ASSERT_TRUE(!MakeParameters(0, 5, 3));
  ASSERT_TRUE(!MakeParameters(1, 1, 0));
  ASSERT_TRUE(!MakeParameters(1, 1, 15));
  ASSERT_TRUE(!MakeParameters(1, 1, 3, MaxRingBits));
  ASSERT_TRUE(!SelectParameters(std::numeric_limits<Int>::max(), 1));
  bool rejected = false;
  try
  {
    MultiplyFermat(Limbs{1}, Limbs{1}, Parameters{3, 8, 1, 1});
  }
  catch (const std::invalid_argument &)
  {
    rejected = true;
  }
  ASSERT_TRUE(rejected);
  // About 20 billion decimal digits in each operand.
  const auto large = SelectParameters(1038102530, 1038102530);
  ASSERT_TRUE(large.has_value());
  ASSERT_EQ(large->count, 1024u);
  ASSERT_TRUE(large->ringLimbs * 64 <= MaxRingBits);
  ASSERT_TRUE(!MakeParameters(1038102530, 1038102530, large->exponent - 1));

  const auto a = Limbs(29087, ~DataT{0});
  std::mt19937_64 rng(0x2a136c);
  const auto b = RandomValue(46019, rng);
  const auto params = SelectParameters(a.size(), b.size());
  ASSERT_TRUE(params.has_value());
  ASSERT_TRUE(MultiplyFermat(a, b, *params) == KaratsubaMultiplication::Multiply(a, b, Base2_64));
}

#if BIGMATH_NTT_WIDE && BIGMATH_NTT_TRUNC && BIGMATH_NTT_MFA && BIGMATH_NTT_MFA_FUSE
REGISTER_TEST(FermatOuter, DispatchOptInAboveWideCrt)
{
  std::mt19937_64 rng(0x426413);
  for (BaseT base : {Base2_64, Base2_32})
  {
    auto a = RandomValue(723, rng), b = RandomValue(633, rng);
    if (base == Base2_32)
    {
      for (auto &x : a) x &= 0xffffffffu;
      for (auto &x : b) x &= 0xffffffffu;
    }
    const auto want = ClassicMultiplication::Multiply(a, b, base);
    constexpr Int h = 512;
    const auto blocked = NttCrt::MultiplyAboveTransform(a, b, base, h, 128, 256, 384, 1536);
    const auto optIn = NttCrt::MultiplyAboveTransform(a, b, base, h, 128, 256, 384, 1536, 1);
    ASSERT_TRUE(blocked == want);
    ASSERT_TRUE(optIn == want);
  }
}
#endif
