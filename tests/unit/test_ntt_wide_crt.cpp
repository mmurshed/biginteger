#include "unit_test_framework.h"

#include <random>

#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"

#if BIGMATH_NTT_TRUNC && BIGMATH_NTT_MFA && BIGMATH_NTT_MFA_FUSE
using namespace BigMath;
using namespace BigMath::NttCrt;
using namespace BigMath::NttCrt::Wide;

namespace
{
  std::vector<DataT> PatternOperand(size_t coeffs, BaseT base, int pattern, std::mt19937_64 &random)
  {
    const size_t perLimb = CoeffsPerLimb(base);
    const DataT mask = base == Base2_64 ? ~DataT{0} : 0xffffffffu;
    std::vector<DataT> v((coeffs + perLimb - 1) / perLimb);
    for (auto &x : v)
      x = pattern == 1 ? mask : (pattern == 2 && random() % 4 ? 0 : random() & mask);
    if (v.back() == 0) v.back() = 1;
    return v;
  }

  template <typename F, UInt G>
  void CheckGenerator(std::initializer_list<UInt> factors)
  {
    for (UInt f : factors)
      ASSERT_TRUE(F::Power(G, (F::Prime - 1) / f) != 1);
  }
}

REGISTER_TEST(NttWideCrt, FieldTransformsAndGenerators)
{
  CheckGenerator<FW1, GW1>({2, 3, 5});
  CheckGenerator<FW2, GW2>({2, 3});
  CheckGenerator<FW3, GW3>({2, 3, 7});

  std::mt19937_64 random(0xa407);
  for (Int n = 2; n <= 1 << 16; n *= 2)
  {
    const auto &plan = GetPlan<FW3, GW3>(n);
    std::vector<UInt> input(n);
    for (UInt &x : input) x = (UInt)(random() % W3);
    for (Int i = 0; i < std::min(n, 8); ++i) input[i] = i % 2 ? W3 - 1 : 0;
    auto reference = input, got = input;
    for (Int length = n; length >= 2; length /= 2)
      ForwardRadix2Layer<FW3>(reference.data(), n, length, plan.forwardRoots.data());
    ForwardPtr<FW3>(got.data(), n, plan);
    ASSERT_TRUE(got == reference);

    for (Int length = 2; length <= n; length *= 2)
      InverseRadix2Layer<FW3>(reference.data(), n, length, plan.inverseRoots.data());
    for (UInt &x : reference) x = FW3::Mul(x, plan.invSize);
    InversePtr<FW3>(got.data(), n, plan, /*scale=*/true);
    ASSERT_TRUE(got == reference);
    ASSERT_TRUE(got == input);
  }
}

REGISTER_TEST(NttWideCrt, FactoredWeightsAndFold)
{
  constexpr Int H = 256;
  std::mt19937_64 random(0x5414);
  const UInt w = FW3::Power(GW3, 11), s = FW3::Power(w, H);
  const PowerTable<FW3> table(w, 1 << 16);
  for (Int test = 0; test < 1000; ++test)
  {
    const UInt i = (UInt)(random() % (1 << 16));
    ASSERT_TRUE(table.At(i) == FW3::Power(w, i));
  }

  for (BaseT base : {Base2_32, Base2_64})
    for (int pattern : {0, 1, 2})
    {
      const auto a = PatternOperand(5 * H - 6, base, pattern, random);
      const BiasedOperandView src(PackedOperandView(a, base));
      std::vector<UInt> plane(H);
      FoldTwist(src, H, s, table, plane.data());
      for (SizeT i = 0; i < H; ++i)
      {
        UInt want = 0;
        for (SizeT u = i; u < src.coeffCount; u += H)
          want = FW3::Add(want, FW3::Mul(FW3::Power(w, u), src.Coeff<FW3>(u)));
        ASSERT_EQ(plane[i], want);
      }
    }
}

REGISTER_TEST(NttWideCrt, ProductsMatchIndependentMultiplication)
{
  constexpr Int H = 512;
  std::mt19937_64 random(0x5733);
  for (BaseT base : {Base2_32, Base2_64})
  {
    const size_t perLimb = CoeffsPerLimb(base);
    auto reference = [base](const auto &a, const auto &b) {
      return a.size() + b.size() <= 2048 ? ClassicMultiplication::Multiply(a, b, base)
                                       : KaratsubaMultiplication::Multiply(a, b, base);
    };
    for (Int k : {1, 2, 3, 7, 16, 31, 53})
      for (Int tail : {2, 63, 64, 65, H / 2, H / 2 + 1, H})
      {
        // On Base2_64, C is odd. Use the adjacent valid tail as needed.
        const size_t sum = ((k * H + tail + 1) / perLimb) * perLimb;
        for (int pattern : {0, 1, 2})
        {
          for (size_t divisor : {2, 6})
          {
            const size_t left = (sum / divisor / perLimb) * perLimb;
            const auto a = PatternOperand(left, base, pattern, random);
            const auto b = PatternOperand(sum - left, base, pattern, random);
            const auto want = reference(a, b);
            ASSERT_TRUE(MultiplyWide(a, b, base, H) == want);
            ASSERT_TRUE(MultiplyWide(b, a, base, H) == want);
          }
          const size_t squareSize = ((sum / 2 + perLimb - 1) / perLimb) * perLimb;
          const auto square = PatternOperand(squareSize, base, pattern, random);
          const ULong count = 2 * squareSize - 1;
          if (count <= (WideMaxResidues + 1) * H)
            ASSERT_TRUE(MultiplyWide(square, square, base, H) == reference(square, square));
        }
      }
  }
}

REGISTER_TEST(NttWideCrt, FusedTransformsMatchKaratsuba)
{
  constexpr Int H = 1 << 15;
  if (!MfaFusedWindow(H)) return;
  std::mt19937_64 random(0x3864);
  for (BaseT base : {Base2_32, Base2_64})
    for (const auto &[na, nb] : {std::pair{40000, 63000}, std::pair{100000, 309600}})
    {
      const auto a = PatternOperand(na, base, 0, random);
      const auto b = PatternOperand(nb, base, 1, random);
      ASSERT_TRUE(MultiplyWide(a, b, base, H) == KaratsubaMultiplication::Multiply(a, b, base));
    }
}

REGISTER_TEST(NttWideCrt, CenteredGarnerAndCarry)
{
  constexpr ULong128 half = (WideModulus - 1) / 2;
  for (ULong128 x : {ULong128{0}, half - 1, half, half + 1, half + 2, WideModulus - 1})
    ASSERT_TRUE(WideCrt::Value(x % W1, x % W2, x % W3, WideCrt::GetInverses()) == x);

  for (BaseT base : {Base2_32, Base2_64})
  {
    const DataT mask = base == Base2_64 ? ~DataT{0} : 0xffffffffu;
    const std::vector<DataT> a(4096, mask), b(3, mask);
    const PackedOperandView va(a, base), vb(b, base);
    const SizeT count = (SizeT)(va.coeffCount + vb.coeffCount - 1);
    const ULong128 shiftedProduct = (((ULong128)1 << 31) - 1) * (((ULong128)1 << 31) - 1);
    const ULong128 rawProduct = (ULong128)0xffffffffu * 0xffffffffu;
    std::vector<UInt> r1(count), r2(count), r3(count);
    // Accumulate each synthetic coefficient in an independent radix-2^32
    // result. The positive values near half exercise carries across chunks.
    std::vector<UInt> digits(count + 4);
    for (SizeT i = 0; i < count; ++i)
    {
      const ULong terms = std::min<ULong>({(ULong)i + 1, vb.coeffCount, (ULong)count - i});
      const bool positive = i % 7 != 0;
      const ULong128 x = positive ? half - i : WideModulus - terms * shiftedProduct;
      r1[i] = x % W1;
      r2[i] = x % W2;
      r3[i] = x % W3;
      ULong128 c = positive ? x + terms * (rawProduct - shiftedProduct)
                            : terms * (rawProduct - 2 * shiftedProduct);
      for (SizeT j = i; c != 0; ++j)
      {
        c += digits[j];
        digits[j] = (UInt)c;
        c >>= 32;
      }
    }
    std::vector<DataT> want;
    const SizeT perLimb = CoeffsPerLimb(base);
    for (SizeT i = 0; i < digits.size(); i += perLimb)
    {
      DataT value = digits[i];
      if (perLimb == 2 && i + 1 < digits.size()) value |= (DataT)digits[i + 1] << 32;
      want.push_back(value);
    }
    TrimZeros(want);
    ASSERT_TRUE(FinalizeBiasedProduct<WideCrt>(r1, r2, r3, va, vb, base, 0) == want);
  }
}

#if BIGMATH_NTT_WIDE
REGISTER_TEST(NttWideCrt, DispatchesStandardWideAndBlockedProducts)
{
  constexpr Int H = 512;
  constexpr ULong plain = 400, biased = 700, wideTerms = 1100, wideCount = 6500;
  std::mt19937_64 random(0x4353);
  for (BaseT base : {Base2_32, Base2_64})
    for (const auto &[na, nb] : {std::pair{600, 640}, std::pair{800, 1400}, std::pair{600, 4100},
                                 std::pair{1800, 3200}, std::pair{1000, 9000}, std::pair{1800, 9000}})
    {
      const auto a = PatternOperand(na, base, 0, random);
      const auto b = PatternOperand(nb, base, 1, random);
      const auto want = ClassicMultiplication::Multiply(a, b, base);
      auto product = [base](const auto &left, const auto &right) {
        return MultiplyAboveTransform(left, right, base, H, plain, biased, wideTerms, wideCount);
      };
      ASSERT_TRUE(product(a, b) == want);
      ASSERT_TRUE(product(b, a) == want);
      ASSERT_TRUE(product(a, a) == ClassicMultiplication::Multiply(a, a, base));
    }
}
#endif
#endif
