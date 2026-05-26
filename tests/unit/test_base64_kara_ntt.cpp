// Direct exercise of Base2_64 sentinel paths in KaratsubaMultiplication and
// NTTMultiplication. Dormant under LIMB_64=0 default; covered here via explicit
// Base2_64 calls. Cross-checked against the schoolbook ClassicMultiplication
// Base2_64 path (covered in test_base64_mult.cpp).

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
#include "biginteger/algorithms/multiplication/NTTMultiplication.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/common/Constants.h"

using namespace BigMath;

static std::vector<DataT> RandomLimbs64(SizeT n, std::mt19937_64 &gen)
{
  std::uniform_int_distribution<uint64_t> dist(0, 0xFFFFFFFFFFFFFFFFULL);
  std::vector<DataT> v(n);
  for (SizeT i = 0; i < n; ++i)
    v[i] = dist(gen);
  if (v.back() == 0)
    v.back() = 1;
  return v;
}

static std::vector<DataT> ClassicProduct(std::vector<DataT> const &a, std::vector<DataT> const &b)
{
  std::vector<DataT> r(a.size() + b.size() + 1, 0);
  ClassicMultiplication::Multiply(a, 0, (SizeT)a.size() - 1, b, 0, (SizeT)b.size() - 1, r, 0, Base2_64);
  while (r.size() > 1 && r.back() == 0)
    r.pop_back();
  return r;
}

static bool LimbVectorsEqual(std::vector<DataT> const &x, std::vector<DataT> const &y)
{
  auto trim = [](std::vector<DataT> v)
  {
    while (v.size() > 1 && v.back() == 0)
      v.pop_back();
    return v;
  };
  auto xt = trim(x);
  auto yt = trim(y);
  if (xt.size() != yt.size())
    return false;
  for (SizeT i = 0; i < xt.size(); ++i)
    if (xt[i] != yt[i])
      return false;
  return true;
}

REGISTER_TEST(Base64Kara, AgainstClassicSmall)
{
  std::mt19937_64 gen(0xC0FFEEULL);
  // 50 limbs > KARATSUBA_THRESHOLD (48), so recursion is exercised.
  auto a = RandomLimbs64(50, gen);
  auto b = RandomLimbs64(50, gen);
  auto k = KaratsubaMultiplication::Multiply(a, b, Base2_64);
  auto c = ClassicProduct(a, b);
  ASSERT_TRUE(LimbVectorsEqual(k, c));
}

REGISTER_TEST(Base64Kara, AgainstClassicMidUnbalanced)
{
  std::mt19937_64 gen(0xBEEFCAFEULL);
  auto a = RandomLimbs64(200, gen);
  auto b = RandomLimbs64(73, gen);
  auto k = KaratsubaMultiplication::Multiply(a, b, Base2_64);
  auto c = ClassicProduct(a, b);
  ASSERT_TRUE(LimbVectorsEqual(k, c));
}

REGISTER_TEST(Base64Kara, AgainstClassicMaxFill)
{
  // Worst-case carry: every limb is 2^64-1.
  std::vector<DataT> a(60, 0xFFFFFFFFFFFFFFFFULL);
  std::vector<DataT> b(60, 0xFFFFFFFFFFFFFFFFULL);
  auto k = KaratsubaMultiplication::Multiply(a, b, Base2_64);
  auto c = ClassicProduct(a, b);
  ASSERT_TRUE(LimbVectorsEqual(k, c));
}

REGISTER_TEST(Base64Ntt, AgainstClassicSmall)
{
  std::mt19937_64 gen(0xDEADBEEFULL);
  auto a = RandomLimbs64(40, gen);
  auto b = RandomLimbs64(40, gen);
  auto n = NTTMultiplication::Multiply(a, b, Base2_64);
  auto c = ClassicProduct(a, b);
  ASSERT_TRUE(LimbVectorsEqual(n, c));
}

REGISTER_TEST(Base64Ntt, AgainstClassicMidUnbalanced)
{
  std::mt19937_64 gen(0xBADDF00DULL);
  auto a = RandomLimbs64(256, gen);
  auto b = RandomLimbs64(97, gen);
  auto n = NTTMultiplication::Multiply(a, b, Base2_64);
  auto c = ClassicProduct(a, b);
  ASSERT_TRUE(LimbVectorsEqual(n, c));
}

REGISTER_TEST(Base64Ntt, AgainstClassicMaxFill)
{
  std::vector<DataT> a(128, 0xFFFFFFFFFFFFFFFFULL);
  std::vector<DataT> b(128, 0xFFFFFFFFFFFFFFFFULL);
  auto n = NTTMultiplication::Multiply(a, b, Base2_64);
  auto c = ClassicProduct(a, b);
  ASSERT_TRUE(LimbVectorsEqual(n, c));
}

REGISTER_TEST(Base64Ntt, KaratsubaVsNttAgree)
{
  std::mt19937_64 gen(0xABCDEFULL);
  auto a = RandomLimbs64(150, gen);
  auto b = RandomLimbs64(150, gen);
  auto k = KaratsubaMultiplication::Multiply(a, b, Base2_64);
  auto n = NTTMultiplication::Multiply(a, b, Base2_64);
  ASSERT_TRUE(LimbVectorsEqual(k, n));
}
