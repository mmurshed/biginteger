// Direct exercise of Base2_64 paths in ClassicDivision + FastDivision.
// BZ and Newton fall back to FastDivision for Base2_64, so testing
// FastDivision directly covers all three dispatched paths.
//
// Cross-checks the quotient/remainder identity `q*b + r == a` and `r < b`
// without relying on a reference Base2_64 divider (none exists yet).

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/division/ClassicDivision.h"
#include "biginteger/algorithms/division/FastDivision.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
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

static std::vector<DataT> Trim(std::vector<DataT> v)
{
  while (v.size() > 1 && v.back() == 0)
    v.pop_back();
  return v;
}

static bool VerifyIdentity(std::vector<DataT> const &a, std::vector<DataT> const &b,
                            std::vector<DataT> const &q, std::vector<DataT> const &r)
{
  // recomposed = q * b + r
  std::vector<DataT> qb(q.size() + b.size() + 1, 0);
  ClassicMultiplication::Multiply(q, 0, (SizeT)q.size() - 1,
                                  b, 0, (SizeT)b.size() - 1,
                                  qb, 0, Base2_64);
  qb = Trim(qb);
  std::vector<DataT> recomposed = Add(qb, r, Base2_64);
  return Compare(Trim(recomposed), Trim(a)) == 0 && Compare(r, b) < 0;
}

REGISTER_TEST(Base64Div, ClassicScalarSmall)
{
  std::vector<DataT> u{100, 200, 300};
  auto qr = ClassicDivision::DivideAndRemainder(u, (DataT)7, Base2_64);
  ASSERT_TRUE(VerifyIdentity(u, {(DataT)7}, qr.first, qr.second));
}

REGISTER_TEST(Base64Div, ClassicScalarMaxDivisor)
{
  std::vector<DataT> u{0xDEADBEEFCAFEBABEULL, 0x0123456789ABCDEFULL, 0x42};
  auto qr = ClassicDivision::DivideAndRemainder(u, (DataT)0xFFFFFFFFFFFFFFFFULL, Base2_64);
  ASSERT_TRUE(VerifyIdentity(u, {(DataT)0xFFFFFFFFFFFFFFFFULL}, qr.first, qr.second));
}

REGISTER_TEST(Base64Div, FastDivisionSmall2x1)
{
  std::vector<DataT> a{100, 200};
  std::vector<DataT> b{7};
  auto qr = FastDivision::DivideAndRemainder(a, b, Base2_64);
  ASSERT_TRUE(VerifyIdentity(a, b, qr.first, qr.second));
}

REGISTER_TEST(Base64Div, FastDivisionMidBalanced)
{
  std::mt19937_64 gen(0xC0FFEEULL);
  auto a = RandomLimbs64(16, gen);
  auto b = RandomLimbs64(8, gen);
  auto qr = FastDivision::DivideAndRemainder(a, b, Base2_64);
  ASSERT_TRUE(VerifyIdentity(a, b, qr.first, qr.second));
}

REGISTER_TEST(Base64Div, FastDivisionLargerBalanced)
{
  std::mt19937_64 gen(0xBADF00DULL);
  auto a = RandomLimbs64(128, gen);
  auto b = RandomLimbs64(64, gen);
  auto qr = FastDivision::DivideAndRemainder(a, b, Base2_64);
  ASSERT_TRUE(VerifyIdentity(a, b, qr.first, qr.second));
}

REGISTER_TEST(Base64Div, FastDivisionSkewed)
{
  std::mt19937_64 gen(0xDEADBEEFULL);
  auto a = RandomLimbs64(256, gen);
  auto b = RandomLimbs64(32, gen);
  auto qr = FastDivision::DivideAndRemainder(a, b, Base2_64);
  ASSERT_TRUE(VerifyIdentity(a, b, qr.first, qr.second));
}

REGISTER_TEST(Base64Div, FastDivisionMaxFillDivisor)
{
  // All-max-limb divisor — exercises the normalize d=1 branch in the new
  // Base2_64 path.
  std::vector<DataT> a{0x1, 0x2, 0x3, 0x4, 0x5, 0x6};
  std::vector<DataT> b(3, 0xFFFFFFFFFFFFFFFFULL);
  auto qr = FastDivision::DivideAndRemainder(a, b, Base2_64);
  ASSERT_TRUE(VerifyIdentity(a, b, qr.first, qr.second));
}

REGISTER_TEST(Base64Div, FastDivisionLeadingFFDivisor)
{
  // Adversarial qhat-correction case: divisor starts with 0xFF..FF in top limb.
  std::mt19937_64 gen(0xABCDEFULL);
  auto a = RandomLimbs64(96, gen);
  std::vector<DataT> b(32, 0xFFFFFFFFFFFFFFFFULL);
  auto qr = FastDivision::DivideAndRemainder(a, b, Base2_64);
  ASSERT_TRUE(VerifyIdentity(a, b, qr.first, qr.second));
}

REGISTER_TEST(Base64Div, FastDivisionAlessB)
{
  std::vector<DataT> a{0x123};
  std::vector<DataT> b{0xFFFF, 0xABC};
  auto qr = FastDivision::DivideAndRemainder(a, b, Base2_64);
  ASSERT_EQ((ULong)0, (ULong)qr.first[0]);
  ASSERT_EQ((SizeT)1, (SizeT)qr.first.size());
  ASSERT_EQ((ULong)0x123, (ULong)qr.second[0]);
}
