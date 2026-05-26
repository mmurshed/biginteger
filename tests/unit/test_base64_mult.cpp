// Direct exercise of the Base2_64 sentinel paths in ClassicMultiplication.
// Dormant under LIMB_64=0 default; covered here via explicit Base2_64 calls.

#include "unit_test_framework.h"

#include <cstdint>
#include <vector>

#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/common/Constants.h"

using namespace BigMath;

REGISTER_TEST(Base64Mult, ScalarSmall)
{
  // (3 + 5*B) * 7 = 21 + 35*B.
  std::vector<DataT> a{3, 5};
  ClassicMultiplication::MultiplyTo(a, 0, (SizeT)a.size() - 1, (DataT)7, Base2_64);
  ASSERT_EQ((ULong)21, (ULong)a[0]);
  ASSERT_EQ((ULong)35, (ULong)a[1]);
}

REGISTER_TEST(Base64Mult, ScalarCarryToNewLimb)
{
  // (2^64 - 1) * 2 = 2*(2^64-1) = 0xFFFFFFFFFFFFFFFE + (1 << 64).
  std::vector<DataT> a{0xFFFFFFFFFFFFFFFFULL};
  ClassicMultiplication::MultiplyTo(a, 0, (SizeT)a.size() - 1, (DataT)2, Base2_64);
  ASSERT_EQ((ULong)0xFFFFFFFFFFFFFFFEULL, (ULong)a[0]);
  ASSERT_EQ((SizeT)2, (SizeT)a.size());
  ASSERT_EQ((ULong)1, (ULong)a[1]);
}

REGISTER_TEST(Base64Mult, ScalarMaxByMax)
{
  // (2^64 - 1) * (2^64 - 1) = 2^128 - 2^65 + 1
  //                         = (1, max-1) in (high, low) base-2^64 = limbs[1]=max-1, limbs[0]=1.
  std::vector<DataT> a{0xFFFFFFFFFFFFFFFFULL};
  ClassicMultiplication::MultiplyTo(a, 0, (SizeT)a.size() - 1, (DataT)0xFFFFFFFFFFFFFFFFULL, Base2_64);
  ASSERT_EQ((ULong)1, (ULong)a[0]);
  ASSERT_EQ((SizeT)2, (SizeT)a.size());
  ASSERT_EQ((ULong)0xFFFFFFFFFFFFFFFEULL, (ULong)a[1]);
}

REGISTER_TEST(Base64Mult, VectorScalarRangeForm)
{
  // a = (10, 20, 30), b = 100, into fresh w.
  std::vector<DataT> a{10, 20, 30};
  std::vector<DataT> w(4, 0);
  ClassicMultiplication::Multiply(a, 0, 2, (ULong)100, w, 0, 3, Base2_64);
  ASSERT_EQ((ULong)1000, (ULong)w[0]);
  ASSERT_EQ((ULong)2000, (ULong)w[1]);
  ASSERT_EQ((ULong)3000, (ULong)w[2]);
  ASSERT_EQ((ULong)0, (ULong)w[3]);
}

REGISTER_TEST(Base64Mult, FullProductSmall)
{
  // (3 + 5*B) * (7 + 11*B) = 21 + 33*B + 35*B + 55*B^2 = 21 + 68*B + 55*B^2.
  std::vector<DataT> a{3, 5};
  std::vector<DataT> b{7, 11};
  std::vector<DataT> r(5, 0);
  ClassicMultiplication::Multiply(a, 0, 1, b, 0, 1, r, 0, Base2_64);
  ASSERT_EQ((ULong)21, (ULong)r[0]);
  ASSERT_EQ((ULong)68, (ULong)r[1]);
  ASSERT_EQ((ULong)55, (ULong)r[2]);
}

REGISTER_TEST(Base64Mult, FullProductCarryChain)
{
  // (max, max) * (1, 0) — single-limb b — should give (max, max, 0, 0) shifted.
  // Actually with b = (1, 0): result = a * 1 = (max, max, 0, 0).
  std::vector<DataT> a{0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL};
  std::vector<DataT> b{1, 0};
  std::vector<DataT> r(5, 0);
  ClassicMultiplication::Multiply(a, 0, 1, b, 0, 1, r, 0, Base2_64);
  ASSERT_EQ((ULong)0xFFFFFFFFFFFFFFFFULL, (ULong)r[0]);
  ASSERT_EQ((ULong)0xFFFFFFFFFFFFFFFFULL, (ULong)r[1]);
  ASSERT_EQ((ULong)0, (ULong)r[2]);
  ASSERT_EQ((ULong)0, (ULong)r[3]);
}

REGISTER_TEST(Base64Mult, FullProductMaxByMax2Limb)
{
  // (max, max) * (max, max) = (2^128 - 1)² = 2^256 - 2^129 + 1.
  // Decomposing in base 2^64: low=1, then -2 (i.e. max-1 with borrow), then max, then max.
  // Concretely: 2^256 - 2^129 + 1 = (..., r3, r2, r1, r0).
  //   r0 = 1
  //   r1 = 0
  //   r2 = max - 1   (= 2^64 - 2)
  //   r3 = max       (= 2^64 - 1)
  // Verify by checking sum: r0 + r1*B + r2*B^2 + r3*B^3
  //   = 1 + 0 + (B-2)*B^2 + (B-1)*B^3
  //   = 1 + B^3 - 2*B^2 + B^4 - B^3
  //   = 1 - 2*B^2 + B^4
  //   = B^4 - 2*B^2 + 1
  //   = (B^2 - 1)² = (2^128 - 1)². Matches.
  std::vector<DataT> a{0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL};
  std::vector<DataT> b{0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL};
  std::vector<DataT> r(5, 0);
  ClassicMultiplication::Multiply(a, 0, 1, b, 0, 1, r, 0, Base2_64);
  ASSERT_EQ((ULong)1, (ULong)r[0]);
  ASSERT_EQ((ULong)0, (ULong)r[1]);
  ASSERT_EQ((ULong)0xFFFFFFFFFFFFFFFEULL, (ULong)r[2]);
  ASSERT_EQ((ULong)0xFFFFFFFFFFFFFFFFULL, (ULong)r[3]);
}

REGISTER_TEST(Base64Mult, ScalarZeroIsZero)
{
  std::vector<DataT> a{0xDEADBEEFCAFEBABEULL, 0x1234};
  ClassicMultiplication::MultiplyTo(a, 0, (SizeT)a.size() - 1, (DataT)0, Base2_64);
  // Implementation sets the range to a single zero limb when b == 0.
  ASSERT_EQ((ULong)0, (ULong)a[0]);
}
