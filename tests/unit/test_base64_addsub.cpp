// Direct exercise of the Base2_64 sentinel paths in Addition / Subtraction.
// These paths are dormant under LIMB_64=0 (BigInteger::Base() still returns
// Base2_32), so without explicit tests the new code is unreachable. Cover the
// branches via direct std::vector<DataT> calls.

#include "unit_test_framework.h"

#include <cstdint>
#include <vector>

#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Subtraction.h"
#include "biginteger/common/Constants.h"

using namespace BigMath;

REGISTER_TEST(Base64Add, ScalarNoCarry)
{
  std::vector<DataT> a{0x100, 0x200};
  AddTo(a, 0, (SizeT)a.size() - 1, 0x50, Base2_64);
  ASSERT_EQ((ULong)0x150, (ULong)a[0]);
  ASSERT_EQ((ULong)0x200, (ULong)a[1]);
}

REGISTER_TEST(Base64Add, ScalarWithCarryPropagate)
{
  // a = 2^64 - 1 (single limb at max) + 1 should propagate carry.
  std::vector<DataT> a{0xFFFFFFFFFFFFFFFFULL};
  AddTo(a, 0, (SizeT)a.size() - 1, 1, Base2_64);
  ASSERT_EQ((ULong)0, (ULong)a[0]);
  ASSERT_EQ((SizeT)2, (SizeT)a.size());
  ASSERT_EQ((ULong)1, (ULong)a[1]);
}

REGISTER_TEST(Base64Add, MultiLimbCarryChain)
{
  // a = (2^64 - 1, 2^64 - 1, 0), b = (1, 0, 0) → (0, 0, 1).
  std::vector<DataT> a{0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL, 0};
  std::vector<DataT> b{1, 0, 0};
  std::vector<DataT> r(4, 0);
  Add(a, 0, 2, b, 0, 2, r, 0, Base2_64);
  ASSERT_EQ((ULong)0, (ULong)r[0]);
  ASSERT_EQ((ULong)0, (ULong)r[1]);
  ASSERT_EQ((ULong)1, (ULong)r[2]);
}

REGISTER_TEST(Base64Add, ScalarIntoEmptyVector)
{
  std::vector<DataT> a;
  AddTo(a, (ULong)0x123456789ABCDEFULL, Base2_64);
  ASSERT_EQ((SizeT)1, (SizeT)a.size());
  ASSERT_EQ((ULong)0x123456789ABCDEFULL, (ULong)a[0]);
}

REGISTER_TEST(Base64Sub, ScalarNoBorrow)
{
  std::vector<DataT> a{0x150, 0x200};
  SubtractFrom(a, 0, (SizeT)a.size() - 1, 0x50, Base2_64);
  ASSERT_EQ((ULong)0x100, (ULong)a[0]);
  ASSERT_EQ((ULong)0x200, (ULong)a[1]);
}

REGISTER_TEST(Base64Sub, ScalarBorrowPropagate)
{
  // a = (0, 1) - 1 = (max, 0) → trimmed to (max).
  std::vector<DataT> a{0, 1};
  SubtractFrom(a, 0, (SizeT)a.size() - 1, 1, Base2_64);
  ASSERT_EQ((SizeT)1, (SizeT)a.size());
  ASSERT_EQ((ULong)0xFFFFFFFFFFFFFFFFULL, (ULong)a[0]);
}

REGISTER_TEST(Base64Sub, MultiLimbBorrowChain)
{
  // a = (0, 0, 1), b = (1, 0, 0) → (max, max, 0).
  std::vector<DataT> a{0, 0, 1};
  std::vector<DataT> b{1, 0, 0};
  std::vector<DataT> r(3, 0);
  Subtract(a, 0, 2, b, 0, 2, r, 0, Base2_64);
  ASSERT_EQ((ULong)0xFFFFFFFFFFFFFFFFULL, (ULong)r[0]);
  ASSERT_EQ((ULong)0xFFFFFFFFFFFFFFFFULL, (ULong)r[1]);
  ASSERT_EQ((ULong)0, (ULong)r[2]);
}

REGISTER_TEST(Base64Sub, EqualOperandsZero)
{
  std::vector<DataT> a{0x123456789ABCDEFULL, 0xFEDCBA9876543210ULL};
  std::vector<DataT> b{0x123456789ABCDEFULL, 0xFEDCBA9876543210ULL};
  std::vector<DataT> r(2, 0);
  Subtract(a, 0, 1, b, 0, 1, r, 0, Base2_64);
  ASSERT_EQ((ULong)0, (ULong)r[0]);
  ASSERT_EQ((ULong)0, (ULong)r[1]);
}

REGISTER_TEST(Base64AddSub, RoundTripRandom)
{
  // For a random pair (a, b), (a + b) - b == a and (a - b) + b == a.
  std::vector<DataT> a{0xDEADBEEFCAFEBABEULL, 0x0123456789ABCDEFULL, 0x42};
  std::vector<DataT> b{0xBADF00DDEADC0DE5ULL, 0xFEDCBA9876543210ULL, 0x07};

  // Compute s = a + b, then s - b should equal a.
  std::vector<DataT> sum(4, 0);
  Add(a, 0, 2, b, 0, 2, sum, 0, Base2_64);
  std::vector<DataT> back(4, 0);
  Subtract(sum, 0, 3, b, 0, 2, back, 0, Base2_64);
  for (SizeT i = 0; i < 3; ++i)
    ASSERT_EQ((ULong)a[i], (ULong)back[i]);
  // back[3] (extra high limb) must be zero.
  ASSERT_EQ((ULong)0, (ULong)back[3]);
}
