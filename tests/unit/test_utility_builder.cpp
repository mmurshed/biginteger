// Limb-vector utilities (TrimZeros, IsZero, Compare, gcd) plus
// BigIntegerBuilder long-double conversion.

#include "unit_test_framework.h"

#include <limits>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/common/Parser.h"
#include "biginteger/common/Util.h"
#include "biginteger/ops/Comparison.h"

using namespace BigMath;

// ─── gcd ─────────────────────────────────────────────────────────────────────

REGISTER_TEST(Util, GcdBasic)
{
  ASSERT_EQ(gcd((DataT)0,  (DataT)5), 5u);
  ASSERT_EQ(gcd((DataT)5,  (DataT)0), 5u);
  ASSERT_EQ(gcd((DataT)6,  (DataT)9), 3u);
  ASSERT_EQ(gcd((DataT)17, (DataT)5), 1u);
  ASSERT_EQ(gcd((DataT)100, (DataT)75), 25u);
}

// ─── TrimZeros / TrimZerosToOne ──────────────────────────────────────────────

REGISTER_TEST(Util, TrimZerosStripsTrailing)
{
  std::vector<DataT> v{1, 2, 3, 0, 0, 0};
  SizeT removed = TrimZeros(v);
  ASSERT_EQ(removed, 3u);
  ASSERT_EQ(v.size(), 3u);
  ASSERT_EQ(v[2], 3u);
}

REGISTER_TEST(Util, TrimZerosToOneNeverEmpties)
{
  std::vector<DataT> v{0, 0, 0};
  TrimZerosToOne(v);
  ASSERT_EQ(v.size(), 1u);
  ASSERT_EQ(v[0], 0u);
}

// ─── IsZero ──────────────────────────────────────────────────────────────────

REGISTER_TEST(Util, IsZeroEmpty)
{
  std::vector<DataT> v;
  ASSERT_TRUE(IsZero(v));
}

REGISTER_TEST(Util, IsZeroSingleLimb)
{
  std::vector<DataT> v{0};
  ASSERT_TRUE(IsZero(v));
  v[0] = 7;
  ASSERT_FALSE(IsZero(v));
}

REGISTER_TEST(Util, IsZeroAllZerosManyLimbs)
{
  std::vector<DataT> v{0, 0, 0, 0};
  ASSERT_TRUE(IsZero(v));
}

REGISTER_TEST(Util, IsZeroAnyNonZero)
{
  std::vector<DataT> v{0, 0, 0, 1};
  ASSERT_FALSE(IsZero(v));
  std::vector<DataT> w{0, 1, 0, 0};
  ASSERT_FALSE(IsZero(w));
}

// ─── Compare ─────────────────────────────────────────────────────────────────

REGISTER_TEST(Util, CompareEqual)
{
  std::vector<DataT> a{1, 2, 3};
  std::vector<DataT> b{1, 2, 3};
  ASSERT_EQ(Compare(a, b), 0);
}

REGISTER_TEST(Util, CompareIgnoresTrailingZeros)
{
  std::vector<DataT> a{1, 2, 3};
  std::vector<DataT> b{1, 2, 3, 0, 0};
  ASSERT_EQ(Compare(a, b), 0);
}

REGISTER_TEST(Util, CompareHighLimbWins)
{
  std::vector<DataT> a{0xFFFFFFFF, 0xFFFFFFFF, 1};
  std::vector<DataT> b{0,          0,          2};
  ASSERT_LT(Compare(a, b), 0);
}

// ─── Builder edge cases ──────────────────────────────────────────────────────

REGISTER_TEST(Builder, LongDoubleLargeRoundTrip)
{
  long double v = 1.0e18L;
  BigInteger x = BigIntegerBuilder::From(v);
  // 1e18 fits exactly in 60 bits → exact decimal representation expected.
  ASSERT_EQ(ToString(x), "1000000000000000000");
}

REGISTER_TEST(Builder, ULongRoundTrip)
{
  for (ULong v : {(ULong)0, (ULong)1, (ULong)42,
                  (ULong)0xFFFFFFFFULL,                  // 2^32 - 1
                  (ULong)0x100000000ULL,                 // 2^32
                  (ULong)0xFFFFFFFFFFFFFFFFULL})         // 2^64 - 1
  {
    BigInteger x = BigIntegerBuilder::From(v);
    ASSERT_EQ(ToString(x), std::to_string(v));
  }
}

REGISTER_TEST(Builder, LongRoundTrip)
{
  for (Long v : {(Long)0, (Long)1, (Long)-1, (Long)42, (Long)-9999999999LL,
                 std::numeric_limits<Long>::min() + 1,   // skip INT64_MIN: |x| not representable in Long
                 std::numeric_limits<Long>::max()})
  {
    BigInteger x = BigIntegerBuilder::From(v);
    ASSERT_EQ(ToString(x), std::to_string(v));
  }
}
