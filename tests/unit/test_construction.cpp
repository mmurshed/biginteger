// Construction, accessors, sign normalization, builder.

#include "unit_test_framework.h"

#include <limits>
#include <string>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/Comparison.h"

using namespace BigMath;

REGISTER_TEST(Construction, DefaultIsZero)
{
  BigInteger x;
  ASSERT_TRUE(x.Zero());
  ASSERT_FALSE(x.IsNegative());
  ASSERT_EQ(x.size(), 1u);
  ASSERT_EQ(x[0], 0u);
}

REGISTER_TEST(Construction, SizedDefaultZero)
{
  BigInteger x(4, false);
  ASSERT_TRUE(x.Zero());
  ASSERT_EQ(x.size(), 4u);
  for (SizeT i = 0; i < x.size(); ++i)
    ASSERT_EQ(x[i], 0u);
}

REGISTER_TEST(Construction, NegativeZeroNormalized)
{
  // Constructing a zero BigInteger with negative=true must canonicalize to non-negative.
  BigInteger x(1, true);
  ASSERT_TRUE(x.Zero());
  ASSERT_FALSE(x.IsNegative());

  std::vector<DataT> zero{0, 0, 0};
  BigInteger y(zero, true);
  ASSERT_TRUE(y.Zero());
  ASSERT_FALSE(y.IsNegative());
}

REGISTER_TEST(Construction, FromVectorTrimmed)
{
  // Leading-zero limbs are trimmed; value preserved.
  std::vector<DataT> v{42, 0, 0, 0};
  BigInteger x(v, false);
  // Trimmed to one limb of value 42.
  ASSERT_EQ(x.size(), 1u);
  ASSERT_EQ(x[0], 42u);

  std::vector<DataT> w{1, 2, 0};
  BigInteger y(w, true);
  ASSERT_EQ(y.size(), 2u);
  ASSERT_EQ(y[0], 1u);
  ASSERT_EQ(y[1], 2u);
  ASSERT_TRUE(y.IsNegative());
}

REGISTER_TEST(Construction, SizedFill)
{
  BigInteger x(3, false, 7);
  ASSERT_FALSE(x.Zero());
  ASSERT_EQ(x.size(), 3u);
  ASSERT_EQ(x[0], 7u);
  ASSERT_EQ(x[1], 7u);
  ASSERT_EQ(x[2], 7u);
}

REGISTER_TEST(Construction, CopyAndAssign)
{
  BigInteger a = BigIntegerBuilder::From("-12345678901234567890");
  BigInteger b(a);
  ASSERT_EQ(a, b);
  ASSERT_EQ(a.IsNegative(), b.IsNegative());

  BigInteger c;
  c = a;
  ASSERT_EQ(a, c);

  // Self-assign safe.
  c = c;
  ASSERT_EQ(a, c);
}

REGISTER_TEST(Construction, UnaryNegationFlipsSign)
{
  BigInteger a = BigIntegerBuilder::From("123456789012345678");
  ASSERT_FALSE(a.IsNegative());
  BigInteger b = -a;
  ASSERT_TRUE(b.IsNegative());
  // Operand is unchanged: operator-() is value-returning, not in-place.
  ASSERT_FALSE(a.IsNegative());
  BigInteger c = -b;
  ASSERT_FALSE(c.IsNegative());
  ASSERT_EQ(a, c);
}

REGISTER_TEST(Construction, UnaryNegationZeroStaysPositive)
{
  BigInteger z;
  BigInteger nz = -z;
  ASSERT_FALSE(nz.IsNegative());
  ASSERT_TRUE(nz.Zero());
}

REGISTER_TEST(Construction, SetSignZeroIgnored)
{
  BigInteger z;
  z.SetSign(true);
  ASSERT_FALSE(z.IsNegative());
}

REGISTER_TEST(Construction, SetSignNonZero)
{
  BigInteger a = BigIntegerBuilder::From("42");
  a.SetSign(true);
  ASSERT_TRUE(a.IsNegative());
  a.SetSign(false);
  ASSERT_FALSE(a.IsNegative());
}

REGISTER_TEST(Builder, FromULongZero)
{
  BigInteger x = BigIntegerBuilder::From((ULong)0);
  ASSERT_TRUE(x.Zero());
  ASSERT_FALSE(x.IsNegative());
}

REGISTER_TEST(Builder, FromULongMax)
{
  ULong v = std::numeric_limits<ULong>::max();
  BigInteger x = BigIntegerBuilder::From(v);
  ASSERT_FALSE(x.IsNegative());
  ASSERT_EQ(ToString(x), std::to_string(v));
}

REGISTER_TEST(Builder, FromLongNegative)
{
  Long v = -123456789012345LL;
  BigInteger x = BigIntegerBuilder::From(v);
  ASSERT_TRUE(x.IsNegative());
  ASSERT_EQ(ToString(x), std::to_string(v));
}

REGISTER_TEST(Builder, FromString)
{
  BigInteger x = BigIntegerBuilder::From(std::string("100000000000000000000"));
  ASSERT_FALSE(x.IsNegative());
  ASSERT_EQ(ToString(x), "100000000000000000000");
}

REGISTER_TEST(Builder, FromCStrLeadingPlus)
{
  BigInteger x = BigIntegerBuilder::From("+999");
  ASSERT_FALSE(x.IsNegative());
  ASSERT_EQ(ToString(x), "999");
}

REGISTER_TEST(Builder, FromCStrLeadingZeros)
{
  BigInteger x = BigIntegerBuilder::From("00012345");
  ASSERT_EQ(ToString(x), "12345");
}

REGISTER_TEST(Builder, FromLongDoubleZero)
{
  BigInteger x = BigIntegerBuilder::From((long double)0.0L);
  ASSERT_TRUE(x.Zero());
  ASSERT_FALSE(x.IsNegative());
}

REGISTER_TEST(Builder, FromLongDoublePositive)
{
  BigInteger x = BigIntegerBuilder::From((long double)4294967296.0L);  // = 2^32
  ASSERT_EQ(ToString(x), "4294967296");
}

REGISTER_TEST(Builder, FromLongDoubleNegative)
{
  BigInteger x = BigIntegerBuilder::From((long double)-1000000.0L);
  ASSERT_TRUE(x.IsNegative());
  ASSERT_EQ(ToString(x), "-1000000");
}

#if !BIGMATH_LIMB_64
// Pins the Base2_32 limb encoding of 2^32. Under LIMB_64=1 the same value
// fits in one 64-bit limb (v = {2^32}); skip rather than assert wrong layout.
REGISTER_TEST(Builder, VectorFromULong)
{
  std::vector<DataT> v = BigIntegerBuilder::VectorFrom((ULong)Base2_32);  // 2^32
  // 2^32 = 1 * 2^32 + 0  →  limbs {0, 1}
  ASSERT_EQ(v.size(), 2u);
  ASSERT_EQ(v[0], 0u);
  ASSERT_EQ(v[1], 1u);
}
#endif

REGISTER_TEST(Construction, MoveLimbVector)
{
  std::vector<DataT> limbs = {12345, 67890};
  BigInteger x(std::move(limbs), true);
  ASSERT_TRUE(x.IsNegative());
  ASSERT_EQ(x.size(), 2u);
  ASSERT_EQ(x[0], 12345u);
  ASSERT_EQ(x[1], 67890u);
  // No assertion on `limbs` here: a moved-from vector is in a valid but
  // unspecified state per the standard.
}

REGISTER_TEST(Construction, ReleaseInteger)
{
  BigInteger x(std::vector<DataT>{100, 200}, true);
  std::vector<DataT> limbs = x.ReleaseInteger();
  ASSERT_EQ(limbs.size(), 2u);
  ASSERT_EQ(limbs[0], 100u);
  ASSERT_EQ(limbs[1], 200u);
  
  // x should be left in a valid positive zero state
  ASSERT_TRUE(x.Zero());
  ASSERT_FALSE(x.IsNegative());
  ASSERT_EQ(x.size(), 1u);
  ASSERT_EQ(x[0], 0u);
}

REGISTER_TEST(Serialization, ByteArrays)
{
  using ByteOrder = BigInteger::ByteOrder;

  // Test value 0
  {
    BigInteger z;
    auto bytes = z.ToByteArray(ByteOrder::BigEndian);
    ASSERT_TRUE(bytes.empty());

    BigInteger z2 = BigInteger::FromByteArray(bytes, false, ByteOrder::BigEndian);
    ASSERT_TRUE(z2.Zero());
    ASSERT_FALSE(z2.IsNegative());
  }

  // Test standard value, big-endian and little-endian
  {
    // 0x1234567890abcdef
    ULong val = 0x1234567890abcdefULL;
    BigInteger x = BigIntegerBuilder::From(val);

    // Big Endian bytes
    std::vector<uint8_t> expectedBE = {
      0x12, 0x34, 0x56, 0x78, 0x90, 0xab, 0xcd, 0xef
    };
    auto bytesBE = x.ToByteArray(ByteOrder::BigEndian);
    ASSERT_EQ(bytesBE.size(), expectedBE.size());
    for (size_t i = 0; i < bytesBE.size(); ++i) {
      ASSERT_EQ(bytesBE[i], expectedBE[i]);
    }

    BigInteger x2 = BigInteger::FromByteArray(bytesBE, true, ByteOrder::BigEndian);
    ASSERT_TRUE(x2.IsNegative());
    ASSERT_EQ(x2, -x);

    // Little Endian bytes
    std::vector<uint8_t> expectedLE = {
      0xef, 0xcd, 0xab, 0x90, 0x78, 0x56, 0x34, 0x12
    };
    auto bytesLE = x.ToByteArray(ByteOrder::LittleEndian);
    ASSERT_EQ(bytesLE.size(), expectedLE.size());
    for (size_t i = 0; i < bytesLE.size(); ++i) {
      ASSERT_EQ(bytesLE[i], expectedLE[i]);
    }

    BigInteger x3 = BigInteger::FromByteArray(bytesLE, false, ByteOrder::LittleEndian);
    ASSERT_FALSE(x3.IsNegative());
    ASSERT_EQ(x3, x);
  }
}

REGISTER_TEST(Serialization, OddLengthByteArrays)
{
  using ByteOrder = BigInteger::ByteOrder;

  // 5 bytes: partial top limb in both limb configurations.
  {
    std::vector<uint8_t> be = {0x01, 0x02, 0x03, 0x04, 0x05};
    BigInteger x = BigInteger::FromByteArray(be, false);
    ASSERT_EQ(x, BigIntegerBuilder::From(0x0102030405ULL));

    auto rt = x.ToByteArray();
    ASSERT_EQ(rt.size(), be.size());
    for (size_t i = 0; i < rt.size(); ++i)
      ASSERT_EQ(rt[i], be[i]);
  }

  // 9 bytes: partial top limb above a full 64-bit limb.
  {
    std::vector<uint8_t> be = {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09};
    BigInteger x = BigInteger::FromByteArray(be, false);
    // 0x010203040506070809
    ASSERT_EQ(ToString(x), std::string("18591708106338011145"));

    auto rt = x.ToByteArray();
    ASSERT_EQ(rt.size(), be.size());
    for (size_t i = 0; i < rt.size(); ++i)
      ASSERT_EQ(rt[i], be[i]);

    BigInteger x2 = BigInteger::FromByteArray(x.ToByteArray(ByteOrder::LittleEndian), false, ByteOrder::LittleEndian);
    ASSERT_EQ(x2, x);
  }
}

REGISTER_TEST(Serialization, MultiLimbRoundTrip)
{
  using ByteOrder = BigInteger::ByteOrder;

  // 17 bytes: 3 limbs at 64-bit limbs, 5 limbs at 32-bit limbs; the top limb
  // is partial in both configurations.
  std::vector<uint8_t> be;
  for (uint8_t v = 1; v <= 17; ++v)
    be.push_back(v);

  BigInteger x = BigInteger::FromByteArray(be, false);
  // 0x0102030405060708090a0b0c0d0e0f1011
  ASSERT_EQ(ToString(x), std::string("342956481330728537355412814650493833233"));
  ASSERT_EQ(x.size(), (be.size() + LimbBits / 8 - 1) / (LimbBits / 8));

  auto rtBE = x.ToByteArray(ByteOrder::BigEndian);
  ASSERT_EQ(rtBE.size(), be.size());
  for (size_t i = 0; i < rtBE.size(); ++i)
    ASSERT_EQ(rtBE[i], be[i]);

  auto le = x.ToByteArray(ByteOrder::LittleEndian);
  ASSERT_EQ(le.size(), be.size());
  for (size_t i = 0; i < le.size(); ++i)
    ASSERT_EQ(le[i], be[be.size() - 1 - i]);

  BigInteger x2 = BigInteger::FromByteArray(le, false, ByteOrder::LittleEndian);
  ASSERT_EQ(x2, x);
}

REGISTER_TEST(Serialization, LeadingZeroBytes)
{
  using ByteOrder = BigInteger::ByteOrder;

  BigInteger expected = BigIntegerBuilder::From(0x0102ULL);

  // Big-endian input with leading zero bytes parses to the same value.
  std::vector<uint8_t> paddedBE = {0x00, 0x00, 0x00, 0x01, 0x02};
  BigInteger x = BigInteger::FromByteArray(paddedBE, false, ByteOrder::BigEndian);
  ASSERT_EQ(x, expected);

  // Little-endian input with trailing zero bytes parses to the same value.
  std::vector<uint8_t> paddedLE = {0x02, 0x01, 0x00, 0x00, 0x00};
  BigInteger y = BigInteger::FromByteArray(paddedLE, false, ByteOrder::LittleEndian);
  ASSERT_EQ(y, expected);

  // Serialization is canonical: zero padding never round-trips.
  auto bytes = x.ToByteArray(ByteOrder::BigEndian);
  ASSERT_EQ(bytes.size(), 2u);
  ASSERT_EQ(bytes[0], 0x01u);
  ASSERT_EQ(bytes[1], 0x02u);
}
