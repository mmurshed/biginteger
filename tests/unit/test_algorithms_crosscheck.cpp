// Cross-checks: every multiplication/squaring/division algorithm must agree
// with the canonical reference (Classic for mult/square, BurnikelZiegler for
// division) on the same input. Random seeded inputs across a wide size matrix.

#include "unit_test_framework.h"

#include <random>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Division.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/Squaring.h"
#include "biginteger/algorithms/division/BurnikelZieglerDivision.h"
#include "biginteger/algorithms/division/ClassicDivision.h"
#include "biginteger/algorithms/division/FastDivision.h"
#include "biginteger/algorithms/division/KnuthDivision.h"
#include "biginteger/algorithms/division/NewtonDivision.h"
#include "biginteger/algorithms/division/ReciprocalDivision.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/algorithms/multiplication/ClassicSquare.h"
#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
#include "biginteger/algorithms/multiplication/KaratsubaSquare.h"
#include "biginteger/algorithms/multiplication/NTTMultiplication.h"
#include "biginteger/algorithms/multiplication/NTTSquare.h"
#include "biginteger/algorithms/multiplication/Toom5Multiplication.h"
#include "biginteger/algorithms/multiplication/ToomCookMultiplication.h"
#include "biginteger/common/Comparator.h"

using namespace BigMath;

static std::vector<DataT> RandomLimbs(SizeT limbs, std::mt19937_64 &gen, bool ensureTopNonZero = true)
{
  std::uniform_int_distribution<uint64_t> d(0, 0xFFFFFFFFULL);
  std::vector<DataT> v(limbs);
  for (auto &x : v) x = (DataT)d(gen);
  if (ensureTopNonZero && !v.empty() && v.back() == 0)
    v.back() = 1 + (d(gen) & 0xFFFF);
  return v;
}

// ─── multiplication: 4 algorithms agree ──────────────────────────────────────

static void CrossMul(SizeT aLimbs, SizeT bLimbs, uint64_t seed)
{
  std::mt19937_64 gen(seed);
  auto a = RandomLimbs(aLimbs, gen);
  auto b = RandomLimbs(bLimbs, gen);
  auto classic = ClassicMultiplication::Multiply(a, b, BigInteger::Base());
  auto kara    = KaratsubaMultiplication::Multiply(a, b, BigInteger::Base());
  auto ntt     = NTTMultiplication::Multiply(a, b, BigInteger::Base());
  ASSERT_EQ(Compare(classic, kara), 0);
  ASSERT_EQ(Compare(classic, ntt),  0);
#if !BIGMATH_LIMB_64
  // Toom variants aren't ported to Base2_64 yet (not in production dispatch).
  auto toom    = ToomCookMultiplication::Multiply(a, b, BigInteger::Base());
  auto toom5   = Toom5Multiplication::Multiply(a, b, BigInteger::Base());
  ASSERT_EQ(Compare(classic, toom), 0);
  ASSERT_EQ(Compare(classic, toom5), 0);
#endif
}

REGISTER_TEST(MulCross, Tiny_2x2)         { CrossMul(2,    2,    0x0001); }
REGISTER_TEST(MulCross, Small_16)         { CrossMul(16,   16,   0x0002); }
REGISTER_TEST(MulCross, Mid_64)           { CrossMul(64,   64,   0x0003); }
REGISTER_TEST(MulCross, Mid_256)          { CrossMul(256,  256,  0x0004); }
REGISTER_TEST(MulCross, Toom3Region_512)  { CrossMul(512,  512,  0x0005); }
REGISTER_TEST(MulCross, NTTRegion_1024)   { CrossMul(1024, 1024, 0x0006); }
REGISTER_TEST(MulCross, Skewed_512x16)    { CrossMul(512,  16,   0x0007); }
REGISTER_TEST(MulCross, Skewed_1024x32)   { CrossMul(1024, 32,   0x0008); }
REGISTER_TEST(MulCross, OneByMany)        { CrossMul(1,    300,  0x0009); }
REGISTER_TEST(MulCross, ManyByOne)        { CrossMul(300,  1,    0x000A); }

REGISTER_TEST(MulCross, MaxCarry257)
{
  // Every limb = 0xFFFFFFFF: longest possible carry chain in classic.
  std::vector<DataT> a(257, 0xFFFFFFFFu);
  std::vector<DataT> b(257, 0xFFFFFFFFu);
  auto classic = ClassicMultiplication::Multiply(a, b, BigInteger::Base());
  auto kara    = KaratsubaMultiplication::Multiply(a, b, BigInteger::Base());
  auto ntt     = NTTMultiplication::Multiply(a, b, BigInteger::Base());
  ASSERT_EQ(Compare(classic, kara), 0);
  ASSERT_EQ(Compare(classic, ntt),  0);
#if !BIGMATH_LIMB_64
  auto toom    = ToomCookMultiplication::Multiply(a, b, BigInteger::Base());
  auto toom5   = Toom5Multiplication::Multiply(a, b, BigInteger::Base());
  ASSERT_EQ(Compare(classic, toom), 0);
  ASSERT_EQ(Compare(classic, toom5), 0);
#endif
}

// ─── squaring: 4 algorithms agree, and equal Multiply(a, a) ──────────────────

static void CrossSquare(SizeT aLimbs, uint64_t seed)
{
  std::mt19937_64 gen(seed);
  auto a = RandomLimbs(aLimbs, gen);
  auto mul   = ClassicMultiplication::Multiply(a, a, BigInteger::Base());
  auto dispatched = Square(a, BigInteger::Base());
  ASSERT_EQ(Compare(mul, dispatched), 0);
#if !BIGMATH_LIMB_64
  // ClassicSquare/KaratsubaSquare/NTTSquare have no Base2_64 paths yet.
  auto csq   = ClassicSquare::Square(a, BigInteger::Base());
  auto ksq   = KaratsubaSquare::Square(a, BigInteger::Base());
  auto nsq   = NTTSquare::Square(a, BigInteger::Base());
  ASSERT_EQ(Compare(mul, csq), 0);
  ASSERT_EQ(Compare(mul, ksq), 0);
  ASSERT_EQ(Compare(mul, nsq), 0);
#endif
}

REGISTER_TEST(SquareCross, Small_4)    { CrossSquare(4,    0x100); }
REGISTER_TEST(SquareCross, Mid_64)     { CrossSquare(64,   0x101); }
REGISTER_TEST(SquareCross, Mid_256)    { CrossSquare(256,  0x102); }
REGISTER_TEST(SquareCross, Large_1024) { CrossSquare(1024, 0x103); }

// ─── division: 5 algorithms agree ────────────────────────────────────────────

static void CrossDiv(SizeT aLimbs, SizeT bLimbs, uint64_t seed)
{
  std::mt19937_64 gen(seed);
  auto a = RandomLimbs(aLimbs, gen);
  auto b = RandomLimbs(bLimbs, gen);

  auto dispatched = DivideAndRemainder(a, b, BigInteger::Base());
  auto bz         = BurnikelZieglerDivision::DivideAndRemainder(a, b, BigInteger::Base());
  auto fd         = FastDivision::DivideAndRemainder(a, b, BigInteger::Base());

  ASSERT_EQ(Compare(dispatched.first,  bz.first),  0);
  ASSERT_EQ(Compare(dispatched.second, bz.second), 0);
  ASSERT_EQ(Compare(dispatched.first,  fd.first),  0);
  ASSERT_EQ(Compare(dispatched.second, fd.second), 0);

  // Newton requires multi-limb b
  if (bLimbs > 1)
  {
    auto nt = NewtonDivision::DivideAndRemainder(a, b, BigInteger::Base());
    ASSERT_EQ(Compare(dispatched.first,  nt.first),  0);
    ASSERT_EQ(Compare(dispatched.second, nt.second), 0);

    auto rd = ReciprocalDivision::DivideAndRemainder(a, b, BigInteger::Base());
    ASSERT_EQ(Compare(dispatched.first,  rd.first),  0);
    ASSERT_EQ(Compare(dispatched.second, rd.second), 0);
  }

  // Identity: q*b + r == a, and r < b.
  auto recomposed = Add(Multiply(dispatched.first, b, BigInteger::Base()),
                        dispatched.second, BigInteger::Base());
  ASSERT_EQ(Compare(recomposed, a), 0);
  ASSERT_LT(Compare(dispatched.second, b), 0);
}

REGISTER_TEST(DivCross, Tiny_8x4)            { CrossDiv(8,    4,    0x200); }
REGISTER_TEST(DivCross, Mid_64x32)           { CrossDiv(64,   32,   0x201); }
REGISTER_TEST(DivCross, Balanced_256)        { CrossDiv(256,  256,  0x202); }
REGISTER_TEST(DivCross, AlmostBalanced_257)  { CrossDiv(257,  256,  0x203); }
REGISTER_TEST(DivCross, Balanced_1024)       { CrossDiv(1024, 1024, 0x204); }
REGISTER_TEST(DivCross, Skewed_4096_x_1100)  { CrossDiv(4096, 1100, 0x205); }
REGISTER_TEST(DivCross, Skewed_8192_x_512)   { CrossDiv(8192, 512,  0x206); }

REGISTER_TEST(DivCross, ScalarDivisor)
{
  std::mt19937_64 gen(0x300);
  auto a = RandomLimbs(512, gen);
  std::vector<DataT> b{0x12345678u};

  auto dispatched = DivideAndRemainder(a, b, BigInteger::Base());
  auto fd         = FastDivision::DivideAndRemainder(a, b, BigInteger::Base());
  auto bz         = BurnikelZieglerDivision::DivideAndRemainder(a, b, BigInteger::Base());

  ASSERT_EQ(Compare(dispatched.first,  fd.first),  0);
  ASSERT_EQ(Compare(dispatched.second, fd.second), 0);
  ASSERT_EQ(Compare(dispatched.first,  bz.first),  0);
  ASSERT_EQ(Compare(dispatched.second, bz.second), 0);

  // Classic scalar form
  auto cq = ClassicDivision::Divide(a, b[0], BigInteger::Base());
  ASSERT_EQ(Compare(dispatched.first, cq), 0);
}

// ─── ReciprocalDivision: cached divisor reused across numerators ─────────────

REGISTER_TEST(DivCross, ReciprocalReuseAcrossNumerators)
{
  std::mt19937_64 gen(0x400);
  auto b = RandomLimbs(128, gen);
  ReciprocalDivision::Divider divider(b, BigInteger::Base());

  for (SizeT na : {SizeT{128}, SizeT{129}, SizeT{256}, SizeT{1024}})
  {
    auto a = RandomLimbs(na, gen);
    auto qr_div = divider.DivideAndRemainder(a);
    auto qr_ref = DivideAndRemainder(a, b, BigInteger::Base());
    ASSERT_EQ(Compare(qr_div.first,  qr_ref.first),  0);
    ASSERT_EQ(Compare(qr_div.second, qr_ref.second), 0);
  }
}
