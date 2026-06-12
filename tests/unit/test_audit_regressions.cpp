// Regression tests for the 2026-06-12 codebase audit (PR #120).
// Each test pins a fixed bug; comments name the original failure mode.

#include "unit_test_framework.h"

#include <sstream>
#include <type_traits>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/ops/Operations.h"
#include "biginteger/ops/IO.h"
#include "biginteger/algorithms/BaseConversion.h"
#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Subtraction.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/division/FastDivision.h"
#include "biginteger/algorithms/division/KnuthDivision.h"
#include "biginteger/algorithms/division/QuotientSizedDivision.h"

using namespace BigMath;

namespace
{
  BigInteger BI(const char *s) { return BigIntegerBuilder::From(s); }

  constexpr DataT kMax = ~(DataT)0;

  // q*b + r == a and |r| < |b| for vector-level results.
  bool IdentityHolds(std::vector<DataT> const &a, std::vector<DataT> const &b,
                     std::pair<std::vector<DataT>, std::vector<DataT>> const &qr)
  {
    auto qb = Multiply(qr.first, b, BigInteger::Base());
    auto back = Add(qb, qr.second, BigInteger::Base());
    return Compare(back, a) == 0 && Compare(qr.second, b) < 0;
  }
}

// ── remainder sign: truncated division, r follows the dividend ──────────────
// Was: r sign computed as (a.IsNegative() || b.IsNegative()); 7 % -3 == -1
// and q*b + r != a for the positive-dividend/negative-divisor quadrant.

REGISTER_TEST(AuditDivision, RemainderSignAllQuadrants)
{
  const char *as[] = {"7", "-7", "7", "-7"};
  const char *bs[] = {"3", "3", "-3", "-3"};
  const char *eq[] = {"2", "-2", "-2", "2"};
  const char *er[] = {"1", "-1", "1", "-1"};
  for (int i = 0; i < 4; i++)
  {
    BigInteger a = BI(as[i]), b = BI(bs[i]);
    auto [q, r] = DivideAndRemainder(a, b);
    ASSERT_EQ(ToString(q), std::string(eq[i]));
    ASSERT_EQ(ToString(r), std::string(er[i]));
    BigInteger back = q * b + r;
    ASSERT_EQ(back.CompareTo(a), 0);
  }
}

REGISTER_TEST(AuditDivision, CachedDivisionRemainderSign)
{
  // Same bug existed in the duplicated CachedDivision copy.
  auto [q, r] = CacheDivision(BI("-3")).DivideAndRemainder(BI("7"));
  ASSERT_EQ(ToString(q), std::string("-2"));
  ASSERT_EQ(ToString(r), std::string("1"));
}

// ── unary minus: value-returning, operand unchanged ──────────────────────────
// Was: operator-() flipped the sign in place and returned an lvalue ref, so
// `BigInteger y = -x;` corrupted x.

REGISTER_TEST(AuditBigInteger, UnaryMinusDoesNotMutate)
{
  BigInteger x = BI("5");
  BigInteger y = -x;
  ASSERT_EQ(ToString(x), std::string("5"));
  ASSERT_EQ(ToString(y), std::string("-5"));
  // -(-y) round-trips.
  ASSERT_EQ((-(-y)).CompareTo(y), 0);
}

// ── explicit size ctor: no implicit int → BigInteger zero ───────────────────
// Was: BigInteger(SizeT size = 0, ...) was implicit, so `a + 5` resolved to
// a + BigInteger(5) == a + 0 and `a == 100` compared against zero.

REGISTER_TEST(AuditBigInteger, SizeCtorIsExplicit)
{
  static_assert(!std::is_convertible<SizeT, BigInteger>::value,
                "BigInteger(SizeT) must be explicit");
  static_assert(!std::is_convertible<int, BigInteger>::value,
                "int must not implicitly convert to BigInteger");
  ASSERT_TRUE(true);
}

// ── FastDivision: MG 3/2 top-equal special case ──────────────────────────────
// Was: when a remainder window's top two limbs equal the divisor's (d1:d0),
// MGQhat's modular arithmetic wrapped and produced a garbage digit (the true
// digit is exactly B-1). Shapes below reproduce the original mismatches.

REGISTER_TEST(AuditDivision, FastDivisionMGTopEqual)
{
  // Original fuzz counterexample (24/4 limbs, d1 = 2^64-1, d0 = 0).
  std::vector<DataT> a = {
      kMax, 8632100396935015419ull, kMax, 4438608630380892118ull, kMax,
      84453682283771477ull, kMax, 10368953306968301788ull,
      11214772461763125535ull, kMax, 5816448965103977582ull, 0,
      18261805052452885511ull, kMax, kMax, 3462971057439967324ull,
      1809960185198586381ull, 0, 0, kMax, 7343587919163703040ull,
      7165561623291364460ull, 0, kMax};
  std::vector<DataT> b = {kMax, 7917626342253556694ull, 0, kMax};
  auto qr = FastDivision::DivideAndRemainder(a, b, BigInteger::Base());
  ASSERT_TRUE(IdentityHolds(a, b, qr));
}

REGISTER_TEST(AuditDivision, FastDivisionMGTopEqualMinimal)
{
  // Minimal trigger: dividend window top two limbs == divisor top two limbs,
  // divisor has a nonzero limb below them.
  std::vector<DataT> b = {1, 5, kMax};           // (d1:d0) = (kMax, 5)
  std::vector<DataT> a = {0, 0, 0, 5, kMax};     // window tops hit (kMax, 5)
  auto qr = FastDivision::DivideAndRemainder(a, b, BigInteger::Base());
  ASSERT_TRUE(IdentityHolds(a, b, qr));
}

REGISTER_TEST(AuditDivision, FastDivisionStructuredLimbSweep)
{
  // Deterministic sweep over structured limb patterns (0 / 1 / max) that the
  // random correctness harness never produces.
  const DataT vals[] = {0, 1, kMax, (DataT)1 << 63, kMax - 1};
  int idx = 0;
  for (int na = 3; na <= 7; ++na)
    for (int nb = 2; nb <= 4 && nb < na; ++nb)
    {
      std::vector<DataT> a(na), b(nb);
      for (int i = 0; i < na; ++i) a[i] = vals[(idx + i) % 5];
      for (int i = 0; i < nb; ++i) b[i] = vals[(idx + 2 * i + 1) % 5];
      ++idx;
      if (b.back() == 0) b.back() = kMax;
      while (a.size() > 1 && a.back() == 0) a.pop_back();
      auto qr = FastDivision::DivideAndRemainder(a, b, BigInteger::Base());
      ASSERT_TRUE(IdentityHolds(a, b, qr));
    }
}

// ── KnuthDivision: works under the default base, agrees with FastDivision ───
// Was: d = base / (b.back()+1) divided by zero under the Base2_64 sentinel;
// qhat*v[i] overflowed signed Long under Base2_32; untrimmed operands
// overflowed the heap because Multiply trims trailing zeros.

REGISTER_TEST(AuditDivision, KnuthMatchesFastDefaultBase)
{
  const std::vector<std::pair<std::vector<DataT>, std::vector<DataT>>> cases = {
      {{1, 2, 3, 4, 5}, {7, 9}},
      {{kMax, kMax, kMax, kMax}, {kMax, kMax}},
      {{0, 0, kMax, kMax}, {1, kMax}},
      {{5, kMax, kMax, kMax, kMax}, {kMax, 1, kMax}},
      {{0, 0, 0, kMax - 1, kMax}, {1, 0, kMax}}, // qhat-clamp shape
  };
  for (auto const &[a, b] : cases)
  {
    auto k = KnuthDivision::DivideAndRemainder(a, b, BigInteger::Base());
    auto f = FastDivision::DivideAndRemainder(a, b, BigInteger::Base());
    ASSERT_EQ(Compare(k.first, f.first), 0);
    ASSERT_EQ(Compare(k.second, f.second), 0);
  }
}

REGISTER_TEST(AuditDivision, KnuthUntrimmedOperands)
{
  // High zero limbs on both operands: previously a heap-buffer-overflow in
  // KnuthDivision (Multiply trims trailing zeros; the loop indexed past).
  // FastDivision's contract requires trimmed inputs, so compare against a
  // trimmed-input reference call.
  std::vector<DataT> a = {7, 9, 11, 0, 0};
  std::vector<DataT> b = {kMax, 3, 0};
  auto k = KnuthDivision::DivideAndRemainder(a, b, BigInteger::Base());
  std::vector<DataT> at = {7, 9, 11};
  std::vector<DataT> bt = {kMax, 3};
  auto f = FastDivision::DivideAndRemainder(at, bt, BigInteger::Base());
  ASSERT_EQ(Compare(k.first, f.first), 0);
  ASSERT_EQ(Compare(k.second, f.second), 0);
}

// ── QuotientSizedDivision: contract violation falls back, no underflow ──────
// Was: delta + GUARD >= nb made `s = nb - t` wrap as unsigned and index UB.

REGISTER_TEST(AuditDivision, QuotientSizedContractFallback)
{
  // na - nb + GUARD >= nb violates the dispatch contract on purpose.
  std::vector<DataT> a(40, 3);
  a.back() = kMax;
  std::vector<DataT> b(8, 7);
  b.back() = kMax;
  auto qr = QuotientSizedDivision::DivideAndRemainder(a, b, BigInteger::Base());
  ASSERT_TRUE(IdentityHolds(a, b, qr));
}

// ── AddTo(vector, vector): grows for the final carry ────────────────────────
// Was: a += b never grew `a`; {2^64-1} + {1} silently truncated to {0}.

REGISTER_TEST(AuditAddition, AddToGrowsOnFinalCarry)
{
  std::vector<DataT> a{kMax};
  std::vector<DataT> b{1};
  AddTo(a, b, BigInteger::Base());
  ASSERT_EQ(a.size(), (size_t)2);
  ASSERT_EQ(a[0], (DataT)0);
  ASSERT_EQ(a[1], (DataT)1);
}

REGISTER_TEST(AuditAddition, AddToCarryChainGrowth)
{
  std::vector<DataT> a{kMax, kMax, kMax};
  std::vector<DataT> b{1};
  AddTo(a, b, BigInteger::Base());
  std::vector<DataT> want{0, 0, 0, 1};
  ASSERT_EQ(Compare(a, want), 0);
}

// ── windowed Add: final carry propagates past one slot ──────────────────────
// Was: `result[rPos] += carry` could overflow the slot without propagating.

REGISTER_TEST(AuditAddition, WindowedAddCarryPropagates)
{
  // result has headroom; the carry out of the window lands on a kMax limb
  // and must ripple one further.
  std::vector<DataT> a{kMax};
  std::vector<DataT> b{1};
  std::vector<DataT> result{0, kMax, 0};
  Add(a, 0, 0, b, 0, 0, result, 0, BigInteger::Base());
  std::vector<DataT> want{0, 0, 1};
  ASSERT_EQ(Compare(result, want), 0);
}

// ── scalar SubtractFrom: multi-limb amounts in small-base mode ───────────────
// Was: assumed a single borrow; b >= base corrupted the result in Base2_32.

REGISTER_TEST(AuditSubtraction, ScalarSubtractMultiLimbBase32)
{
  // a = {5, 4} base-2^32 = 4*2^32 + 5; subtract 2^33 + 3 (spans two limbs).
  // Result: 2*2^32 + 2 = {2, 2}.
  std::vector<DataT> c{5, 4};
  SubtractFrom(c, 0, 1, ((ULong)1 << 33) + 3, Base2_32);
  std::vector<DataT> want{2, 2};
  ASSERT_EQ(Compare(c, want), 0);
}

REGISTER_TEST(AuditSubtraction, ScalarSubtractBorrowChainBase64)
{
  std::vector<DataT> a{0, 0, 1}; // 2^128
  SubtractFrom(a, 0, 2, 1, Base2_64);
  std::vector<DataT> want{kMax, kMax};
  ASSERT_EQ(Compare(a, want), 0);
}

// ── operator>>: failbit on non-numeric input, "0" stays valid ───────────────
// Was: "abc" / lone "-" silently produced 0 with the stream still good(),
// so `while (cin >> x)` loops spun forever on bad input.

REGISTER_TEST(AuditIO, ExtractionSetsFailbit)
{
  std::istringstream bad("abc");
  BigInteger v;
  bad >> v;
  ASSERT_TRUE(bad.fail());

  std::istringstream dash("-");
  BigInteger v2;
  dash >> v2;
  ASSERT_TRUE(dash.fail());

  std::istringstream zero("0");
  BigInteger v3;
  zero >> v3;
  ASSERT_FALSE(zero.fail());
  ASSERT_TRUE(v3.Zero());

  std::istringstream zeros("0000");
  BigInteger v4;
  zeros >> v4;
  ASSERT_FALSE(zeros.fail());
  ASSERT_TRUE(v4.Zero());

  std::istringstream ok("-42xyz");
  BigInteger v5;
  ok >> v5;
  ASSERT_FALSE(ok.fail());
  ASSERT_EQ(ToString(v5), std::string("-42"));
}

// ── ConvertBase: subrange honored when base1 == base2 ────────────────────────
// Was: returned a copy of the WHOLE vector regardless of [start, end].

REGISTER_TEST(AuditConversion, ConvertBaseSameBaseRange)
{
  std::vector<DataT> v{10, 20, 30, 40, 50};
  auto out = ConvertBase(v, 1, 3, Base2_32, Base2_32);
  std::vector<DataT> want{20, 30, 40};
  ASSERT_EQ(out.size(), want.size());
  for (size_t i = 0; i < want.size(); ++i)
    ASSERT_EQ(out[i], want[i]);
}

// ── Compare(vector, scalar): untrimmed vectors compare by value ──────────────
// Was: {0, 0} (untrimmed zero) compared greater than any scalar because the
// raw size() > 1 check ran before trimming.

REGISTER_TEST(AuditComparison, ScalarCompareIgnoresLeadingZeros)
{
  std::vector<DataT> zeroes{0, 0};
  ASSERT_EQ(Compare(zeroes, (DataT)0), 0);
  ASSERT_EQ(Compare(zeroes, (DataT)5), -1);

  std::vector<DataT> five{5, 0, 0};
  ASSERT_EQ(Compare(five, (DataT)5), 0);
  ASSERT_EQ(Compare(five, (DataT)4), 1);
  ASSERT_EQ(Compare(five, (DataT)6), -1);
}
