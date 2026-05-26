// Evaluator: precedence, unary chaining, modulo, power, hex/bin literals,
// variable lookup, error positioning.

#include "unit_test_framework.h"

#include <map>
#include <string>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/Comparison.h"

#include "../../calculator/ExpressionEvaluator.h"

using namespace BigMath;

static std::string EvalDec(std::string const &src)
{
  ExpressionEvaluator e;
  return ToString(e.Eval(src));
}

REGISTER_TEST(Eval, Basic)
{
  ASSERT_EQ(EvalDec("1+2"),       "3");
  ASSERT_EQ(EvalDec("2*3+4"),     "10");
  ASSERT_EQ(EvalDec("2+3*4"),     "14");
  ASSERT_EQ(EvalDec("(2+3)*4"),   "20");
  ASSERT_EQ(EvalDec("10-3-2"),    "5");      // left-assoc
  ASSERT_EQ(EvalDec("100/10/2"),  "5");
  ASSERT_EQ(EvalDec("100 % 7"),   "2");
}

REGISTER_TEST(Eval, UnaryMinus)
{
  ASSERT_EQ(EvalDec("-5"),         "-5");
  ASSERT_EQ(EvalDec("-(-5)"),      "5");
  ASSERT_EQ(EvalDec("--5"),        "5");      // was broken pre-rewrite
  ASSERT_EQ(EvalDec("---5"),       "-5");
  ASSERT_EQ(EvalDec("3 - -4"),     "7");
  ASSERT_EQ(EvalDec("3 + -4"),     "-1");
}

REGISTER_TEST(Eval, Power)
{
  ASSERT_EQ(EvalDec("2^10"),       "1024");
  ASSERT_EQ(EvalDec("2^0"),        "1");
  ASSERT_EQ(EvalDec("0^0"),        "1");
  ASSERT_EQ(EvalDec("0^5"),        "0");
  ASSERT_EQ(EvalDec("2^3^2"),      "512");   // right-assoc: 2^(3^2)
  ASSERT_EQ(EvalDec("(-2)^3"),     "-8");
  ASSERT_EQ(EvalDec("3 * 2^4"),    "48");    // ^ tighter than *
}

REGISTER_TEST(Eval, Modulo)
{
  ASSERT_EQ(EvalDec("100 % 7"),       "2");
  ASSERT_EQ(EvalDec("2^100 % 1000"),  "376");
}

REGISTER_TEST(Eval, DigitSeparator)
{
  ASSERT_EQ(EvalDec("1_000_000"),     "1000000");
  ASSERT_EQ(EvalDec("1_000 * 1_000"), "1000000");
}

REGISTER_TEST(Eval, HexLiteral)
{
  ASSERT_EQ(EvalDec("0xFF"),       "255");
  ASSERT_EQ(EvalDec("0xff"),       "255");
  ASSERT_EQ(EvalDec("0xDEADBEEF"), "3735928559");
  ASSERT_EQ(EvalDec("0x1_0000_0000"), "4294967296");
}

REGISTER_TEST(Eval, BinaryLiteral)
{
  ASSERT_EQ(EvalDec("0b0"),        "0");
  ASSERT_EQ(EvalDec("0b1010"),     "10");
  ASSERT_EQ(EvalDec("0b1111_1111"),"255");
}

REGISTER_TEST(Eval, Whitespace)
{
  ASSERT_EQ(EvalDec("   1  +  2  "),  "3");
  ASSERT_EQ(EvalDec("\t2\t*\t3\t"),   "6");
}

REGISTER_TEST(Eval, CommentToEndOfLine)
{
  ASSERT_EQ(EvalDec("1 + 2 # ignored"), "3");
}

REGISTER_TEST(Eval, LargePower)
{
  std::string s = EvalDec("2^100");
  ASSERT_EQ(s, "1267650600228229401496703205376");
}

REGISTER_TEST(Eval, Variables)
{
  std::map<std::string, BigInteger> env;
  env["x"] = BigIntegerBuilder::From("10");
  env["y"] = BigIntegerBuilder::From("7");
  ExpressionEvaluator e(
      [&](std::string const &n, BigInteger &out) {
        auto it = env.find(n);
        if (it == env.end()) return false;
        out = it->second;
        return true;
      });
  ASSERT_EQ(ToString(e.Eval("x * y + 1")), "71");
  ASSERT_EQ(ToString(e.Eval("x^2 - y^2")), "51");
}

REGISTER_TEST(Eval, ErrorUnknownIdentifier)
{
  ExpressionEvaluator e;
  bool threw = false;
  try { (void)e.Eval("1 + foo"); }
  catch (EvalError const &ex) {
    threw = true;
    ASSERT_EQ(ex.Position(), 4);       // caret at 'foo'
  }
  ASSERT_TRUE(threw);
}

REGISTER_TEST(Eval, ErrorDivByZero)
{
  ExpressionEvaluator e;
  bool threw = false;
  try { (void)e.Eval("1 / 0"); }
  catch (EvalError const &) { threw = true; }
  ASSERT_TRUE(threw);
}

REGISTER_TEST(Eval, ErrorNegativeExponent)
{
  ExpressionEvaluator e;
  bool threw = false;
  try { (void)e.Eval("2 ^ -1"); }
  catch (EvalError const &) { threw = true; }
  ASSERT_TRUE(threw);
}

REGISTER_TEST(Eval, ErrorTrailingJunk)
{
  ExpressionEvaluator e;
  bool threw = false;
  try { (void)e.Eval("1 + 2 garbage"); }
  catch (EvalError const &) { threw = true; }
  ASSERT_TRUE(threw);
}

REGISTER_TEST(Eval, ErrorParenMismatch)
{
  ExpressionEvaluator e;
  for (char const *bad : {"(1+2", "1+2)", "((1+2)"})
  {
    bool threw = false;
    try { (void)e.Eval(bad); }
    catch (EvalError const &) { threw = true; }
    ASSERT_TRUE(threw);
  }
}
