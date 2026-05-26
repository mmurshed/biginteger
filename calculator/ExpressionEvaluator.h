/**
 * BigMath expression evaluator.
 *
 * Recursive-descent parser over an arithmetic grammar:
 *
 *   summands := factors  ( ( '+' | '-' ) factors )*
 *   factors  := power    ( ( '*' | '/' | '%' ) power )*
 *   power    := unary    ( '^' power )?          // right-associative
 *   unary    := ( '+' | '-' ) unary | atom
 *   atom     := number | identifier | '(' summands ')'
 *
 * Variables: lookup function passed in via constructor. Unknown identifiers
 * raise EvalError with the source position.
 *
 * Operators: +, -, *, /, %, ^   (^ is integer exponent; non-negative).
 * Whitespace skipped between tokens. '#' to end of line is a comment.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef EXPRESSION_EVALUATOR
#define EXPRESSION_EVALUATOR

#include <cctype>
#include <functional>
#include <stdexcept>
#include <string>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/Operations.h"

namespace BigMath
{
  // Exception with the source-relative byte offset of the offending token.
  class EvalError : public std::runtime_error
  {
  public:
    EvalError(std::string const &msg, int position)
        : std::runtime_error(msg), pos(position) {}
    int Position() const { return pos; }

  private:
    int pos;
  };

  class ExpressionEvaluator
  {
  public:
    using Lookup = std::function<bool(std::string const &, BigInteger &)>;

    ExpressionEvaluator() = default;
    explicit ExpressionEvaluator(Lookup lookup) : lookup_(std::move(lookup)) {}

    BigInteger Eval(char const *expr)
    {
      origin_ = expr;
      cur_ = expr;
      BigInteger res = ParseSummands();
      SkipSpaces();
      if (*cur_ == ')')
        Fail("Unmatched closing parenthesis");
      if (*cur_ != '\0')
        Fail("Unexpected trailing input");
      return res;
    }

    BigInteger Eval(std::string const &expr) { return Eval(expr.c_str()); }

  private:
    char const *origin_ = nullptr;
    char const *cur_ = nullptr;
    Lookup lookup_;

    [[noreturn]] void Fail(std::string const &msg) const
    {
      throw EvalError(msg, (int)(cur_ - origin_));
    }

    void SkipSpaces()
    {
      while (*cur_ == ' ' || *cur_ == '\t')
        ++cur_;
      // '#' to end of line is a comment.
      if (*cur_ == '#')
      {
        while (*cur_ != '\0' && *cur_ != '\n')
          ++cur_;
      }
    }

    // ── grammar ───────────────────────────────────────────────────────────────

    BigInteger ParseSummands()
    {
      BigInteger a = ParseFactors();
      while (true)
      {
        SkipSpaces();
        char op = *cur_;
        if (op != '+' && op != '-')
          return a;
        ++cur_;
        BigInteger b = ParseFactors();
        a = (op == '+') ? (a + b) : (a - b);
      }
    }

    BigInteger ParseFactors()
    {
      BigInteger a = ParsePower();
      while (true)
      {
        SkipSpaces();
        char op = *cur_;
        if (op != '*' && op != '/' && op != '%')
          return a;
        ++cur_;
        BigInteger b = ParsePower();
        if (op == '*')
        {
          a = a * b;
        }
        else
        {
          if (b.Zero())
            Fail(op == '/' ? "Division by zero" : "Modulo by zero");
          if (op == '/')
            a = a / b;
          else
            a = a % b;
        }
      }
    }

    BigInteger ParsePower()
    {
      BigInteger base = ParseUnary();
      SkipSpaces();
      if (*cur_ != '^')
        return base;
      ++cur_;
      // Right-associative: 2^3^2 == 2^(3^2).
      BigInteger exp = ParsePower();
      return Pow(base, exp);
    }

    BigInteger ParseUnary()
    {
      SkipSpaces();
      if (*cur_ == '+')
      {
        ++cur_;
        return ParseUnary();
      }
      if (*cur_ == '-')
      {
        ++cur_;
        BigInteger x = ParseUnary();
        if (!x.Zero())
          -x;  // in-place negate
        return x;
      }
      return ParseAtom();
    }

    BigInteger ParseAtom()
    {
      SkipSpaces();

      if (*cur_ == '(')
      {
        ++cur_;
        BigInteger res = ParseSummands();
        SkipSpaces();
        if (*cur_ != ')')
          Fail("Expected ')'");
        ++cur_;
        return res;
      }

      if (*cur_ >= '0' && *cur_ <= '9')
      {
        // Hex (0x...) and binary (0b...) literals: parsed via repeated
        // multiply-and-add since the library only natively parses decimal.
        if (*cur_ == '0' && (cur_[1] == 'x' || cur_[1] == 'X'))
        {
          cur_ += 2;
          return ParseRadixLiteral(16);
        }
        if (*cur_ == '0' && (cur_[1] == 'b' || cur_[1] == 'B'))
        {
          cur_ += 2;
          return ParseRadixLiteral(2);
        }

        // Decimal: skip underscore separators if present.
        char const *start = cur_;
        char const *probe = cur_;
        bool hasUnderscore = false;
        while ((*probe >= '0' && *probe <= '9') || *probe == '_')
        {
          if (*probe == '_') hasUnderscore = true;
          ++probe;
        }
        if (hasUnderscore)
        {
          std::string clean;
          clean.reserve((size_t)(probe - start));
          for (char const *p = start; p < probe; ++p)
            if (*p != '_') clean.push_back(*p);
          cur_ = probe;
          return Parse(clean.c_str());
        }
        Int consumed = 0;
        BigInteger res = Parse(cur_, &consumed);
        cur_ += consumed;
        return res;
      }

      if (IsIdentifierStart(*cur_))
      {
        char const *start = cur_;
        while (IsIdentifierCont(*cur_))
          ++cur_;
        std::string name(start, cur_);
        BigInteger out;
        if (!lookup_ || !lookup_(name, out))
        {
          // Roll back to the start of the identifier so error position points
          // at the name, not past it.
          cur_ = start;
          Fail("Unknown identifier: " + name);
        }
        return out;
      }

      if (*cur_ == '\0')
        Fail("Unexpected end of input");
      Fail(std::string("Invalid character: '") + *cur_ + "'");
    }

    // ── radix literal parser (hex / binary) ───────────────────────────────────
    BigInteger ParseRadixLiteral(int radix)
    {
      char const *start = cur_;
      BigInteger result = BigIntegerBuilder::From("0");
      BigInteger base = BigIntegerBuilder::From(std::to_string(radix));
      int digits = 0;
      while (true)
      {
        char c = *cur_;
        int d;
        if (c == '_') { ++cur_; continue; }
        if (c >= '0' && c <= '9') d = c - '0';
        else if (c >= 'a' && c <= 'f') d = 10 + (c - 'a');
        else if (c >= 'A' && c <= 'F') d = 10 + (c - 'A');
        else break;
        if (d >= radix)
        {
          Fail(std::string("digit '") + c + "' out of range for base " + std::to_string(radix));
        }
        result = result * base + BigIntegerBuilder::From(std::to_string(d));
        ++cur_;
        ++digits;
      }
      if (digits == 0)
      {
        cur_ = start;
        Fail(std::string("expected ") + (radix == 16 ? "hex" : "binary") + " digits");
      }
      return result;
    }

    static bool IsIdentifierStart(char c)
    {
      return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_';
    }
    static bool IsIdentifierCont(char c)
    {
      return IsIdentifierStart(c) || (c >= '0' && c <= '9');
    }

    // ── integer power: square-and-multiply, exponent must be non-negative ────
    BigInteger Pow(BigInteger base, BigInteger exp)
    {
      if (exp.IsNegative())
        Fail("Negative exponent not supported for integer power");
      if (exp.Zero())
      {
        // 0^0 conventionally = 1.
        return BigIntegerBuilder::From("1");
      }
      if (base.Zero())
        return BigIntegerBuilder::From("0");

      // Reject absurd exponents that would never fit in memory (>= 2^32 bits
      // would already be terabytes-scale).
      if (exp.size() > 1 || exp[0] > (DataT)20000000)
        Fail("Exponent too large for integer power");

      uint32_t e = (uint32_t)exp[0];
      BigInteger result = BigIntegerBuilder::From("1");
      BigInteger b = base;
      while (e > 0)
      {
        if (e & 1u)
          result = result * b;
        e >>= 1;
        if (e > 0)
          b = b * b;
      }
      return result;
    }
  };
}

#endif
