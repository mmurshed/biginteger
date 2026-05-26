/**
 * BigDecimal implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "bigdecimal/BigDecimal.h"

#include <algorithm>
#include <cctype>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/Operations.h"

namespace BigMath
{
  namespace
  {
    // Negation that does not mutate: BigInteger's operator- is in-place.
    inline BigInteger Negated(BigInteger const &x)
    {
      BigInteger r(x);
      if (!r.Zero()) -r;
      return r;
    }

    inline BigInteger AbsCopy(BigInteger const &x)
    {
      BigInteger r(x);
      r.SetSign(false);
      return r;
    }

    inline BigInteger BIOne()  { return BigIntegerBuilder::From("1"); }
    inline BigInteger BIZero() { return BigIntegerBuilder::From("0"); }
    inline BigInteger BITen()  { return BigIntegerBuilder::From("10"); }

    // 10^n as BigInteger. n must be >= 0.
    //
    // The underlying Pow10(d) returns a thread-local cached vector<DataT>; we
    // additionally cache the BigInteger wrapper here so repeated Multiply/Add
    // calls at the same scale don't pay for the vector copy + sign assignment
    // every time. Wrappers reference the same cached value vector.
    BigInteger const &Pow10Bi(int n)
    {
      if (n < 0)
        throw std::invalid_argument("Pow10Bi: negative exponent");
      static thread_local std::unordered_map<int, BigInteger> cache;
      auto it = cache.find(n);
      if (it != cache.end()) return it->second;
      BigInteger v = (n == 0) ? BIOne() : BigInteger(Pow10((SizeT)n), false);
      return cache.emplace(n, std::move(v)).first->second;
    }

    // True ⇔ rounded result magnitude should increase by 1 (away from zero).
    // `r` and `den` are non-negative magnitudes; `q` is the truncated quotient
    // magnitude. `resultNeg` is the sign of the *final* signed quotient.
    bool RoundAwayFromZero(BigInteger const &r, BigInteger const &den,
                           BigInteger const &q, RoundingMode mode, bool resultNeg)
    {
      if (r.Zero()) return false;
      switch (mode)
      {
        case RoundingMode::DOWN:    return false;
        case RoundingMode::UP:      return true;
        case RoundingMode::CEILING: return !resultNeg;
        case RoundingMode::FLOOR:   return  resultNeg;
        case RoundingMode::HALF_UP:
        case RoundingMode::HALF_DOWN:
        case RoundingMode::HALF_EVEN:
        {
          BigInteger twoR = r + r;
          Int c = twoR.CompareTo(den);
          if (c >  0) return true;
          if (c <  0) return false;
          if (mode == RoundingMode::HALF_UP)   return true;
          if (mode == RoundingMode::HALF_DOWN) return false;
          // HALF_EVEN: tie → round to even ⇒ round away iff q is odd
          if (q.size() == 0) return false;
          return (q[0] & 1u) != 0;
        }
        case RoundingMode::UNNECESSARY:
          throw std::runtime_error("BigDecimal: rounding necessary");
      }
      return false;
    }
  }

  // ── construction ──────────────────────────────────────────────────────────

  BigDecimal::BigDecimal() : unscaled_(BIZero()), scale_(0) {}

  BigDecimal::BigDecimal(BigInteger unscaled, int scale)
      : unscaled_(std::move(unscaled)), scale_(scale)
  {
    if (unscaled_.Zero())
      scale_ = scale;  // keep caller-requested scale even for zero (Java does)
  }

  BigDecimal::BigDecimal(long n) : unscaled_(BigIntegerBuilder::From((Long)n)), scale_(0) {}

  BigDecimal::BigDecimal(char const *s)        : BigDecimal(FromString(s)) {}
  BigDecimal::BigDecimal(std::string const &s) : BigDecimal(FromString(s)) {}

  BigDecimal BigDecimal::FromLong(long n) { return BigDecimal(n); }
  BigDecimal BigDecimal::Zero() { return BigDecimal(); }
  BigDecimal BigDecimal::One()  { return BigDecimal(BIOne(),  0); }
  BigDecimal BigDecimal::Ten()  { return BigDecimal(BITen(),  0); }

  // Grammar (java.math.BigDecimal style):
  //   [+-]? ( digits ('.' digits?)? | '.' digits ) ([eE] [+-]? digits)?
  BigDecimal BigDecimal::FromString(std::string const &s)
  {
    size_t i = 0, n = s.size();
    if (n == 0)
      throw std::invalid_argument("BigDecimal: empty string");

    bool neg = false;
    if (s[i] == '+' || s[i] == '-')
    {
      neg = (s[i] == '-');
      ++i;
    }

    std::string digits;        // concatenated integer + fractional, no '.'
    digits.reserve(n);
    int fracLen = 0;
    bool seenDigit = false;
    bool seenDot = false;

    while (i < n)
    {
      char c = s[i];
      if (c >= '0' && c <= '9') { digits.push_back(c); if (seenDot) ++fracLen; seenDigit = true; ++i; }
      else if (c == '.' && !seenDot) { seenDot = true; ++i; }
      else break;
    }
    if (!seenDigit)
      throw std::invalid_argument("BigDecimal: no digits");

    long expPart = 0;
    if (i < n && (s[i] == 'e' || s[i] == 'E'))
    {
      ++i;
      bool eNeg = false;
      if (i < n && (s[i] == '+' || s[i] == '-'))
      {
        eNeg = (s[i] == '-');
        ++i;
      }
      if (i >= n || !(s[i] >= '0' && s[i] <= '9'))
        throw std::invalid_argument("BigDecimal: exponent missing digits");
      while (i < n && s[i] >= '0' && s[i] <= '9')
      {
        expPart = expPart * 10 + (s[i] - '0');
        if (expPart > 1'000'000'000L)
          throw std::out_of_range("BigDecimal: exponent overflow");
        ++i;
      }
      if (eNeg) expPart = -expPart;
    }
    if (i != n)
      throw std::invalid_argument("BigDecimal: trailing junk in '" + s + "'");

    // Strip leading zeros from digits, but keep at least one.
    size_t firstNon = digits.find_first_not_of('0');
    if (firstNon == std::string::npos)
      digits = "0";
    else if (firstNon > 0)
      digits.erase(0, firstNon);

    BigInteger unscaled = BigIntegerBuilder::From(digits);
    if (neg && !unscaled.Zero())
      -unscaled;

    long finalScale = (long)fracLen - expPart;
    if (finalScale > 2'000'000'000L || finalScale < -2'000'000'000L)
      throw std::out_of_range("BigDecimal: scale out of int range");

    return BigDecimal(std::move(unscaled), (int)finalScale);
  }

  // ── accessors ─────────────────────────────────────────────────────────────

  int BigDecimal::Signum() const
  {
    if (unscaled_.Zero()) return 0;
    return unscaled_.IsNegative() ? -1 : 1;
  }

  int BigDecimal::Precision() const
  {
    if (unscaled_.Zero()) return 1;
    std::string s = BigMath::ToString(unscaled_);
    int n = (int)s.size();
    if (!s.empty() && s[0] == '-') --n;
    return n;
  }

  // ── helpers ───────────────────────────────────────────────────────────────

  BigInteger BigDecimal::MulPow10(BigInteger const &value, int n)
  {
    if (n == 0) return value;
    if (value.Zero()) return BIZero();
    return value * Pow10Bi(n);
  }

  void BigDecimal::Align(BigDecimal const &a, BigDecimal const &b,
                         BigInteger &outA, BigInteger &outB, int &outScale)
  {
    if (a.scale_ == b.scale_)
    {
      outA = a.unscaled_;
      outB = b.unscaled_;
      outScale = a.scale_;
      return;
    }
    if (a.scale_ < b.scale_)
    {
      outA = MulPow10(a.unscaled_, b.scale_ - a.scale_);
      outB = b.unscaled_;
      outScale = b.scale_;
    }
    else
    {
      outA = a.unscaled_;
      outB = MulPow10(b.unscaled_, a.scale_ - b.scale_);
      outScale = a.scale_;
    }
  }

  // ── arithmetic ────────────────────────────────────────────────────────────

  BigDecimal BigDecimal::Add(BigDecimal const &o) const
  {
    BigInteger a, b; int s;
    Align(*this, o, a, b, s);
    return BigDecimal(a + b, s);
  }

  BigDecimal BigDecimal::Subtract(BigDecimal const &o) const
  {
    BigInteger a, b; int s;
    Align(*this, o, a, b, s);
    return BigDecimal(a - b, s);
  }

  BigDecimal BigDecimal::Multiply(BigDecimal const &o) const
  {
    long ns = (long)scale_ + (long)o.scale_;
    if (ns > 2'000'000'000L || ns < -2'000'000'000L)
      throw std::out_of_range("BigDecimal: multiply scale overflow");
    return BigDecimal(unscaled_ * o.unscaled_, (int)ns);
  }

  BigDecimal BigDecimal::Negate() const
  {
    return BigDecimal(Negated(unscaled_), scale_);
  }

  BigDecimal BigDecimal::Abs() const
  {
    if (!unscaled_.IsNegative()) return *this;
    return BigDecimal(AbsCopy(unscaled_), scale_);
  }

  BigDecimal BigDecimal::Divide(BigDecimal const &o, int newScale, RoundingMode mode) const
  {
    if (o.IsZero())
      throw std::invalid_argument("BigDecimal: division by zero");

    bool resultNeg = unscaled_.IsNegative() != o.unscaled_.IsNegative();

    BigInteger absA = AbsCopy(unscaled_);
    BigInteger absO = AbsCopy(o.unscaled_);

    // result.unscaled ≈ a/o · 10^newScale
    //   = (a.unscaled · 10^(newScale + o.scale - a.scale)) / o.unscaled
    long shiftL = (long)newScale + (long)o.scale_ - (long)scale_;
    if (shiftL > 2'000'000'000L || shiftL < -2'000'000'000L)
      throw std::out_of_range("BigDecimal: divide shift overflow");
    int shift = (int)shiftL;

    BigInteger num, den;
    if (shift >= 0)
    {
      num = MulPow10(absA, shift);
      den = absO;
    }
    else
    {
      num = absA;
      den = MulPow10(absO, -shift);
    }

    auto [q, r] = DivideAndRemainder(num, den);

    if (RoundAwayFromZero(r, den, q, mode, resultNeg))
      q = q + BIOne();

    BigInteger signed_q = q;
    if (resultNeg && !signed_q.Zero())
      -signed_q;

    return BigDecimal(std::move(signed_q), newScale);
  }

  // ── scale manipulation ────────────────────────────────────────────────────

  BigDecimal BigDecimal::SetScale(int newScale, RoundingMode mode) const
  {
    if (newScale == scale_) return *this;
    if (newScale > scale_)
    {
      // pure rescale up: multiply unscaled by 10^(newScale - scale)
      return BigDecimal(MulPow10(unscaled_, newScale - scale_), newScale);
    }
    // newScale < scale: divide unscaled by 10^(scale - newScale) with rounding
    int drop = scale_ - newScale;
    BigInteger divisor = Pow10Bi(drop);
    bool resultNeg = unscaled_.IsNegative();
    BigInteger absU = AbsCopy(unscaled_);
    auto [q, r] = DivideAndRemainder(absU, divisor);
    if (RoundAwayFromZero(r, divisor, q, mode, resultNeg))
      q = q + BIOne();
    BigInteger signed_q = q;
    if (resultNeg && !signed_q.Zero())
      -signed_q;
    return BigDecimal(std::move(signed_q), newScale);
  }

  BigDecimal BigDecimal::StripTrailingZeros() const
  {
    if (unscaled_.Zero()) return BigDecimal(BIZero(), 0);
    BigInteger u = unscaled_;
    int s = scale_;
    BigInteger ten = BITen();
    while (s > 0)
    {
      auto [q, r] = DivideAndRemainder(u, ten);
      if (!r.Zero()) break;
      u = q;
      --s;
    }
    return BigDecimal(std::move(u), s);
  }

  // ── comparison ────────────────────────────────────────────────────────────

  int BigDecimal::CompareTo(BigDecimal const &o) const
  {
    int sa = Signum(), sb = o.Signum();
    if (sa != sb) return sa < sb ? -1 : 1;
    if (sa == 0) return 0;
    BigInteger a, b; int s;
    Align(*this, o, a, b, s);
    Int c = a.CompareTo(b);
    return c < 0 ? -1 : (c > 0 ? 1 : 0);
  }

  bool BigDecimal::Equals(BigDecimal const &o) const
  {
    return scale_ == o.scale_ && unscaled_.CompareTo(o.unscaled_) == 0;
  }

  // ── format ────────────────────────────────────────────────────────────────

  std::string BigDecimal::ToPlainString() const
  {
    std::string s = BigMath::ToString(unscaled_);
    bool neg = (!s.empty() && s[0] == '-');
    std::string digits = neg ? s.substr(1) : s;

    if (scale_ == 0)
      return s;

    if (scale_ < 0)
    {
      // value = unscaled · 10^(-scale); append (-scale) zeros
      std::string out;
      out.reserve(s.size() + (size_t)(-scale_));
      out += s;
      out.append((size_t)(-scale_), '0');
      return out;
    }

    // scale_ > 0: insert decimal point or pad with leading zeros
    int dlen = (int)digits.size();
    std::string out;
    out.reserve(digits.size() + (size_t)scale_ + 4);
    if (neg) out.push_back('-');
    if (dlen > scale_)
    {
      out.append(digits, 0, (size_t)(dlen - scale_));
      out.push_back('.');
      out.append(digits, (size_t)(dlen - scale_), (size_t)scale_);
    }
    else
    {
      out += "0.";
      out.append((size_t)(scale_ - dlen), '0');
      out += digits;
    }
    return out;
  }

  // ── operators ─────────────────────────────────────────────────────────────

  BigDecimal operator+(BigDecimal const &a, BigDecimal const &b) { return a.Add(b); }
  BigDecimal operator-(BigDecimal const &a, BigDecimal const &b) { return a.Subtract(b); }
  BigDecimal operator-(BigDecimal const &a)                       { return a.Negate(); }
  BigDecimal operator*(BigDecimal const &a, BigDecimal const &b) { return a.Multiply(b); }

  std::ostream &operator<<(std::ostream &os, BigDecimal const &x)
  {
    return os << x.ToPlainString();
  }
}
