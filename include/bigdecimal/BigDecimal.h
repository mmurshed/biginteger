/**
 * BigDecimal — arbitrary-precision signed decimal numbers.
 *
 * Java-style fixed-point model: value = unscaled · 10^(-scale).
 * Positive scale = fractional digits; negative scale = trailing zeros.
 *
 * Add / Subtract / Multiply are exact. Divide requires an explicit target
 * scale and RoundingMode because the quotient is generally non-terminating.
 *
 * Internally backed by BigInteger for the unscaled value, so every arithmetic
 * cost reduces to BigInteger arithmetic plus one Pow10 multiplication during
 * rescaling.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGDECIMAL
#define BIGDECIMAL

#include <iosfwd>
#include <string>

#include "biginteger/BigInteger.h"
#include "bigdecimal/RoundingMode.h"

namespace BigMath
{
  class BigDecimal
  {
  public:
    // Constructors
    BigDecimal();                                              // 0, scale = 0
    BigDecimal(BigInteger unscaled, int scale);
    BigDecimal(long n);
    explicit BigDecimal(char const *s);
    explicit BigDecimal(std::string const &s);

    // Factories
    static BigDecimal FromString(std::string const &s);
    static BigDecimal FromLong(long n);
    static BigDecimal Zero();
    static BigDecimal One();
    static BigDecimal Ten();

    // Accessors
    BigInteger const &UnscaledValue() const { return unscaled_; }
    int  Scale()      const { return scale_; }
    int  Signum()     const;                       // -1 / 0 / +1
    bool IsNegative() const { return unscaled_.IsNegative(); }
    bool IsZero()     const { return unscaled_.Zero(); }
    int  Precision()  const;                       // number of significant digits

    // Arithmetic
    BigDecimal Add(BigDecimal const &o) const;
    BigDecimal Subtract(BigDecimal const &o) const;
    BigDecimal Multiply(BigDecimal const &o) const;
    BigDecimal Divide(BigDecimal const &o, int scale, RoundingMode mode) const;
    BigDecimal Negate() const;
    BigDecimal Abs() const;

    // Scale manipulation
    BigDecimal SetScale(int newScale, RoundingMode mode) const;
    BigDecimal StripTrailingZeros() const;

    // Comparison
    int  CompareTo(BigDecimal const &o) const;
    bool Equals(BigDecimal const &o) const;          // strict: unscaled + scale
    bool ValueEquals(BigDecimal const &o) const { return CompareTo(o) == 0; }

    // String form: plain (no scientific notation), with sign and decimal point.
    std::string ToPlainString() const;
    std::string ToString() const { return ToPlainString(); }

  private:
    BigInteger unscaled_;
    int        scale_;

    // Multiplies `value` by 10^n. n must be >= 0.
    static BigInteger MulPow10(BigInteger const &value, int n);

    // Aligns scales of two BigDecimals by scaling up the one with the smaller
    // scale; returns the common scale and the rescaled unscaled values.
    static void Align(BigDecimal const &a, BigDecimal const &b,
                      BigInteger &outA, BigInteger &outB, int &outScale);
  };

  // Operator overloads
  BigDecimal operator+(BigDecimal const &a, BigDecimal const &b);
  BigDecimal operator-(BigDecimal const &a, BigDecimal const &b);
  BigDecimal operator-(BigDecimal const &a);
  BigDecimal operator*(BigDecimal const &a, BigDecimal const &b);

  inline bool operator==(BigDecimal const &a, BigDecimal const &b) { return a.CompareTo(b) == 0; }
  inline bool operator!=(BigDecimal const &a, BigDecimal const &b) { return a.CompareTo(b) != 0; }
  inline bool operator< (BigDecimal const &a, BigDecimal const &b) { return a.CompareTo(b) <  0; }
  inline bool operator<=(BigDecimal const &a, BigDecimal const &b) { return a.CompareTo(b) <= 0; }
  inline bool operator> (BigDecimal const &a, BigDecimal const &b) { return a.CompareTo(b) >  0; }
  inline bool operator>=(BigDecimal const &a, BigDecimal const &b) { return a.CompareTo(b) >= 0; }

  std::ostream &operator<<(std::ostream &os, BigDecimal const &x);
}

#endif
