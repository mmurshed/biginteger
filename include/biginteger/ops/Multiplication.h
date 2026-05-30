/**
 * BigMath: BigInteger × BigInteger multiplication operator.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_MULTIPLICATION
#define BIGINTEGER_MULTIPLICATION

#include "../BigInteger.h"
#include "../algorithms/Multiplication.h"
#include "../algorithms/multiplication/NTTMultiplication.h"

namespace BigMath
{
  class PreparedMultiplication
  {
  private:
    BigInteger operand;
    NTTMultiplication::PreparedOperand prepared;

  public:
    PreparedMultiplication(BigInteger const &a, SizeT maxOtherLimbs)
        : operand(a),
          prepared(NTTMultiplication::PrepareOperand(a.GetInteger(), maxOtherLimbs, BigInteger::Base()))
    {
    }

    BigInteger Multiply(BigInteger const &b) const
    {
      auto result = NTTMultiplication::Multiply(prepared, b.GetInteger());
      return BigInteger(result, operand.IsNegative() != b.IsNegative());
    }

    BigInteger operator*(BigInteger const &b) const
    {
      return Multiply(b);
    }

    BigInteger const &Operand() const
    {
      return operand;
    }

    NTTMultiplication::PreparedOperand const &Prepared() const
    {
      return prepared;
    }
  };

  inline PreparedMultiplication PrepareMultiplication(BigInteger const &a, SizeT maxOtherLimbs)
  {
    return PreparedMultiplication(a, maxOtherLimbs);
  }

  BigInteger Multiply(BigInteger const &a, BigInteger const &b);
  BigInteger operator*(BigInteger const &a, BigInteger const &b);
}

#endif
