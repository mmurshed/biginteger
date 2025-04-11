/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_MULTIPLICATION
#define BIGINTEGER_MULTIPLICATION

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../algorithms/classic/ClassicMultiplication.h"
#include "../algorithms/karatsuba/KaratsubaMultiplication.h"
#include "../algorithms/toomcook/ToomCookMultiplication2.h"
#include "../algorithms/stonehagestrassen/FFTMultiplication.h"

namespace BigMath
{
  class BigIntegerMultiplication
  {
    // Karatsuba performs better for over 128 digits for the result
    static const SizeT KARATSUBA_THRESHOLD = 128;
    static const SizeT TOOM_COOK_THRESHOLD = 512;
    static const SizeT FFT_THRESHOLD = 1024;

  public:
    static vector<DataT> MultiplyUnsigned(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
    {
      SizeT size = a.size() + b.size();

      if (size <= KARATSUBA_THRESHOLD)
      {
        return ClassicMultiplication::Multiply(a, b, base);
      }
      else if (size <= TOOM_COOK_THRESHOLD)
      {
        return KaratsubaMultiplication::Multiply(a, b, base);
      }
      else if (size <= FFT_THRESHOLD)
      {
        ToomCookMultiplication2 toomCook;
        return toomCook.Multiply(a, b, base);
      }

      return FFTMultiplication::Multiply(a, b, base);
    }

    static BigInteger Multiply(BigInteger const &a, BigInteger const &b)
    {
      if (a.IsZero() || b.IsZero())
        return BigInteger(); // 0 times anything is zero

      BigInteger result = MultiplyUnsigned(a.GetInteger(), b.GetInteger(), BigInteger::Base());
      if (a.IsNegative() != b.IsNegative())
        result.SetSign(true);

      return result;
    }

    static BigInteger Multiply(BigInteger const &a, DataT b)
    {
      if (a.IsZero() || b == 0)
        return BigInteger(); // 0 times anything is zero

      vector<DataT> m = ClassicMultiplication::Multiply(a.GetInteger(), b, BigInteger::Base());
      BigInteger result = BigInteger(m, false);
      if (a.IsNegative() != b < 0)
        result.SetSign(true);

      return result;
    }
  };

  // Multiplies Two BigInteger
  BigInteger operator*(BigInteger const &a, BigInteger const &b)
  {
    return BigIntegerMultiplication::Multiply(a, b);
  }

  BigInteger operator*(BigInteger const &a, DataT b)
  {
    return BigIntegerMultiplication::Multiply(a, b);
  }

}

#endif