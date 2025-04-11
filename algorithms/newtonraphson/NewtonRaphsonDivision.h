#include <vector>
#include <cassert>

#include "../../BigInteger.h"
#include "../BigIntegerUtil.h"
#include "../../ops/BigIntegerAddition.h"
#include "../../ops/BigIntegerSubtraction.h"
#include "../../ops/BigIntegerMultiplication.h"
#include "../../ops/BigIntegerComparison.h"
#include "../../ops/BigIntegerShift.h"

using namespace std;

#ifndef NEWTONRAPHSON_DIVISION
#define NEWTONRAPHSON_DIVISION


namespace BigMath
{  
    class NewtonRaphsonDivision
    {
  // Helper: determine a normalization shift.
  // For example, if our BigInteger is internally stored in base 2^16,
  // we want to shift the divisor so that its most significant digit
  // uses as many bits as possible.
  static SizeT getNormalizationShift(BigInteger const& d) {
      // Assume BigInteger has a method size() that returns the count of base–2^16 digits,
      // and that digits() returns a vector in little‑endian order.
      // We normalize by ensuring that the highest digit (in the most–significant limb) is large.
      DataT highest = d.GetInteger().back();  // access the most significant digit
      SizeT shift = 0;
      // For base 2^16, the maximum value of a digit is 65535.
      // We shift until the high bit of the highest digit is set.
      while (highest < (1 << 15)) {  // less than 32768 means top bit not set
          highest <<= 1;
          shift++;
      }
      return shift;
  }

  // Helper: obtain an initial guess for the reciprocal using a double approximation.
  // The guess is computed for fixed point reciprocal with scaling factor R = 2^k.
  static BigInteger initialGuess(BigInteger const& dNorm, SizeT k) {
      // Convert dNorm to double.
      long double dDouble = dNorm.ToDouble();  // You must implement toDouble() appropriately.
      long double xDouble = 1.0 / dDouble;
      // Create fixed-point scaling factor R = 2^k.
      BigInteger R = ONE << k;  // assumes << operator works
      // Multiply xDouble by R and convert to BigInteger.
      // For example, if BigInteger can be constructed from a double, do:
      BigInteger x = BigInteger(xDouble * R.ToDouble());
      return x;
  }

public:
  // Newton–Raphson division: Returns a pair {Q, R} such that Q = floor(N/D) and R = N - Q*D.
  // 'iterations' controls the number of Newton iterations.
  static pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const& N, BigInteger const& D, int iterations = 5) {
      
      // Normalize divisor and dividend.
      SizeT shift = getNormalizationShift(D);
      BigInteger dNorm = D << shift;
      BigInteger nNorm = N << shift;
      
      // Choose fixed-point scaling precision. For instance, if dNorm has L limbs (each being 16 bits):
      SizeT L = dNorm.size();
      SizeT k = 2 * L * 16;
      BigInteger R_val = ONE << k;
      
      // Compute an initial guess for x ≈ R/dNorm.
      BigInteger x = initialGuess(dNorm, k);
      
      // Newton–Raphson iterations: x_{i+1} = (x_i * (2R - dNorm * x_i)) / R.
      for (int i = 0; i < iterations; i++) {
          BigInteger prod = dNorm * x;            // Compute dNorm * x.
          BigInteger correction = (R_val << (SizeT)1) - prod; // Compute 2R - dNorm*x.
          x = (x * correction) >> k;                // Update: dividing by R via shifting.
      }
      
      // Compute quotient. With x ≈ R/dNorm, the normalized quotient is:
      BigInteger Q = (x * nNorm) >> k;
      
      // Correct Q in case of slight over- or under-estimation.
      while (Q * D > N)
          Q = Q - ONE;

      BigInteger Qp1 = Q + ONE;
      while (Qp1 * D <= N) {
          Q = Qp1;
          Qp1 = Q + ONE;
      }
      
      // Once Q is computed, the remainder is calculated as:
      BigInteger Remainder = N - (Q * D);
      
      return {Q, Remainder};
  }
    };
}
#endif