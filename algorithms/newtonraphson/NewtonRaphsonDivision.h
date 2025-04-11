#ifndef NEWTONRAPHSON_DIVISION
#define NEWTONRAPHSON_DIVISION

#include <vector>
#include <cassert>

#include "../../BigInteger.h"
#include "../../BigIntegerBuilder.h"
#include "../BigIntegerUtil.h"
#include "../../ops/BigIntegerAddition.h"
#include "../../ops/BigIntegerSubtraction.h"
#include "../../ops/BigIntegerMultiplication.h"
#include "../../ops/BigIntegerComparison.h"
#include "../../ops/BigIntegerShift.h"

using namespace std;

namespace BigMath
{
    class NewtonRaphsonDivision
    {
        // Helper: determine a normalization shift.
        // For example, if our BigInteger is internally stored in base 2^16,
        // we want to shift the divisor so that its most significant digit
        // uses as many bits as possible.
        static SizeT getNormalizationShift(BigInteger const &d)
        {
            // Assume BigInteger has a method size() that returns the count of base–2^16 digits,
            // and that digits() returns a vector in little‑endian order.
            // We normalize by ensuring that the highest digit (in the most–significant limb) is large.
            DataT highest = d.GetInteger().back(); // access the most significant digit
            SizeT shift = 0;
            // For base 2^16, the maximum value of a digit is 65535.
            // We shift until the high bit of the highest digit is set.
            while (highest < (1 << 15))
            { // less than 32768 means top bit not set
                highest <<= 1;
                shift++;
            }
            return shift;
        }

        static long double toDouble(BigInteger const &x)
        {
            // Assume Base() returns 65536.0L (i.e. 2^16), and the integer is stored
            // in little-endian order in the vector theInteger.
            const long double base = BigInteger::Base();
            const SizeT totalLimbs = x.size();

            // Choose K such that we use only K most-significant limbs.
            // For example, if long double has an 80-bit significand and each limb holds 16 bits,
            // then K = 5 is a reasonable choice.
            const SizeT K = (totalLimbs < 5 ? totalLimbs : 5);

            // Compute the fraction from the top K limbs.
            // Since theInteger is stored in little-endian order, the most-significant limbs
            // are at the end.
            long double fraction = 0.0L;
            vector<DataT> theInteger = x.GetInteger();
            // Start from index = totalLimbs - K up to totalLimbs - 1.
            for (size_t i = totalLimbs - K; i < totalLimbs; i++)
            {
                fraction = fraction * base + static_cast<long double>(theInteger[i]);
            }

            // Each limb represents 16 bits, so the total number of bits represented by the skipped limbs is:
            const Long limbBits = 16;
            Long computedExponent = static_cast<Long>((totalLimbs - K) * limbBits);

            // Check against the maximum exponent for long double.
            // std::numeric_limits<long double>::max_exponent gives the maximum exponent value,
            // i.e. the power of 2 below which numbers are representable.
            Long maxExp = std::numeric_limits<long double>::max_exponent;

            long double result;
            if (computedExponent > maxExp)
            {
                // The exponent is too high; clamp the result.
                // We use the maximum finite long double value.
                result = (fraction < 0.0L ? -numeric_limits<long double>::max()
                                          : numeric_limits<long double>::max());
            }
            else
            {
                // Use ldexpl to compute fraction * 2^(computedExponent).
                result = ldexpl(fraction, computedExponent);
            }

            // Apply the sign.
            if (x.IsNegative())
                result = -result;

            return result;
        }

        // Assume toDouble(BigInteger const &) is implemented elsewhere.
        static BigInteger initialGuess(BigInteger const &dNorm, SizeT k)
        {
            // Convert dNorm to long double (this gives an approximate magnitude).
            long double dDouble = toDouble(dNorm);
            // Compute the reciprocal in double precision.
            long double xDouble = 1.0L / dDouble;

            // Instead of computing R = ONE << k and then toDouble(R),
            // use ldexpl to scale xDouble by 2^k directly.
            long double scaled = ldexpl(xDouble, static_cast<int>(k));

            // Optionally, you can add a check to clamp the result if necessary.
            if (scaled == 0.0L)
            {
                // If underflow occurs, clamp to the smallest normalized value.
                scaled = std::numeric_limits<long double>::min();
            }
            else if (!std::isfinite(scaled))
            {
                // If the value is too high, clamp to the maximum finite value.
                scaled = std::numeric_limits<long double>::max();
            }

            // Construct and return a BigInteger from the scaled value.
            BigInteger x = BigIntegerBuilder::From(scaled);
            return x;
        }

    public:
        // Newton–Raphson division: Returns a pair {Q, R} such that Q = floor(N/D) and R = N - Q*D.
        // 'iterations' controls the number of Newton iterations.
        static pair<BigInteger, BigInteger> DivideAndRemainder(BigInteger const &N, BigInteger const &D, int iterations = 5)
        {
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
            for (int i = 0; i < iterations; i++)
            {
                BigInteger prod = dNorm * x;                        // Compute dNorm * x.
                BigInteger correction = (R_val << (SizeT)1) - prod; // Compute 2R - dNorm*x.
                x = (x * correction) >> k;                          // Update: dividing by R via shifting.
            }

            // Compute quotient. With x ≈ R/dNorm, the normalized quotient is:
            BigInteger Q = (x * nNorm) >> k;

            // Correct Q in case of slight over- or under-estimation.
            while (Q * D > N)
                Q = Q - ONE;

            BigInteger Qp1 = Q + ONE;
            while (Qp1 * D <= N)
            {
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