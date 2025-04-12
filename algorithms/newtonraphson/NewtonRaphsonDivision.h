#ifndef NEWTONRAPHSON_DIVISION
#define NEWTONRAPHSON_DIVISION

#include <vector>
#include <cassert>

#include "../../BigIntegerBuilder.h"
#include "../../BigIntegerUtil.h"
#include "../classic/CommonAlgorithms.h"
#include "../classic/ClassicAddition.h"
#include "../classic/ClassicSubtraction.h"

using namespace std;

namespace BigMath
{
    class NewtonRaphsonDivision
    {
        // Helper: determine a normalization shift.
        // For example, if our BigInteger is internally stored in base 2^16,
        // we want to shift the divisor so that its most significant digit
        // uses as many bits as possible.
        static SizeT getNormalizationShift(vector<DataT> const &d)
        {
            // Assume BigInteger has a method size() that returns the count of base–2^16 digits,
            // and that digits() returns a vector in little‑endian order.
            // We normalize by ensuring that the highest digit (in the most–significant limb) is large.
            DataT highest = d.back(); // access the most significant digit
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

        static long double toDouble(vector<DataT> const &x, BaseT B)
        {
            // Assume Base() returns 65536.0L (i.e. 2^16), and the integer is stored
            // in little-endian order in the vector theInteger.
            const long double base = B;
            const SizeT totalLimbs = x.size();

            // Choose K such that we use only K most-significant limbs.
            // For example, if long double has an 80-bit significand and each limb holds 16 bits,
            // then K = 5 is a reasonable choice.
            const SizeT K = (totalLimbs < 5 ? totalLimbs : 5);

            // Compute the fraction from the top K limbs.
            // Since theInteger is stored in little-endian order, the most-significant limbs
            // are at the end.
            long double fraction = 0.0L;
            vector<DataT> theInteger = x;
            // Start from index = totalLimbs - K up to totalLimbs - 1.
            for (size_t i = totalLimbs - K; i < totalLimbs; i++)
            {
                fraction = fraction * base + static_cast<long double>(theInteger[i]);
            }

            // Each limb represents 16 bits, so the total number of bits represented by the skipped limbs is:
            const Long limbBits = 16;
            Long computedExponent = static_cast<Long>((totalLimbs - K) * limbBits);

            // Check against the maximum exponent for long double.
            // numeric_limits<long double>::max_exponent gives the maximum exponent value,
            // i.e. the power of 2 below which numbers are representable.
            Long maxExp = numeric_limits<long double>::max_exponent;

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

            return result;
        }

        // Assume toDouble(BigInteger const &) is implemented elsewhere.
        static vector<DataT> initialGuess(vector<DataT> const &dNorm, BaseT base, SizeT k)
        {
            // Convert dNorm to long double (this gives an approximate magnitude).
            long double dDouble = toDouble(dNorm, base);
            // Compute the reciprocal in double precision.
            long double xDouble = 1.0L / dDouble;

            // Instead of computing R = ONE << k and then toDouble(R),
            // use ldexpl to scale xDouble by 2^k directly.
            long double scaled = ldexpl(xDouble, static_cast<int>(k));

            // Optionally, you can add a check to clamp the result if necessary.
            if (scaled == 0.0L)
            {
                // If underflow occurs, clamp to the smallest normalized value.
                scaled = numeric_limits<long double>::min();
            }
            else if (!isfinite(scaled))
            {
                // If the value is too high, clamp to the maximum finite value.
                scaled = numeric_limits<long double>::max();
            }

            // Construct and return a BigInteger from the scaled value.
            vector<DataT> x = BigIntegerBuilder::VectorFrom(scaled).first;
            return x;
        }

    public:
        // Newton–Raphson division: Returns a pair {q, r} such that q = floor(a/b) and r = a - q*b.
        // 'iterations' controls the number of Newton iterations.
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(vector<DataT> const &a, vector<DataT> const &b, BaseT base, int iterations = 5)
        {
            // Normalize divisor and dividend.
            SizeT shift = getNormalizationShift(b);
            vector<DataT> dNorm = CommonAlgorithms::ShiftLeft(b, shift);
            vector<DataT> nNorm = CommonAlgorithms::ShiftLeft(a, shift);
            const vector<DataT> ONE(1, 1); // size 1, value 1

            // Choose fixed-point scaling precision. For instance, if dNorm has L limbs (each being 16 bits):
            SizeT L = dNorm.size();
            SizeT k = 2 * L * 16;
            vector<DataT> R_val = CommonAlgorithms::ShiftLeft(ONE, k);

            // Compute an initial guess for x ≈ r/dNorm.
            vector<DataT> x = initialGuess(dNorm, base, k);

            // Newton–Raphson iterations: x_{i+1} = (x_i * (2R - dNorm * x_i)) / R.
            for (int i = 0; i < iterations; i++)
            {
                vector<DataT> prod = BigIntegerMultiplication::MultiplyUnsigned(dNorm, x, base);                        // Compute dNorm * x.
                vector<DataT> rshift = ClassicMultiplication::Multiply(R_val, (DataT)2, base);
                vector<DataT> correction = ClassicSubtraction::Subtract(rshift, prod, base); // Compute 2R - dNorm*x.
                vector<DataT> xcor = BigIntegerMultiplication::MultiplyUnsigned(x, correction, base);
                x = CommonAlgorithms::ShiftRight(xcor, k);                          // Update: dividing by R via shifting.
            }

            // Compute quotient. With x ≈ R/dNorm, the normalized quotient is:
            vector<DataT> xn = BigIntegerMultiplication::MultiplyUnsigned(x, nNorm, base);
            vector<DataT> q = CommonAlgorithms::ShiftLeft(xn, k);

            // Correct Q in case of slight over- or under-estimation.
            vector<DataT> qb = BigIntegerMultiplication::MultiplyUnsigned(q, b, base);
            while (qb > a) {
                q = ClassicSubtraction::Subtract(q, ONE, base);
                qb = BigIntegerMultiplication::MultiplyUnsigned(q, b, base);
            }

            vector<DataT> Qp1 = ClassicAddition::Add(q, ONE, base);
            qb = BigIntegerMultiplication::MultiplyUnsigned(Qp1, b, base);
            while (qb <= a)
            {
                q = ClassicAddition::Add(q, ONE, base);
                qb = BigIntegerMultiplication::MultiplyUnsigned(q, b, base);
            }

            // Once q is computed, the remainder is calculated as:
            qb = BigIntegerMultiplication::MultiplyUnsigned(q, b, base);
            vector<DataT> r = ClassicSubtraction::Subtract(a, qb, base);

            return {q, r};
        }
    };
}
#endif