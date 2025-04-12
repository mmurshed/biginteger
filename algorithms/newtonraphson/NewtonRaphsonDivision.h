#ifndef NEWTONRAPHSON_DIVISION
#define NEWTONRAPHSON_DIVISION

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "../../BigIntegerBuilder.h"
#include "../../BigIntegerUtil.h"
#include "../../ops/BigIntegerIO.h"
#include "../classic/CommonAlgorithms.h"
#include "../classic/ClassicAddition.h"
#include "../classic/ClassicSubtraction.h"
#include "../classic/ClassicDivision.h"

using namespace std;

namespace BigMath
{
    class NewtonRaphsonDivision
    {
        // Helper: determine a normalization shift.
        // For example, if our BigInteger is internally stored in base 2^16,
        // we want to shift the divisor so that its most significant digit
        // uses as many bits as possible.
        static SizeT getNormalizationShift(const vector<DataT> &d, BaseT base)
        {
            if (d.empty()) return 0;
        
            DataT highest = d.back(); // Most significant digit
            SizeT shift = 0;
        
            // Find how many bit shifts (left) are needed to bring highest close to base
            // We'll shift until the digit is >= base / 2 (top half of the range)
            // So that its top bit(s) are significant when viewed in binary
            while (static_cast<ULong>(highest) < static_cast<ULong>(base) / 2)
            {
                highest <<= 1;
                shift++;
            }

            cerr << "[getNormalizationShift] shift: " << shift << endl;
        
            return shift;
        }

        static long double toDouble(vector<DataT> const &x, long double base, long double log2base)
        {
            // Assume Base() returns 65536.0L (i.e. 2^16), and the integer is stored
            // in little-endian order in the vector theInteger.
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
            Long computedExponent = static_cast<Long>((totalLimbs - K) * log2base);

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

            cerr << "[toDouble] result: " << result << endl;

            return result;
        }

        // Assume toDouble(BigInteger const &) is implemented elsewhere.
        static vector<DataT> initialGuess(vector<DataT> const &dNorm, long double base, long double log2base, SizeT k)
        {
            // Convert dNorm to long double (this gives an approximate magnitude).
            long double dDouble = toDouble(dNorm, base, log2base);
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

            cerr << "[initialGuess] x: " << x << endl;
            return x;
        }

    public:
        // Newton–Raphson division: Returns a pair {q, r} such that q = floor(a/b) and r = a - q*b.
        // 'iterations' controls the number of Newton iterations.
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(vector<DataT> const &a, vector<DataT> const &b, BaseT base, int iterations = 5)
        {
            // Normalize divisor and dividend.
            SizeT shift = getNormalizationShift(b, base);
            vector<DataT> dNorm = CommonAlgorithms::ShiftLeft(b, shift);
            vector<DataT> nNorm = CommonAlgorithms::ShiftLeft(a, shift);
            const vector<DataT> ONE(1, 1); // size 1, value 1

            // Choose fixed-point scaling precision. For instance, if dNorm has L limbs (each being 16 bits):
            SizeT L = dNorm.size();
            long double bitsPerDigit = log2(static_cast<long double>(base));
            SizeT k = static_cast<SizeT>(2 * L * bitsPerDigit);
            cerr << "[Newton] k: " << k << endl;
            // Shift R by k bits = shift left by k / bitsPerDigit digits
            SizeT digitShift = k / bitsPerDigit;
            cerr << "[Newton] digitShift: " << digitShift << endl;
            vector<DataT> R = CommonAlgorithms::ShiftLeft(ONE, digitShift);
            cerr << "[Newton] R: " << R << endl;

            // vector<DataT> twoR = ClassicMultiplication::Multiply(R_val, (DataT)2, base);
            vector<DataT> twoR = CommonAlgorithms::ShiftLeft(R, (SizeT)1); // Compute 2R.
            cerr << "[Newton] 2R: " << twoR << endl;

            // Compute an initial guess for x ≈ r/dNorm.
            vector<DataT> x = initialGuess(dNorm, base, bitsPerDigit, k);
            // vector<DataT> x = ClassicDivision::DivideAndRemainder(twoR, dNorm, base).first;

            // Newton–Raphson iterations: x_{i+1} = (x_i * (2R - dNorm * x_i)) / R.
            for (int i = 0; i < iterations; i++)
            {
                cerr << "[Newton] Iteration " << i + 1 << endl;

                vector<DataT> prod = BigIntegerMultiplication::MultiplyUnsigned(dNorm, x, base);                        // Compute dNorm * x.
                cerr << "[Newton] dNorm*x: " << prod << endl;

                vector<DataT> correction = ClassicSubtraction::Subtract(twoR, prod, base); // Compute 2R - dNorm*x.
                cerr << "[Newton] 2R - dNorm*x: " << correction << endl;
                
                vector<DataT> xcor = BigIntegerMultiplication::MultiplyUnsigned(x, correction, base);
                cerr << "[Newton] x_i * (2R - dNorm*x): " << xcor << endl;

                x = CommonAlgorithms::ShiftRight(xcor, digitShift);                          // Update: dividing by R via shifting.
                cerr << "[Newton] (x_i * (2R - dNorm*x)) / R: " << x << endl;
            }

            // Compute quotient. With x ≈ R/dNorm, the normalized quotient is:
            cerr << "[Quotient] Computing a * x and shifting right by " << k << " bits" << endl;
            vector<DataT> xn = BigIntegerMultiplication::MultiplyUnsigned(x, nNorm, base);
            vector<DataT> q = CommonAlgorithms::ShiftRight(xn, digitShift);
            cerr << "[Quotient] q " << q << endl;
            
            cerr << "[Correction] Starting downward correction if q*b > a" << endl;
            vector<DataT> qb = BigIntegerMultiplication::MultiplyUnsigned(q, b, base);
            while (BigIntegerComparator::Compare(qb, a) > 0)
            {
                ClassicSubtraction::SubtractFrom(q, 1, base);
                qb = BigIntegerMultiplication::MultiplyUnsigned(q, b, base);
            }
            cerr << "[Quotient] q " << q << endl;
            
            cerr << "[Correction] Starting upward correction if (q+1)*b <= a" << endl;
            
            do
            {
                ClassicAddition::AddTo(q, 1, base);
                qb = BigIntegerMultiplication::MultiplyUnsigned(q, b, base);
                if (BigIntegerComparator::Compare(qb, a) > 0)
                    break;                
            } while (true);
            ClassicSubtraction::SubtractFrom(q, 1, base);

            cerr << "[Quotient] q " << q << endl;

            // Once q is computed, the remainder is calculated as:
            qb = BigIntegerMultiplication::MultiplyUnsigned(q, b, base);
            vector<DataT> r = ClassicSubtraction::Subtract(a, qb, base);

            cerr << "[Result] Final quotient: " << q << ", remainder: " << r << endl;
            
            return {q, r};
        }
    };
}
#endif