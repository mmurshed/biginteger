#ifndef BURNIKEL_ZIEGLER_DIVISION
#define BURNIKEL_ZIEGLER_DIVISION

#include <stdexcept>
#include <utility>
using namespace std;

#include "../../BigInteger.h"
#include "../../ops/BigIntegerAddition.h"
#include "../../ops/BigIntegerShift.h"
#include "../../ops/BigIntegerClassicDivision.h"
#include "../../ops/BigIntegerComparison.h"

namespace BigMath
{
    class BurnikelZieglerDivision
    {
        // A threshold (in digits or limbs) below which we use classical long division.
        // You may choose this value based on experiments.
        static const SizeT DIVISION_THRESHOLD = 4;        

        ///
        /// Recursive Burnikel–Ziegler division algorithm.
        /// Divides U by V, returning {Q, R} with Q = floor(U/V) and R = U - Q*V.
        /// This is a simplified skeleton version.
        ///
        /// Requirements on BigInteger:
        ///   - numDigits(): returns the number of limbs/digits.
        ///   - getBlock(index, blockSize): extracts a block (limb group) as a BigInteger.
        ///       (Blocks are numbered from 0 for the least-significant; you may need to adjust
        ///        for the algorithm below which processes from the most-significant block.)
        ///   - shiftLeft(digitCount): multiplies the number by base^(digitCount).
        ///   - operator+(const BigInteger &)
        ///
        public:
        static pair<BigInteger, BigInteger> DivideAndRemainder(const BigInteger &U, const BigInteger &V)
        {
            if (V.IsZero())
            {
                // Division by zero is undefined.
                throw runtime_error("Division by zero");
            }

            // If the dividend U is smaller than V, the quotient is zero.
            if (U < V)
            {
                return {BigInteger(), U};
            }

            // For small divisors, use the classical division algorithm.
            if (V.size() <= DIVISION_THRESHOLD)
            {
                return BigIntegerClassicDivision::DivideAndRemainder(U, V);
            }

            // --- Recursive case ---
            // Determine block size.
            // A common choice is to let n = ceil(m/2), where m = V.numDigits().
            SizeT m = V.size();
            SizeT n = (m + 1) / 2;

            // We partition the dividend U into blocks of length 2*n (i.e. 2n digits).
            // Let t be the number of blocks:
            SizeT t = (U.size() + (2 * n) - 1) / (2 * n);

            // We will compute the quotient Q and remainder R iteratively.
            BigInteger Q(0);
            BigInteger R(0);

            // Process blocks from the most-significant block (index t-1) to the least-significant (index 0).
            // (This loop assumes that U is represented so that block 0 corresponds to the least-significant block.
            //  Adjust if your internal representation differs.)
            for (Int i = t; i > 0; i--)
            {
                // Extract the (i-1)-th block of 2*n digits.
                BigInteger block = U.GetBlock(i - 1, 2 * n);

                // Incorporate this block into the current remainder:
                // R = R * (base^(2*n)) + block.
                R = R << 2 * n;
                R = R + block;

                // Now compute the quotient digit for this block.
                // Recursively divide R by V.
                auto qr = DivideAndRemainder(R, V);
                BigInteger q_i = qr.first; // This is the quotient block (could be more than one digit).
                R = qr.second;             // Remainder after subtracting q_i * V.

                // Accumulate the quotient:
                Q = Q << 2 * n; // Shift Q left by 2*n digits.
                Q = Q + q_i;
            }

            // At this point Q holds the quotient and R the remainder.
            return {Q, R};
        }
    };

} // namespace BigMath

#endif // BURNIKEL_ZIEGLER_DIVISION