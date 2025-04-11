#ifndef BURNIKEL_ZIEGLER_DIVISION
#define BURNIKEL_ZIEGLER_DIVISION

#include <stdexcept>
#include <utility>
using namespace std;

#include "../../BigInteger.h"
#include "../classic/CommonAlgorithms.h"
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

    private:
        // Normalize the divisor a so that its most significant limb has the highest bit set.
        // Both dividend and divisor b are shifted left by the same number of bits.
        // Returns the shift amount.
        static int Normalize(BigInteger &a, BigInteger &b)
        {
            if (b.IsZero())
                throw runtime_error("Division by zero in normalization");

            // Determine the number of bits per limb.
            const Int limbBits = sizeof(DataT) * CHAR_BIT; // For DataT = uint16_t, limbBits = 16

            // In little-endian order, the most significant limb is at the end.
            vector<DataT> dvec = b.GetInteger();
            DataT msd = dvec.back();

            // Compute s such that (msd << s) has its most significant bit set.
            Int s = 0;
            // We want: msd << s >= 1 << (limbBits - 1)
            // (Assume msd is nonzero.)
            while (((UInt)msd << s) < (1U << (limbBits - 1)))
            {
                s++;
            }
            // Shift both dividend and divisor left by s bits.
            a = BigInteger(
                CommonAlgorithms::ShiftLeftBits(a.GetInteger(), s));
            b = BigInteger(
                CommonAlgorithms::ShiftLeftBits(b.GetInteger(), s));
            return s;
        }

        // After performing division on normalized operands, denormalize the remainder r
        // by shifting it right by s bits.
        static BigInteger Denormalize(const BigInteger &r, int s)
        {
            return BigInteger(
                CommonAlgorithms::ShiftRightBits(r.GetInteger(), s));
        }

        ///
        /// Recursive Burnikel–Ziegler division algorithm.
        /// Divides a by b, returning {q, r} with q = floor(a/b) and r = a - q*b.
        ///
    private:
        static pair<BigInteger, BigInteger> DivideAndRemainderRecursive(const BigInteger &a, const BigInteger &b)
        {
            if (b.IsZero())
            {
                // Division by zero is undefined.
                throw runtime_error("Division by zero");
            }

            // If the dividend a is smaller than b, the quotient is zero.
            if (a < b)
            {
                return {BigInteger(), a};
            }

            // For small divisors, use the classical division algorithm.
            if (b.size() <= DIVISION_THRESHOLD)
            {
                return BigIntegerClassicDivision::DivideAndRemainder(a, b);
            }

            // --- Recursive case ---
            // Determine block size.
            // A common choice is to let n = ceil(m/2), where m = b.numDigits().
            SizeT m = b.size();
            SizeT n = (m + 1) / 2;

            // We partition the dividend U into blocks of length 2*n (i.e. 2n digits).
            // Let t be the number of blocks:
            SizeT t = (a.size() + (2 * n) - 1) / (2 * n);

            // We will compute the quotient Q and remainder R iteratively.
            BigInteger q;
            BigInteger r;

            // Process blocks from the most-significant block (index t-1) to the least-significant (index 0).
            // (This loop assumes that a is represented so that block 0 corresponds to the least-significant block.
            //  Adjust if your internal representation differs.)
            for (Int i = t; i > 0; i--)
            {
                // Extract the (i-1)-th block of 2*n digits.
                BigInteger block = a.GetBlock(i - 1, 2 * n);

                // Incorporate this block into the current remainder:
                // r = r * (base^(2*n)) + block.
                r = r << (2 * n);
                r = r + block;

                // Now compute the quotient digit for this block.
                // Recursively divide r by b.
                auto qr = DivideAndRemainderRecursive(r, b);
                BigInteger q_i = qr.first; // This is the quotient block (could be more than one digit).
                r = qr.second;             // Remainder after subtracting q_i * V.

                // *** Correction step: ensure r < b by subtracting multiples of b if needed.
                while (r >= b)
                {
                    r = r - b;
                    q_i = q_i + ONE;
                }

                // Accumulate the quotient:
                q = q << (2 * n); // Shift Q left by 2*n digits.
                q = q + q_i;
            }

            // At this point q holds the quotient and r the remainder.
            return {q, r};
        }

    public:
        static pair<BigInteger, BigInteger> DivideAndRemainder(const BigInteger &a, const BigInteger &b)
        {
            if (b.IsZero())
            {
                // Division by zero is undefined.
                throw runtime_error("Division by zero");
            }

            // If the dividend a is smaller than b, the quotient is zero.
            if (a < b)
            {
                return {BigInteger(), a};
            }

            // For small divisors, use the classical division algorithm.
            if (b.size() <= DIVISION_THRESHOLD)
            {
                return BigIntegerClassicDivision::DivideAndRemainder(a, b);
            }

            BigInteger q = a;
            BigInteger r = b;

            int s = Normalize(q, r);
            auto qr = DivideAndRemainderRecursive(q, r);
            q = qr.first;
            r = qr.second;
            r = Denormalize(r, s);
            return {q, r};
        }
    };

} // namespace BigMath

#endif // BURNIKEL_ZIEGLER_DIVISION