#ifndef MONTGOMERY_DIVISION
#define MONTGOMERY_DIVISION

#include <vector>
#include <cassert>
#include <cmath>

#include "../../common/Builder.h"
#include "../../common/Util.h"
#include "../../ops/IO.h"
#include "../Multiplication.h"
#include "../Shift.h"
#include "../Addition.h"
#include "../Subtraction.h"
#include "KnuthDivision.h"

using namespace std;

namespace BigMath
{
    class MontgomeryDivision
    {
    public:
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            if (IsZero(b))
            {
                throw invalid_argument("Division by zero");
            }

            int m = b.size();
            if (m == 0)
            {
                return {vector<DataT>(), a};
            }

            DataT b0 = b[0];
            if (gcd(b0, base) != 1)
            {
                // Fall back to classic division if b0 and base are not coprime.
                return KnuthDivision::DivideAndRemainder(a, b, base);
            }

            // The Montgomery reduction algorithm as implemented here is efficient only if the dividend
            // is roughly of size m or m+1 digits (i.e. in the range [b, b * base)).
            // If a is much larger, then the final subtraction loop might run unacceptably long.
            if (a.size() > m + 1)
            {
                return KnuthDivision::DivideAndRemainder(a, b, base);
            }

            // ----- Normalization -----
            vector<DataT> A, B;

            // Multiply both norm_a and norm_b by d to ensure that the most significant digit
            // of norm_b (i.e. norm_b.back()) is large (at least base/2).
            DataT d = base / (B.back() + 1);
            if (d > 1)
            {
                A = Multiply(A, d, base);
                B = Multiply(B, d, base);
            }
            else
            {
                A = a;
                B = b;
            }

            // ----- Precompute mu -----
            // Compute mu = - (b0^{-1}) mod base.
            DataT mu = 0;
            bool found = false;
            for (SizeT i = 1; i < base; ++i)
            {
                if ((b0 * i) % base == 1)
                {
                    mu = (base - i) % base;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                // Should never get here under the coprimality condition, but fallback if needed.
                return KnuthDivision::DivideAndRemainder(a, b, base);
            }

            // ----- Montgomery Reduction Loop -----
            // The idea is to “eliminate” the lower m digits one by one.
            // This loop requires that norm_a (the normalized dividend) is at most m+1 digits long.
            vector<DataT> r = A;
            for (int i = 0; i < m; ++i)
            {
                // Access the least significant digit.
                DataT r0 = (r.empty() ? 0 : r[0]);
                // Compute the correction digit q = (r0 * mu) mod base.
                DataT q = (r0 * mu) % base;

                // Multiply the entire divisor by the single digit q.
                vector<DataT> qb = Multiply(B, q, base);
                // In the Montgomery reduction formula this term would be shifted by the digit index.
                // However, since we always operate on the least significant digit,
                // no extra shift is needed.
                AddTo(r, qb, base);

                // "Divide" r by base by removing the least significant digit.
                if (!r.empty())
                    r.erase(r.begin());
                else
                    r = {0};
            }

            // ----- Final Correction -----
            // At this point, r should be less than norm_b modulo base^m.
            // If it is not, subtract norm_b (possibly several times) to bring it into the proper range.
            // Because norm_b is now normalized and r has at most m digits, the following loop will not iterate many times.
            while (Compare(r, B) >= 0)
            {
                SubtractFrom(r, B, base);
            }
            TrimZeros(r);

            // ----- De-normalization -----
            // Since we multiplied the numbers by d at the start, we must de-normalize the remainder.
            // For small d (a single digit), we can recover the true remainder by performing a division by d.
            if (d > 1)
            {
                auto denorm = KnuthDivision::DivideAndRemainder(r, d, base);
                r = denorm.second;
            }

            // ----- Compute Quotient -----
            // The quotient is computed via the identity: a = q * b + r, that is, q = (a - r) / b.
            vector<DataT> a_minus_r = Subtract(a, r, base);
            auto q = KnuthDivision::Divide(a_minus_r, b, base);
            TrimZeros(q);

            return {q, r};
        }
    };
}
#endif