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
        // Helper function to compute the greatest common divisor (GCD)
        static DataT gcd(DataT a, DataT b)
        {
            while (b != 0)
            {
                DataT temp = b;
                b = a % b;
                a = temp;
            }
            return a;
        }

    public:
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            if (IsZero(b))
                throw invalid_argument("Division by zero");

            // Copy and trim zeros.
            vector<DataT> normA = a;
            vector<DataT> normB = b;
            TrimZeros(normA);
            TrimZeros(normB);

            int m = normB.size();
            if (m == 0)
                return {vector<DataT>(), normA};

            // Ensure coprimality of b[0] with base.
            DataT b0 = normB[0];
            if (gcd(b0, base) != 1)
            {
                return KnuthDivision::DivideAndRemainder(normA, normB, base);
            }

            // ----- Normalization -----
            DataT d = base / (normB.back() + 1);
            if (d > 1)
            {
                normA = Multiply(normA, d, base);
                normB = Multiply(normB, d, base);
                TrimZeros(normA);
                TrimZeros(normB);
            }

            int n = normA.size();
            int loopCount = n - m + 1;

            // Precompute mu = - (b0^{-1}) mod base.
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
                return KnuthDivision::DivideAndRemainder(normA, normB, base);
            }

            // ----- Montgomery Reduction Loop -----
            vector<DataT> r = normA;
            for (int i = 0; i < loopCount; ++i)
            {
                DataT r0 = r.empty() ? 0 : r[0];
                DataT q_digit = (r0 * mu) % base;

                vector<DataT> qb = Multiply(normB, q_digit, base);
                vector<DataT> shifted_qb(i, 0);
                shifted_qb.insert(shifted_qb.end(), qb.begin(), qb.end());

                r = Add(r, shifted_qb, base);

                if (!r.empty())
                    r.erase(r.begin());
                else
                    r = {0};
            }

            // ----- Final Correction using Division -----
            if (Compare(r, normB) >= 0)
            {
                auto divRes = KnuthDivision::DivideAndRemainder(r, normB, base);
                r = divRes.second;
            }
            TrimZeros(r);

            // ----- De-normalization of the Remainder -----
            if (d > 1)
            {
                auto denorm = KnuthDivision::DivideAndRemainder(r, vector<DataT>{d}, base);
                r = denorm.second;
            }

            // ----- Compute the Quotient and Recompute the True Remainder -----
            vector<DataT> a_minus_r = Subtract(a, r, base);
            auto qr = KnuthDivision::DivideAndRemainder(a_minus_r, b, base);
            TrimZeros(qr.first);
            vector<DataT> quotient = qr.first;

            vector<DataT> prod = Multiply(quotient, b, base);

            vector<DataT> true_remainder = Subtract(a, prod, base);
            TrimZeros(true_remainder);

            return {quotient, true_remainder};
        }
    };
}
#endif