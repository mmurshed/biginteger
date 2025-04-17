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
#include "ClassicDivision.h"
#include "KnuthDivision.h"

using namespace std;

namespace BigMath
{
    class MontgomeryDivision
    {
        static DataT mod_inverse(DataT a, DataT base) {
            DataT t = 0, newt = 1, r = base, newr = a;
            while (newr != 0) {
                DataT quotient = r / newr;
                std::tie(t, newt) = std::make_tuple(newt, t - quotient * newt);
                std::tie(r, newr) = std::make_tuple(newr, r - quotient * newr);
            }
            if (r > 1) throw runtime_error("No inverse");
            return (t < 0) ? t + base : (base - t) % base;
        }
    public:
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            if (IsZero(b))
                throw invalid_argument("Division by zero");

            vector<DataT> A = a, B = b;
            TrimZeros(B);
            int m = B.size();
            if (m == 0)
                return {vector<DataT>(), a};

            DataT b0 = B[0];
            if (gcd(b0, base) != 1)
                return ClassicDivision::DivideAndRemainder(a, b, base);

            // Normalization
            DataT d = (B.back() + 1 < base) ? base / (B.back() + 1) : 1;
            if (d > 1)
            {
                A = Multiply(A, d, base);
                B = Multiply(B, d, base);
                TrimZeros(B);
                m = B.size();
            }

            // Compute mu using Extended Euclidean Algorithm
            DataT mu = mod_inverse(B[0], base);
            mu = (base - mu) % base;

            // Montgomery Reduction with Index Tracking
            vector<DataT> r = A;
            size_t start = 0;
            for (int i = 0; i < m; ++i)
            {
                DataT r0 = (start < r.size()) ? r[start] : 0;
                DataT q = (r0 * mu) % base;
                vector<DataT> qb = Multiply(B, q, base);
                AddToAt(r, qb, start, base);
                start++;
            }
            r = vector<DataT>(r.begin() + start, r.end());
            TrimZeros(r);

            // Final Subtraction (At Most Once)
            if (Compare(r, B) >= 0)
                SubtractFrom(r, B, base);

            // Denormalize Remainder
            if (d > 1)
            {
                DataT carry = 0;
                vector<DataT> denorm;
                for (auto it = r.rbegin(); it != r.rend(); ++it)
                {
                    DataT temp = carry * base + *it;
                    denorm.push_back(temp / d);
                    carry = temp % d;
                }
                reverse(denorm.begin(), denorm.end());
                r = denorm;
                TrimZeros(r);
            }

            // Compute Quotient
            vector<DataT> a_minus_r = Subtract(a, r, base);
            auto q = KnuthDivision::Divide(a_minus_r, b, base);
            TrimZeros(q);

            return {q, r};
        }
    };
}
#endif