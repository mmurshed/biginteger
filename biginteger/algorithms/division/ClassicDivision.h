/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef CLASSIC_DIVISION
#define CLASSIC_DIVISION

#include <cstddef>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;

#include "../../common/Util.h"
#include "../Shift.h"
#include "../Addition.h"
#include "../Subtraction.h"
#include "../multiplication/ClassicMultiplication.h"

namespace BigMath
{
    class ClassicDivision
    {
    public:
        static void DivideTo(vector<DataT> &u, DataT d, BaseT base)
        {
            Divide(u, 0, u.size() - 1, d, u, 0, u.size() - 1, base);
        }

        static void DivideTo(
            vector<DataT> &u, SizeT uStart, SizeT uEnd,
            DataT d,
            BaseT base)
        {
            Divide(u, uStart, uEnd, d, u, uStart, uEnd, base);
        }

        static vector<DataT> Divide(
            vector<DataT> const &u,
            DataT d,
            BaseT base)
        {
            SizeT n = (SizeT)u.size();
            // Divide (u_n−1 . . . u_1 u_0)_b by d.
            vector<DataT> w(n);

            Divide(u, 0, u.size() - 1, d, w, 0, w.size() - 1, base);

            return w;
        }

        static void Divide(
            vector<DataT> const &u, SizeT uStart, SizeT uEnd,
            DataT d,
            vector<DataT> &w, SizeT wStart, SizeT wEnd,
            BaseT base)
        {
            // Number of limbs in [uStart..uEnd]
            SizeT n = (uEnd >= uStart) ? (uEnd - uStart + 1) : 0;
            if (n == 0)
            {
                // Nothing to divide
                w.clear();
                return;
            }

            // Make sure 'w' can hold at least n limbs for the quotient
            SizeT neededSize = wStart + n;
            if (w.size() < neededSize)
                w.resize(neededSize, 0);

            // 'r' is the running remainder.
            ULong r = 0;

            // Process from the highest index (uEnd) down to the lowest (uStart).
            //   For little-endian, u[uEnd] is the most significant limb.
            for (Int i = (Int)uEnd; i >= (Int)uStart; i--)
            {
                ULong val = r * base + u[i]; // Combine remainder + next limb
                DataT q = (DataT)(val / d);  // Quotient digit
                r = val % d;                 // Remainder

                // The quotient digit for this place goes at wPos (also in little-endian).
                SizeT wPos = wStart + (i - uStart);
                if (wPos <= wEnd)
                    w[wPos] = q;
                else
                    w.push_back(q);
            }

            // Remove any leading zero-limbs in the quotient.
            TrimZeros(w);
        }

        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
            vector<DataT> const &u,
            DataT d,
            BaseT base)
        {
            SizeT n = (SizeT)u.size();
            // Divide (u_n−1 . . . u_1 u_0)_b by d.
            vector<DataT> w(n);
            vector<DataT> v(n);

            DivideAndRemainder(
                u, 0, u.size() - 1,
                d,
                w, 0, w.size() - 1,
                v, 0, v.size() - 1,
                base);

            return {w, v};
        }

        static void DivideAndRemainder(
            vector<DataT> const &u, SizeT uStart, SizeT uEnd,
            DataT d,
            vector<DataT> &w, SizeT wStart, SizeT wEnd,
            vector<DataT> &v, SizeT vStart, SizeT vEnd,
            BaseT base)
        {
            // Number of limbs in [uStart..uEnd]
            SizeT n = (uEnd >= uStart) ? (uEnd - uStart + 1) : 0;
            if (n == 0)
            {
                // Nothing to divide
                return;
            }

            // Make sure 'w' can hold at least n limbs for the quotient
            SizeT neededSize = wStart + n;
            if (w.size() < neededSize)
                w.resize(neededSize, 0);

            // 'r' is the running remainder.
            ULong r = 0;

            // Process from the highest index (uEnd) down to the lowest (uStart).
            //   For little-endian, u[uEnd] is the most significant limb.
            for (Int i = (Int)uEnd; i >= (Int)uStart; i--)
            {
                ULong val = r * base + u[i]; // Combine remainder + next limb
                DataT q = (DataT)(val / d);  // Quotient digit
                r = val % d;                 // Remainder

                // The quotient digit for this place goes at wPos (also in little-endian).
                SizeT wPos = wStart + (i - uStart);
                if (wPos <= wEnd)
                    w[wPos] = q;
                else
                    w.push_back(q);

                // The quotient digit for this place goes at wPos (also in little-endian).
                SizeT vPos = vStart + (i - uStart);
                if (vPos <= vEnd)
                    v[vPos] = r;
                else
                    v.push_back(r);
            }

            // Remove any leading zero-limbs in the quotient.
            TrimZeros(w);
            TrimZeros(v);
        }

        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            if (IsZero(b))
            {
                throw invalid_argument("Division by zero");
            }

            int cmp = Compare(a, b);
            if (cmp < 0)
            {
                return {vector<DataT>(), a};
            }
            else if (cmp == 0)
            {
                return {vector<DataT>{1}, vector<DataT>()};
            }

            vector<DataT> remainder = a;
            vector<DataT> quotient;

            int max_d = a.size() - b.size();
            for (int d = max_d; d >= 0; --d)
            {
                vector<DataT> shifted_divisor;
                shifted_divisor.insert(shifted_divisor.end(), d, 0);
                shifted_divisor.insert(shifted_divisor.end(), b.begin(), b.end());

                if (Compare(shifted_divisor, remainder) > 0)
                {
                    if (quotient.size() > static_cast<size_t>(d))
                    {
                        quotient[d] = 0;
                    }
                    continue;
                }

                DataT q_low = 1, q_high = base - 1;
                DataT q = 0;
                while (q_low <= q_high)
                {
                    DataT candidate = (q_low + q_high) / 2;
                    vector<DataT> product = ClassicMultiplication::Multiply(shifted_divisor, candidate, base);

                    int cmp = Compare(product, remainder);
                    if (cmp <= 0)
                    {
                        q = candidate;
                        q_low = candidate + 1;
                    }
                    else
                    {
                        q_high = candidate - 1;
                    }
                }

                if (q == 0)
                {
                    continue;
                }

                vector<DataT> product = ClassicMultiplication::Multiply(shifted_divisor, q, base);
                remainder = Subtract(remainder, product, base);

                if (quotient.size() <= static_cast<size_t>(d))
                {
                    quotient.resize(d + 1, 0);
                }
                quotient[d] = q;
            }

            TrimZeros(quotient);
            if (quotient.empty())
            {
                quotient.push_back(0);
            }

            TrimZeros(remainder);

            return {quotient, remainder};
        }
    };
}

#endif
