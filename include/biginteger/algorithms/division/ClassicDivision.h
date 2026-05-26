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
#include <stdexcept>
#include <algorithm>
using namespace std;

#include "../../common/Comparator.h"
#include "../../common/Util.h"
#include "../Addition.h"
#include "../Subtraction.h"
#include "../multiplication/ClassicMultiplication.h"

namespace BigMath
{
    class ClassicDivision
    {
    public:
        // Möller-Granlund "div2by1" reciprocal for an invariant single-limb divisor.
        // Replaces compiler __udivmodti4 (128/64 → 64) with UMULH + adds + 1-2 fixups.
        // Reference: Möller, Granlund 2010, "Improved Division by Invariant Integers", Algorithm 4.
        struct GranlundMollerDivider
        {
            ULong d;     // original divisor
            ULong dn;    // d << shift (normalized, top bit set)
            ULong v;     // floor((2^128 - 1) / dn) - 2^64
            int shift;   // clz(d)

            explicit GranlundMollerDivider(ULong divisor)
                : d(divisor), shift(__builtin_clzll(divisor))
            {
                dn = divisor << shift;
                ULong128 num = ~(ULong128)0; // 2^128 - 1
                ULong128 q = num / dn;       // q ∈ [2^64, 2^65)
                v = (ULong)q;                // low 64 bits = q - 2^64
            }

            // Divide (h*2^64 + l) by d. Quotient returned; remainder via rem.
            // Precondition: h < d (so shifted h_n < dn).
            inline ULong DivMod(ULong h, ULong l, ULong &rem) const
            {
                ULong128 combined = ((ULong128)h << 64) | l;
                ULong128 shifted = combined << shift;
                ULong u1 = (ULong)(shifted >> 64);
                ULong u0 = (ULong)shifted;

                ULong128 q = (ULong128)v * u1;
                q += ((ULong128)u1 << 64) | u0;
                ULong q1 = (ULong)(q >> 64) + 1;
                ULong q0 = (ULong)q;
                ULong r = u0 - q1 * dn;
                if (r > q0) { q1--; r += dn; }
                if (r >= dn) { q1++; r -= dn; }
                rem = r >> shift;
                return q1;
            }
        };

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

            if (w.empty())
                w.push_back(0);

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

            // Base-2^32 fast path. Möller-Granlund invariant-divisor reciprocal:
            // precompute once, then each limb step is UMULH + adds (no UDIV).
            if (base == Base2_32)
            {
                GranlundMollerDivider gm(d);
                ULong r = 0;
                for (Int i = (Int)uEnd; i >= (Int)uStart; i--)
                {
                    ULong128 acc = ((ULong128)r << 32) | u[i];
                    ULong h = (ULong)(acc >> 64);
                    ULong l = (ULong)acc;
                    DataT q = (DataT)gm.DivMod(h, l, r);
                    SizeT wPos = wStart + (i - uStart);
                    if (wPos <= wEnd)
                        w[wPos] = q;
                    else
                        w.push_back(q);
                }
                TrimZeros(w);
                return;
            }

            // Base-2^64 fast path. Same GM trick: dividend = (r * 2^64 + u[i]).
            if (base == Base2_64)
            {
                GranlundMollerDivider gm(d);
                ULong r = 0;
                for (Int i = (Int)uEnd; i >= (Int)uStart; i--)
                {
                    DataT q = (DataT)gm.DivMod(r, u[i], r);
                    SizeT wPos = wStart + (i - uStart);
                    if (wPos <= wEnd)
                        w[wPos] = q;
                    else
                        w.push_back(q);
                }
                TrimZeros(w);
                return;
            }

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

        // In-place divmod by single-limb d. Returns remainder.
        // Quotient overwrites u; trailing zero limbs trimmed.
        static DataT DivModTo(vector<DataT> &u, DataT d, BaseT base)
        {
            if (u.empty())
                return 0;

            if (base == Base2_32)
            {
                GranlundMollerDivider gm(d);
                ULong r = 0;
                for (Int i = (Int)u.size() - 1; i >= 0; --i)
                {
                    ULong128 acc = ((ULong128)r << 32) | u[i];
                    ULong h = (ULong)(acc >> 64);
                    ULong l = (ULong)acc;
                    u[i] = (DataT)gm.DivMod(h, l, r);
                }
                while (u.size() > 1 && u.back() == 0)
                    u.pop_back();
                return (DataT)r;
            }

            if (base == Base2_64)
            {
                GranlundMollerDivider gm(d);
                ULong r = 0;
                for (Int i = (Int)u.size() - 1; i >= 0; --i)
                {
                    u[i] = (DataT)gm.DivMod(r, u[i], r);
                }
                while (u.size() > 1 && u.back() == 0)
                    u.pop_back();
                return (DataT)r;
            }

            ULong r = 0;
            for (Int i = (Int)u.size() - 1; i >= 0; --i)
            {
                ULong val = r * base + u[i];
                u[i] = (DataT)(val / d);
                r = val % d;
            }
            while (u.size() > 1 && u.back() == 0)
                u.pop_back();
            return (DataT)r;
        }

        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
            vector<DataT> const &u,
            DataT d,
            BaseT base)
        {
            SizeT n = (SizeT)u.size();
            // Divide (u_n−1 . . . u_1 u_0)_b by d.
            vector<DataT> w(n);
            vector<DataT> v(1, 0);

            if (base == Base2_32)
            {
                GranlundMollerDivider gm(d);
                ULong r = 0;
                for (Int i = (Int)u.size() - 1; i >= 0; i--)
                {
                    ULong128 acc = ((ULong128)r << 32) | u[i];
                    ULong h = (ULong)(acc >> 64);
                    ULong l = (ULong)acc;
                    w[i] = (DataT)gm.DivMod(h, l, r);
                }
                v[0] = (DataT)r;
            }
            else if (base == Base2_64)
            {
                GranlundMollerDivider gm(d);
                ULong r = 0;
                for (Int i = (Int)u.size() - 1; i >= 0; i--)
                {
                    w[i] = (DataT)gm.DivMod(r, u[i], r);
                }
                v[0] = (DataT)r;
            }
            else
            {
                ULong r = 0;
                for (Int i = (Int)u.size() - 1; i >= 0; i--)
                {
                    ULong val = r * base + u[i];
                    w[i] = (DataT)(val / d);
                    r = val % d;
                }
                v[0] = (DataT)r;
            }

            TrimZerosToOne(w);
            TrimZerosToOne(v);

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

            TrimZerosToOne(quotient);

            TrimZerosToOne(remainder);

            return {quotient, remainder};
        }
    };
}

#endif
