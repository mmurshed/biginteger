/**
 * BigInteger Class
 * Version 10.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef KARATSUBA_MULTIPLICATION
#define KARATSUBA_MULTIPLICATION

#include <vector>
#include <algorithm>
#include <cstring>
#include <memory>
using namespace std;

#include "../../common/Util.h"
#include "../multiplication/ClassicMultiplication.h"

namespace BigMath
{
    class KaratsubaMultiplication
    {
    private:
#ifndef BIGMATH_KARATSUBA_THRESHOLD
#define BIGMATH_KARATSUBA_THRESHOLD 48
#endif
        static const SizeT KARATSUBA_THRESHOLD = BIGMATH_KARATSUBA_THRESHOLD;

        // Adds a[0..lenA-1] and b[0..lenB-1] and writes to r[0..lenR-1].
        static void AddPtr(
            const DataT* a, SizeT lenA,
            const DataT* b, SizeT lenB,
            DataT* r, SizeT lenR,
            BaseT base)
        {
            ULong carry = 0;
            for (SizeT i = 0; i < lenR; ++i)
            {
                ULong sum = carry;
                if (i < lenA) sum += a[i];
                if (i < lenB) sum += b[i];
                
                if (base == Base2_32)
                {
                    r[i] = (DataT)(sum & 0xFFFFFFFFULL);
                    carry = sum >> 32;
                }
                else
                {
                    r[i] = (DataT)(sum % base);
                    carry = sum / base;
                }
            }
        }

        // dest[0..destLen-1] += src[0..srcLen-1]
        static void AddToPtr(
            DataT* dest, SizeT destLen,
            const DataT* src, SizeT srcLen,
            BaseT base)
        {
            ULong carry = 0;
            for (SizeT i = 0; i < destLen; ++i)
            {
                ULong sum = dest[i] + carry;
                if (i < srcLen) sum += src[i];
                
                if (base == Base2_32)
                {
                    dest[i] = (DataT)(sum & 0xFFFFFFFFULL);
                    carry = sum >> 32;
                }
                else
                {
                    dest[i] = (DataT)(sum % base);
                    carry = sum / base;
                }
                if (carry == 0 && i >= srcLen) break;
            }
        }

        // dest[0..destLen-1] -= src[0..srcLen-1]
        // Assumes dest >= src (non-negative result)
        static void SubtractFromPtr(
            DataT* dest, SizeT destLen,
            const DataT* src, SizeT srcLen,
            BaseT base)
        {
            Long borrow = 0;
            for (SizeT i = 0; i < destLen; ++i)
            {
                Long diff = (Long)dest[i] - borrow;
                if (i < srcLen) diff -= (Long)src[i];
                
                if (diff < 0)
                {
                    diff += base;
                    borrow = 1;
                }
                else
                {
                    borrow = 0;
                }
                dest[i] = (DataT)diff;
                if (borrow == 0 && i >= srcLen) break;
            }
        }

        // Schoolbook multiplication for base cases: r = a * b.
        static void MultiplyClassicPtr(
            const DataT* a, SizeT lenA,
            const DataT* b, SizeT lenB,
            DataT* r,
            BaseT base)
        {
            std::memset(r, 0, (lenA + lenB) * sizeof(DataT));

            if (base == Base2_32)
            {
                for (SizeT i = 0; i < lenB; ++i)
                {
                    if (b[i] == 0) continue;
                    ULong carry = 0;
                    for (SizeT j = 0; j < lenA; ++j)
                    {
                        ULong prod = (ULong)a[j] * b[i] + r[i + j] + carry;
                        r[i + j] = (DataT)(prod & 0xFFFFFFFFULL);
                        carry = prod >> 32;
                    }
                    r[i + lenA] = (DataT)carry;
                }
            }
            else
            {
                for (SizeT i = 0; i < lenB; ++i)
                {
                    if (b[i] == 0) continue;
                    ULong carry = 0;
                    for (SizeT j = 0; j < lenA; ++j)
                    {
                        ULong prod = (ULong)a[j] * b[i] + r[i + j] + carry;
                        r[i + j] = (DataT)(prod % base);
                        carry = prod / base;
                    }
                    r[i + lenA] = (DataT)carry;
                }
            }
        }

        static void MultiplyRecursive(
            const DataT* a, SizeT lenA,
            const DataT* b, SizeT lenB,
            DataT* c, // Output buffer (size >= lenA + lenB)
            DataT* w, // Workspace buffer
            BaseT base)
        {
            SizeT la = lenA;
            SizeT lb = lenB;

            if (la <= KARATSUBA_THRESHOLD || lb <= KARATSUBA_THRESHOLD)
            {
                MultiplyClassicPtr(a, la, b, lb, c, base);
                return;
            }

            SizeT m = (max(la, lb) + 1) / 2;
            if (m >= la) m = la - 1;
            if (m >= lb) m = lb - 1;

            SizeT lenAl = m;
            SizeT lenAh = la - m;
            SizeT lenBl = m;
            SizeT lenBh = lb - m;

            SizeT lenWl = max(lenAl, lenAh) + 1;
            SizeT lenWh = max(lenBl, lenBh) + 1;
            
            DataT* wl = w;
            DataT* wh = w + lenWl;
            DataT* nextW = wh + lenWh;

            AddPtr(a, lenAl, a + m, lenAh, wl, lenWl, base);
            AddPtr(b, lenBl, b + m, lenBh, wh, lenWh, base);

            // c = al*bl + B^m * ((al+ah)*(bl+bh) - al*bl - ah*bh) + B^2m * ah*bh
            std::memset(c, 0, (la + lb) * sizeof(DataT));

            // 1. T1 = al * bl -> c[0..2m-1]
            MultiplyRecursive(a, lenAl, b, lenBl, c, nextW, base);

            // 2. T2 = ah * bh -> c[2m..la+lb-1]
            MultiplyRecursive(a + m, lenAh, b + m, lenBh, c + 2 * m, nextW, base);

            // 3. T3 = wl * wh -> computed in workspace after wl/wh
            SizeT lenT1 = 2 * m;
            SizeT lenT2 = lenAh + lenBh;
            SizeT lenT3 = lenWl + lenWh;

            // Safe offset: place t3 at w + lenT1 + lenT2 to avoid any overlap with t1Copy/t2Copy
            DataT* t3 = w + lenT1 + lenT2;
            DataT* nextW2 = t3 + lenT3;
            std::memset(t3, 0, lenT3 * sizeof(DataT));

            MultiplyRecursive(wl, lenWl, wh, lenWh, t3, nextW2, base);

            // Now safely copy T1 and T2 to t1Copy and t2Copy (which are at w and w + lenT1)
            DataT* t1Copy = w;
            DataT* t2Copy = w + lenT1;
            std::memcpy(t1Copy, c, lenT1 * sizeof(DataT));
            std::memcpy(t2Copy, c + 2 * m, lenT2 * sizeof(DataT));

            // Perform subtractions on t3 directly in the workspace.
            // Since T3 >= T1 + T2 is mathematically guaranteed, t3 remains non-negative.
            SubtractFromPtr(t3, lenT3, t1Copy, lenT1, base);
            SubtractFromPtr(t3, lenT3, t2Copy, lenT2, base);

            // c[m..la+lb-1] += (T3 - T1 - T2)
            AddToPtr(c + m, la + lb - m, t3, lenT3, base);
        }

    public:
        static vector<DataT> Multiply(
            vector<DataT> const &a,
            vector<DataT> const &b,
            BaseT base)
        {
            if (IsZero(a) || IsZero(b))
                return vector<DataT>();

            if (b.size() == 1)
                return ClassicMultiplication::Multiply(a, b[0], base);
            if (a.size() == 1)
                return ClassicMultiplication::Multiply(b, a[0], base);

            SizeT n = (SizeT)max(a.size(), b.size());
            vector<DataT> c(a.size() + b.size(), 0);
            
            // Workspace is overwritten before use at every recursive level; avoid
            // zero-initializing up to 8*n limbs on every Karatsuba call.
            unique_ptr<DataT[]> w(new DataT[8 * n]);

            MultiplyRecursive(
                a.data(), a.size(),
                b.data(), b.size(),
                c.data(),
                w.get(),
                base);

            TrimZeros(c);
            return c;
        }
    };
}

#endif
