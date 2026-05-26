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
        // Base check hoisted out of loop; per-position branches replaced by phase split:
        //   phase 1: both a and b contribute
        //   phase 2: only the longer one contributes
        //   phase 3: just propagate carry to remaining lenR slots.
        static void AddPtr(
            const DataT* a, SizeT lenA,
            const DataT* b, SizeT lenB,
            DataT* r, SizeT lenR,
            BaseT base)
        {
            SizeT minLen = std::min(lenA, lenB);
            SizeT maxLen = std::max(lenA, lenB);
            const DataT* longer = (lenA >= lenB) ? a : b;

            if (base == Base2_32)
            {
                ULong carry = 0;
                for (SizeT i = 0; i < minLen; ++i)
                {
                    ULong sum = (ULong)a[i] + b[i] + carry;
                    r[i] = (DataT)(sum & 0xFFFFFFFFULL);
                    carry = sum >> 32;
                }
                for (SizeT i = minLen; i < maxLen; ++i)
                {
                    ULong sum = (ULong)longer[i] + carry;
                    r[i] = (DataT)(sum & 0xFFFFFFFFULL);
                    carry = sum >> 32;
                }
                for (SizeT i = maxLen; i < lenR; ++i)
                {
                    r[i] = (DataT)carry;
                    carry = 0;
                }
            }
            else if (base == Base2_64)
            {
                // ULong128 sum to capture the carry bit on full 64-bit limbs.
                ULong128 carry = 0;
                for (SizeT i = 0; i < minLen; ++i)
                {
                    ULong128 sum = (ULong128)a[i] + b[i] + carry;
                    r[i] = (DataT)(sum & 0xFFFFFFFFFFFFFFFFULL);
                    carry = sum >> 64;
                }
                for (SizeT i = minLen; i < maxLen; ++i)
                {
                    ULong128 sum = (ULong128)longer[i] + carry;
                    r[i] = (DataT)(sum & 0xFFFFFFFFFFFFFFFFULL);
                    carry = sum >> 64;
                }
                for (SizeT i = maxLen; i < lenR; ++i)
                {
                    r[i] = (DataT)(carry & 0xFFFFFFFFFFFFFFFFULL);
                    carry = 0;
                }
            }
            else
            {
                ULong carry = 0;
                for (SizeT i = 0; i < minLen; ++i)
                {
                    ULong sum = (ULong)a[i] + b[i] + carry;
                    r[i] = (DataT)(sum % base);
                    carry = sum / base;
                }
                for (SizeT i = minLen; i < maxLen; ++i)
                {
                    ULong sum = (ULong)longer[i] + carry;
                    r[i] = (DataT)(sum % base);
                    carry = sum / base;
                }
                for (SizeT i = maxLen; i < lenR; ++i)
                {
                    r[i] = (DataT)(carry % base);
                    carry /= base;
                }
            }
        }

        // dest[0..destLen-1] += src[0..srcLen-1]
        static void AddToPtr(
            DataT* dest, SizeT destLen,
            const DataT* src, SizeT srcLen,
            BaseT base)
        {
            SizeT shared = std::min(srcLen, destLen);

            if (base == Base2_32)
            {
                ULong carry = 0;
                SizeT i = 0;
                for (; i < shared; ++i)
                {
                    ULong sum = (ULong)dest[i] + src[i] + carry;
                    dest[i] = (DataT)(sum & 0xFFFFFFFFULL);
                    carry = sum >> 32;
                }
                for (; i < destLen && carry; ++i)
                {
                    ULong sum = (ULong)dest[i] + carry;
                    dest[i] = (DataT)(sum & 0xFFFFFFFFULL);
                    carry = sum >> 32;
                }
            }
            else if (base == Base2_64)
            {
                ULong128 carry = 0;
                SizeT i = 0;
                for (; i < shared; ++i)
                {
                    ULong128 sum = (ULong128)dest[i] + src[i] + carry;
                    dest[i] = (DataT)(sum & 0xFFFFFFFFFFFFFFFFULL);
                    carry = sum >> 64;
                }
                for (; i < destLen && carry; ++i)
                {
                    ULong128 sum = (ULong128)dest[i] + carry;
                    dest[i] = (DataT)(sum & 0xFFFFFFFFFFFFFFFFULL);
                    carry = sum >> 64;
                }
            }
            else
            {
                ULong carry = 0;
                SizeT i = 0;
                for (; i < shared; ++i)
                {
                    ULong sum = (ULong)dest[i] + src[i] + carry;
                    dest[i] = (DataT)(sum % base);
                    carry = sum / base;
                }
                for (; i < destLen && carry; ++i)
                {
                    ULong sum = (ULong)dest[i] + carry;
                    dest[i] = (DataT)(sum % base);
                    carry = sum / base;
                }
            }
        }

        // dest[0..destLen-1] -= src[0..srcLen-1]   (assumes dest >= src)
        static void SubtractFromPtr(
            DataT* dest, SizeT destLen,
            const DataT* src, SizeT srcLen,
            BaseT base)
        {
            SizeT shared = std::min(srcLen, destLen);

            if (base == Base2_32)
            {
                Long borrow = 0;
                SizeT i = 0;
                for (; i < shared; ++i)
                {
                    Long diff = (Long)dest[i] - (Long)src[i] - borrow;
                    if (diff < 0) { diff += (Long)1 << 32; borrow = 1; }
                    else borrow = 0;
                    dest[i] = (DataT)diff;
                }
                for (; i < destLen && borrow; ++i)
                {
                    Long diff = (Long)dest[i] - borrow;
                    if (diff < 0) { diff += (Long)1 << 32; borrow = 1; }
                    else borrow = 0;
                    dest[i] = (DataT)diff;
                }
            }
            else if (base == Base2_64)
            {
                // 64-bit limbs don't fit in signed Long; use unsigned with
                // two-step borrow detection (matches Subtraction.cpp pattern).
                ULong borrow = 0;
                SizeT i = 0;
                for (; i < shared; ++i)
                {
                    ULong ai = dest[i];
                    ULong bi = src[i];
                    ULong t1 = ai - borrow;
                    ULong borrow1 = (ai < borrow) ? 1 : 0;
                    ULong diff = t1 - bi;
                    ULong borrow2 = (t1 < bi) ? 1 : 0;
                    borrow = borrow1 + borrow2;
                    dest[i] = (DataT)diff;
                }
                for (; i < destLen && borrow; ++i)
                {
                    ULong ai = dest[i];
                    ULong diff = ai - borrow;
                    borrow = (ai < borrow) ? 1 : 0;
                    dest[i] = (DataT)diff;
                }
            }
            else
            {
                Long borrow = 0;
                SizeT i = 0;
                for (; i < shared; ++i)
                {
                    Long diff = (Long)dest[i] - (Long)src[i] - borrow;
                    if (diff < 0) { diff += base; borrow = 1; }
                    else borrow = 0;
                    dest[i] = (DataT)diff;
                }
                for (; i < destLen && borrow; ++i)
                {
                    Long diff = (Long)dest[i] - borrow;
                    if (diff < 0) { diff += base; borrow = 1; }
                    else borrow = 0;
                    dest[i] = (DataT)diff;
                }
            }
        }

        // Schoolbook multiplication for base cases: r = a * b.
        // Base2_32 path packs pairs of 32-bit limbs into 64-bit values and runs the
        // schoolbook in 64-bit limb space — 4× fewer multiplies, each is UMULL+UMULH
        // on ARM64 (≈ 2× cost per multiply). Net ≈ 2× over scalar 32-bit schoolbook.
        // Buffers up to 64 packed limbs (covers Karatsuba leaf threshold of 48) live
        // on the stack to avoid heap allocation; larger inputs fall back to heap.
        static void MultiplyClassicPtr(
            const DataT* a, SizeT lenA,
            const DataT* b, SizeT lenB,
            DataT* r,
            BaseT base)
        {
            std::memset(r, 0, (lenA + lenB) * sizeof(DataT));

            if (base == Base2_32)
            {
                SizeT na64 = (lenA + 1) >> 1;
                SizeT nb64 = (lenB + 1) >> 1;
                SizeT nr64 = na64 + nb64;

                constexpr SizeT STACK_CAP = 64;
                ULong a64_stack[STACK_CAP];
                ULong b64_stack[STACK_CAP];
                ULong r64_stack[2 * STACK_CAP];
                std::unique_ptr<ULong[]> heap;
                ULong* a64 = a64_stack;
                ULong* b64 = b64_stack;
                ULong* r64 = r64_stack;
                if (na64 > STACK_CAP || nb64 > STACK_CAP)
                {
                    heap.reset(new ULong[na64 + nb64 + nr64]);
                    a64 = heap.get();
                    b64 = a64 + na64;
                    r64 = b64 + nb64;
                }

                // Pack a, b: a64[k] = a[2k] | (a[2k+1] << 32)
                for (SizeT k = 0; k < na64; ++k)
                {
                    SizeT lo = 2 * k;
                    ULong v = a[lo];
                    if (lo + 1 < lenA) v |= (ULong)a[lo + 1] << 32;
                    a64[k] = v;
                }
                for (SizeT k = 0; k < nb64; ++k)
                {
                    SizeT lo = 2 * k;
                    ULong v = b[lo];
                    if (lo + 1 < lenB) v |= (ULong)b[lo + 1] << 32;
                    b64[k] = v;
                }

                std::memset(r64, 0, nr64 * sizeof(ULong));

                // 64-bit schoolbook: prod is 128-bit, carry is 64-bit.
                for (SizeT i = 0; i < nb64; ++i)
                {
                    ULong bi = b64[i];
                    if (bi == 0) continue;
                    ULong carry = 0;
                    for (SizeT j = 0; j < na64; ++j)
                    {
                        ULong128 prod = (ULong128)a64[j] * bi + r64[i + j] + carry;
                        r64[i + j] = (ULong)prod;
                        carry = (ULong)(prod >> 64);
                    }
                    r64[i + na64] = carry;
                }

                // Unpack r64 → r. r has length lenA+lenB; r64 may carry one extra
                // limb whose high 32 bits are always zero (since |product| ≤ B^(lenA+lenB)).
                SizeT rLen = lenA + lenB;
                for (SizeT k = 0; k < nr64; ++k)
                {
                    SizeT lo = 2 * k;
                    if (lo < rLen) r[lo] = (DataT)(r64[k] & 0xFFFFFFFFULL);
                    if (lo + 1 < rLen) r[lo + 1] = (DataT)(r64[k] >> 32);
                }
            }
            else if (base == Base2_64)
            {
                // Plain 64×64→128 schoolbook. No pack-to-larger trick available
                // (would need 256-bit primitives), so just one multiply per cell.
                for (SizeT i = 0; i < lenB; ++i)
                {
                    if (b[i] == 0) continue;
                    ULong carry = 0;
                    for (SizeT j = 0; j < lenA; ++j)
                    {
                        ULong128 prod = (ULong128)a[j] * b[i] + r[i + j] + carry;
                        r[i + j] = (DataT)(prod & 0xFFFFFFFFFFFFFFFFULL);
                        carry = (ULong)(prod >> 64);
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
            //
            // No memset(c) here — recursive calls fully overwrite c[0..la+lb-1]:
            //   call 1 writes c[0..2m-1] (= T1)
            //   call 2 writes c[2m..la+lb-1] (= T2)
            // At leaf, MultiplyClassicPtr memsets its own output region.

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
            // No memset(t3) — recursive call fully overwrites t3[0..lenT3-1].

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

            // Workspace usage per recursion level: lenWl + lenWh (~n+2) for the
            // sum operands, then lenT1 + lenT2 + lenT3 (~3n) for the three sub-
            // products plus T3 = (wl*wh) which is sized 2*(n/2+1). Sub-recursion
            // is dispatched with `nextW` pointing past these.
            //
            // The geometric sum across log₂(n) levels works out to ~10n in the
            // worst case (highly unbalanced split — one half ~all of n, the
            // other ~1 — produces wl/wh sizes near n, pushing the per-level
            // usage above the balanced-case 5n).
            //
            // 16n is a comfortable upper bound that keeps the workspace within
            // a single allocation and never triggers a heap-buffer-overflow
            // even on adversarially-skewed inputs (e.g. 1200×1024).
            unique_ptr<DataT[]> w(new DataT[16 * n]);

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
