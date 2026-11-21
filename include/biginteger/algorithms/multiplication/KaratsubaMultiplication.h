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

#include "../../common/LimbKernels.h"
#include "../../common/Util.h"
#include "../multiplication/ClassicMultiplication.h"

namespace BigMath
{
    class KaratsubaMultiplication
    {
    private:
// Leaf cutoff: operands at or below this go to the schoolbook kernel.
// 64-bit limbs: re-swept 32 -> 24 on 2026-11-20 with the flag-chain limb
// kernels (common/LimbKernels.h) — the add/sub passes got ~2x cheaper, so one
// more recursion level pays (3-7% at 32/64/128/256/512 limbs; 24 and 28 tie).
#ifndef BIGMATH_KARATSUBA_THRESHOLD
#if BIGMATH_LIMB_64
#define BIGMATH_KARATSUBA_THRESHOLD 24
#else
#define BIGMATH_KARATSUBA_THRESHOLD 32
#endif
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
                // Flag-chain kernels — see common/LimbKernels.h.
                LimbKernel::AddExt(r, lenR, a, lenA, b, lenB);
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
                LimbKernel::AddInPlace(dest, destLen, src, srcLen);
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
                LimbKernel::SubInPlace(dest, destLen, src, srcLen);
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
        // Buffers up to 64 packed limbs (covers the Karatsuba leaf threshold of 32) live
        // on the stack to avoid heap allocation; larger inputs fall back to heap.
        static void MultiplyClassicPtr(
            const DataT* a, SizeT lenA,
            const DataT* b, SizeT lenB,
            DataT* r,
            BaseT base)
        {
            if (base == Base2_64 && lenA != 0 && lenB != 0)
            {
                // mul_1 for the first row, addmul_1 for the rest: writes every
                // limb of r, so no zero-fill. See common/LimbKernels.h.
                LimbKernel::MulBasecase(r, a, lenA, b, lenB);
                return;
            }

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
                // Only reached with an empty operand; r is already zeroed.
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

            // Unbalanced operands: slice the longer one into pieces the length
            // of the shorter and multiply piece by piece.
            //
            // The split below clamps m to min(la, lb) - 1. When one operand is
            // at least twice the other, that clamp makes every level peel just
            // lb - 1 limbs off the long operand: the (al+ah)*(bl+bh) product is
            // (la - lb + 2) x lb, barely smaller than the input, so the
            // recursion runs la / lb levels deep instead of log2. Each level
            // parks ~2*la + lb limbs of workspace, the total grows like
            // la^2 / lb, and from about la > 16*lb it ran past the 16n
            // workspace allocated in Multiply() — a heap overflow on shapes the
            // dispatcher really sends here (e.g. 1214 x 65 limbs). Slicing
            // keeps every sub-product balanced, bounds the workspace by ~4*la,
            // and removes the O(la^2 / lb) add/copy overhead as well.
            if (la >= 2 * lb || lb >= 2 * la)
            {
                const DataT* x = a;  // longer operand
                SizeT lx = la;
                const DataT* y = b;  // shorter operand
                SizeT ly = lb;
                if (lx < ly)
                {
                    std::swap(x, y);
                    std::swap(lx, ly);
                }

                // First piece straight into c[0 .. 2*ly); the rest of c starts at zero.
                MultiplyRecursive(x, ly, y, ly, c, w, base);
                std::memset(c + 2 * ly, 0, (lx - ly) * sizeof(DataT));

                // Remaining pieces: product into the workspace, then add at its offset.
                DataT* t = w;
                for (SizeT off = ly; off < lx; off += ly)
                {
                    SizeT len = std::min(ly, lx - off);
                    MultiplyRecursive(x + off, len, y, ly, t, t + len + ly, base);
                    AddToPtr(c + off, lx + ly - off, t, len + ly, base);
                }
                return;
            }

            SizeT m = (std::max(la, lb) + 1) / 2;
            if (m >= la) m = la - 1;
            if (m >= lb) m = lb - 1;

            SizeT lenAl = m;
            SizeT lenAh = la - m;
            SizeT lenBl = m;
            SizeT lenBh = lb - m;

            SizeT lenWl = std::max(lenAl, lenAh) + 1;
            SizeT lenWh = std::max(lenBl, lenBh) + 1;
            
            DataT* wl = w;
            DataT* wh = w + lenWl;
            DataT* nextW = wh + lenWh;

            AddPtr(a, lenAl, a + m, lenAh, wl, lenWl, base);
            AddPtr(b, lenBl, b + m, lenBh, wh, lenWh, base);

            // c = al*bl + B^m * ((al+ah)*(bl+bh) - al*bl - ah*bh) + B^2m * ah*bh
            //
            // No std::memset(c) here — recursive calls fully overwrite c[0..la+lb-1]:
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

            DataT* t3 = nextW; // after both live sums
            DataT* nextW2 = t3 + lenT3;
            // No std::memset(t3) — recursive call fully overwrites t3[0..lenT3-1].

            MultiplyRecursive(wl, lenWl, wh, lenWh, t3, nextW2, base);

            // Both output slices remain intact until the final shifted add.
            SubtractFromPtr(t3, lenT3, c, lenT1, base);
            SubtractFromPtr(t3, lenT3, c + 2 * m, lenT2, base);

            // c[m..la+lb-1] += (T3 - T1 - T2)
            AddToPtr(c + m, la + lb - m, t3, lenT3, base);
        }

    public:
        // The sliced recursion needs scratch proportional to its shorter input.
        static size_t ScratchLimbs(SizeT la, SizeT lb)
        {
            const size_t shorter = std::min(la, lb), longer = std::max(la, lb);
            if (shorter <= KARATSUBA_THRESHOLD) return 0;
            return longer / 2 >= shorter ? 10 * shorter + 1024 : 16 * longer;
        }

        // Non-overlapping output (la + lb limbs) and ScratchLimbs(la,lb) scratch.
        static void MultiplyInto(const DataT *a, SizeT la, const DataT *b, SizeT lb,
                                 DataT *out, DataT *scratch, BaseT base)
        {
            MultiplyRecursive(a, la, b, lb, out, scratch, base);
        }

        static std::vector<DataT> Multiply(
            std::vector<DataT> const &a,
            std::vector<DataT> const &b,
            BaseT base)
        {
            if (IsZero(a) || IsZero(b))
                return std::vector<DataT>();

            if (b.size() == 1)
                return ClassicMultiplication::Multiply(a, b[0], base);
            if (a.size() == 1)
                return ClassicMultiplication::Multiply(b, a[0], base);

            std::vector<DataT> c(a.size() + b.size());
            const size_t count = ScratchLimbs(a.size(), b.size());
            std::unique_ptr<DataT[]> scratch(count ? new DataT[count] : nullptr);
            MultiplyInto(a.data(), a.size(), b.data(), b.size(), c.data(), scratch.get(), base);
            TrimZeros(c);
            return c;
        }
    };
}

#endif
