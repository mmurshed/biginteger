#ifndef BIGMATH_NTT_CORE
#define BIGMATH_NTT_CORE

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "../../common/Util.h"
#include "../../common/Parallel.h"

namespace BigMath
{
    // High-performance Number Theoretic Transform (NTT) framework over
    // P = 2^64 - 2^32 + 1.
    struct ModularField
    {
        static constexpr ULong P = 0xFFFFFFFF00000001ULL;
        static constexpr ULong EPSILON = 0xFFFFFFFFULL;

        static inline ULong Reduce(ULong128 x)
        {
            ULong x_lo = (ULong)x;
            ULong x_hi = (ULong)(x >> 64);
            UInt x_hi_lo = (UInt)x_hi;
            UInt x_hi_hi = (UInt)(x_hi >> 32);

            ULong t0;
            ULong borrow = (ULong)__builtin_sub_overflow(x_lo, (ULong)x_hi_hi, &t0);
            t0 -= EPSILON & -borrow;

            ULong delta = ((ULong)x_hi_lo << 32) - (ULong)x_hi_lo;

            ULong result;
            ULong carry = (ULong)__builtin_add_overflow(t0, delta, &result);
            result += EPSILON & -carry;

            ULong over = (ULong)(result >= P);
            result -= P & -over;
            return result;
        }

        static inline ULong Mul(ULong a, ULong b)
        {
            return Reduce((ULong128)a * b);
        }

        static inline ULong Add(ULong a, ULong b)
        {
            ULong sum;
            ULong overflow = (ULong)__builtin_add_overflow(a, b, &sum);
            ULong over_p = (ULong)(sum >= P);
            ULong mask = -(overflow | over_p);
            sum -= P & mask;
            return sum;
        }

        static inline ULong Sub(ULong a, ULong b)
        {
            ULong diff;
            ULong borrow = (ULong)__builtin_sub_overflow(a, b, &diff);
            diff += P & -borrow;
            return diff;
        }

        static ULong Power(ULong base, ULong exp)
        {
            ULong res = 1;
            base %= P;
            while (exp > 0)
            {
                if (exp & 1) res = Mul(res, base);
                base = Mul(base, base);
                exp >>= 1;
            }
            return res;
        }

        static ULong Inv(ULong a)
        {
            return Power(a, P - 2);
        }
    };

    struct NTTPlan
    {
        std::vector<ULong> forwardRoots;
        std::vector<ULong> inverseRoots;
        ULong inverseSize = 1;
    };

    class NTTCore
    {
    private:
        static std::vector<ULong> BuildRoots(Int n, bool invert)
        {
            ULong root = ModularField::Power(7, (ModularField::P - 1) / n);
            if (invert)
                root = ModularField::Inv(root);

            Int half = n / 2;
            std::vector<ULong> roots(half);
            ULong *rPtr = roots.data();
            // Parallel chunk fill mirroring NttCrt::BuildRoots — each thread
            // jumps ahead via Power then sequential Mul within its range.
            if ((SizeT)half >= ParallelMinSize())
            {
                ParallelFor(half, [rPtr, root](Int start, Int end) {
                    if (end <= start) return;
                    ULong cur = (start == 0) ? (ULong)1 : ModularField::Power(root, (ULong)start);
                    rPtr[start] = cur;
                    for (Int k = start + 1; k < end; ++k)
                    {
                        cur = ModularField::Mul(cur, root);
                        rPtr[k] = cur;
                    }
                });
            }
            else
            {
                rPtr[0] = 1;
                for (Int k = 1; k < half; ++k)
                    rPtr[k] = ModularField::Mul(rPtr[k - 1], root);
            }
            return roots;
        }

        static NTTPlan BuildPlan(Int n)
        {
            NTTPlan plan;
            plan.forwardRoots = BuildRoots(n, false);
            plan.inverseRoots = BuildRoots(n, true);
            plan.inverseSize = ModularField::Inv((ULong)n);
            return plan;
        }

    public:
        static const NTTPlan &GetPlan(Int n)
        {
            static thread_local std::unordered_map<Int, NTTPlan> cache;
            auto it = cache.find(n);
            if (it != cache.end())
                return it->second;
            return cache.emplace(n, BuildPlan(n)).first->second;
        }

        // Forward DIF and inverse DIT pair. Forward leaves data bit-reversed;
        // pointwise multiplication stays aligned; inverse DIT restores natural order.
        //
        // Radix-4 fused butterflies handle two adjacent radix-2 layers in one
        // memory pass: load 4 elements, do 4 muls + 8 add/subs, store 4 — vs
        // 2 loads + 2 stores per element across two separate radix-2 layers.
        // Halves bandwidth for paired layers; same modular op count. Straggler
        // radix-2 layer runs when log2(n) is odd (Forward at len=2, Inverse at
        // the bottom before the radix-4 pairs start at outer_len=4 or 8).
        //
        // Parallelism: per-layer, ParallelFor over i-blocks (one block = `len`
        // elements). Top layer (1 i-block) falls back to serial.
    private:
        static inline void ForwardRadix2Layer(ULong *a, Int n, Int len, const ULong *roots)
        {
            Int halflen = len >> 1;
            Int stride = n / len;
            Int numBlocks = n / len;
            auto body = [a, halflen, len, stride, roots](Int bStart, Int bEnd) {
                for (Int b = bStart; b < bEnd; ++b)
                {
                    Int i = b * len;
                    for (Int j = 0; j < halflen; ++j)
                    {
                        ULong u = a[i + j];
                        ULong v = a[i + j + halflen];
                        a[i + j] = ModularField::Add(u, v);
                        a[i + j + halflen] = ModularField::Mul(ModularField::Sub(u, v), roots[j * stride]);
                    }
                }
            };
            if ((SizeT)(n / 2) >= ParallelMinSize() && numBlocks > 1)
                ParallelFor(numBlocks, body);
            else
                body(0, numBlocks);
        }

        static inline void ForwardRadix4Layer(ULong *a, Int n, Int outerLen, const ULong *roots)
        {
            Int halflen = outerLen >> 1;
            Int qlen = outerLen >> 2;
            Int stride = n / outerLen;
            Int numBlocks = n / outerLen;
            Int omega4_off = n / 4;
            auto body = [a, halflen, qlen, outerLen, stride, omega4_off, roots](Int bStart, Int bEnd) {
                for (Int b = bStart; b < bEnd; ++b)
                {
                    Int i = b * outerLen;
                    for (Int j = 0; j < qlen; ++j)
                    {
                        ULong x0 = a[i + j];
                        ULong x1 = a[i + j + qlen];
                        ULong x2 = a[i + j + halflen];
                        ULong x3 = a[i + j + halflen + qlen];

                        ULong g  = roots[j * stride];
                        ULong gw = roots[j * stride + omega4_off];
                        ULong g2 = roots[2 * j * stride];

                        ULong t0 = ModularField::Add(x0, x2);
                        ULong t1 = ModularField::Add(x1, x3);
                        ULong t2 = ModularField::Mul(ModularField::Sub(x0, x2), g);
                        ULong t3 = ModularField::Mul(ModularField::Sub(x1, x3), gw);

                        a[i + j]                  = ModularField::Add(t0, t1);
                        a[i + j + qlen]           = ModularField::Mul(ModularField::Sub(t0, t1), g2);
                        a[i + j + halflen]        = ModularField::Add(t2, t3);
                        a[i + j + halflen + qlen] = ModularField::Mul(ModularField::Sub(t2, t3), g2);
                    }
                }
            };
            if ((SizeT)(n / 2) >= ParallelMinSize() && numBlocks > 1)
                ParallelFor(numBlocks, body);
            else
                body(0, numBlocks);
        }

        static inline void InverseRadix2Layer(ULong *a, Int n, Int len, const ULong *roots)
        {
            Int halflen = len >> 1;
            Int stride = n / len;
            Int numBlocks = n / len;
            auto body = [a, halflen, len, stride, roots](Int bStart, Int bEnd) {
                for (Int b = bStart; b < bEnd; ++b)
                {
                    Int i = b * len;
                    for (Int j = 0; j < halflen; ++j)
                    {
                        ULong u = a[i + j];
                        ULong v = ModularField::Mul(a[i + j + halflen], roots[j * stride]);
                        a[i + j] = ModularField::Add(u, v);
                        a[i + j + halflen] = ModularField::Sub(u, v);
                    }
                }
            };
            if ((SizeT)(n / 2) >= ParallelMinSize() && numBlocks > 1)
                ParallelFor(numBlocks, body);
            else
                body(0, numBlocks);
        }

        static inline void InverseRadix4Layer(ULong *a, Int n, Int outerLen, const ULong *roots)
        {
            Int halflen = outerLen >> 1;
            Int qlen = outerLen >> 2;
            Int stride = n / outerLen;
            Int numBlocks = n / outerLen;
            Int omega4_off = n / 4;
            auto body = [a, halflen, qlen, outerLen, stride, omega4_off, roots](Int bStart, Int bEnd) {
                for (Int b = bStart; b < bEnd; ++b)
                {
                    Int i = b * outerLen;
                    for (Int j = 0; j < qlen; ++j)
                    {
                        ULong x0 = a[i + j];
                        ULong x1 = a[i + j + qlen];
                        ULong x2 = a[i + j + halflen];
                        ULong x3 = a[i + j + halflen + qlen];

                        ULong g  = roots[j * stride];
                        ULong gw = roots[j * stride + omega4_off];
                        ULong g2 = roots[2 * j * stride];

                        ULong y1 = ModularField::Mul(x1, g2);
                        ULong y3 = ModularField::Mul(x3, g2);
                        ULong t0 = ModularField::Add(x0, y1);
                        ULong t1 = ModularField::Sub(x0, y1);
                        ULong t2 = ModularField::Mul(ModularField::Add(x2, y3), g);
                        ULong t3 = ModularField::Mul(ModularField::Sub(x2, y3), gw);

                        a[i + j]                  = ModularField::Add(t0, t2);
                        a[i + j + halflen]        = ModularField::Sub(t0, t2);
                        a[i + j + qlen]           = ModularField::Add(t1, t3);
                        a[i + j + halflen + qlen] = ModularField::Sub(t1, t3);
                    }
                }
            };
            if ((SizeT)(n / 2) >= ParallelMinSize() && numBlocks > 1)
                ParallelFor(numBlocks, body);
            else
                body(0, numBlocks);
        }

        // Radix-8 fused DIF butterfly. Combines three radix-2 layers (len, len/2,
        // len/4) into a single load/store across 8 elements. Internally a pair of
        // radix-4 fused butterflies (slots {0,2,4,6} and {1,3,5,7}) followed by
        // four radix-2 butterflies on the resulting pairs.
        static inline void ForwardRadix8Layer(ULong *a, Int n, Int outerLen, const ULong *roots)
        {
            Int halflen = outerLen >> 1;
            Int qlen4 = outerLen >> 2;
            Int qlen8 = outerLen >> 3;
            Int stride = n / outerLen;
            Int stride4 = stride << 2;
            Int numBlocks = n / outerLen;
            Int omega4_off = n / 4;
            auto body = [a, halflen, qlen4, qlen8, outerLen, stride, stride4, omega4_off, roots](Int bStart, Int bEnd) {
                for (Int b = bStart; b < bEnd; ++b)
                {
                    Int i = b * outerLen;
                    for (Int j = 0; j < qlen8; ++j)
                    {
                        ULong x0 = a[i + j];
                        ULong x1 = a[i + j + qlen8];
                        ULong x2 = a[i + j + qlen4];
                        ULong x3 = a[i + j + qlen4 + qlen8];
                        ULong x4 = a[i + j + halflen];
                        ULong x5 = a[i + j + halflen + qlen8];
                        ULong x6 = a[i + j + halflen + qlen4];
                        ULong x7 = a[i + j + halflen + qlen4 + qlen8];

                        ULong g_a  = roots[j * stride];
                        ULong gw_a = roots[j * stride + omega4_off];
                        ULong g2_a = roots[2 * j * stride];
                        ULong g_b  = roots[(j + qlen8) * stride];
                        ULong gw_b = roots[(j + qlen8) * stride + omega4_off];
                        ULong g2_b = roots[2 * (j + qlen8) * stride];
                        ULong g3   = roots[j * stride4];

                        ULong t0 = ModularField::Add(x0, x4);
                        ULong t1 = ModularField::Add(x2, x6);
                        ULong t2 = ModularField::Mul(ModularField::Sub(x0, x4), g_a);
                        ULong t3 = ModularField::Mul(ModularField::Sub(x2, x6), gw_a);
                        ULong y0 = ModularField::Add(t0, t1);
                        ULong y2 = ModularField::Mul(ModularField::Sub(t0, t1), g2_a);
                        ULong y4 = ModularField::Add(t2, t3);
                        ULong y6 = ModularField::Mul(ModularField::Sub(t2, t3), g2_a);

                        ULong u0 = ModularField::Add(x1, x5);
                        ULong u1 = ModularField::Add(x3, x7);
                        ULong u2 = ModularField::Mul(ModularField::Sub(x1, x5), g_b);
                        ULong u3 = ModularField::Mul(ModularField::Sub(x3, x7), gw_b);
                        ULong y1 = ModularField::Add(u0, u1);
                        ULong y3 = ModularField::Mul(ModularField::Sub(u0, u1), g2_b);
                        ULong y5 = ModularField::Add(u2, u3);
                        ULong y7 = ModularField::Mul(ModularField::Sub(u2, u3), g2_b);

                        ULong z0 = ModularField::Add(y0, y1);
                        ULong z1 = ModularField::Mul(ModularField::Sub(y0, y1), g3);
                        ULong z2 = ModularField::Add(y2, y3);
                        ULong z3 = ModularField::Mul(ModularField::Sub(y2, y3), g3);
                        ULong z4 = ModularField::Add(y4, y5);
                        ULong z5 = ModularField::Mul(ModularField::Sub(y4, y5), g3);
                        ULong z6 = ModularField::Add(y6, y7);
                        ULong z7 = ModularField::Mul(ModularField::Sub(y6, y7), g3);

                        a[i + j]                          = z0;
                        a[i + j + qlen8]                  = z1;
                        a[i + j + qlen4]                  = z2;
                        a[i + j + qlen4 + qlen8]          = z3;
                        a[i + j + halflen]                = z4;
                        a[i + j + halflen + qlen8]        = z5;
                        a[i + j + halflen + qlen4]        = z6;
                        a[i + j + halflen + qlen4 + qlen8] = z7;
                    }
                }
            };
            if ((SizeT)(n / 2) >= ParallelMinSize() && numBlocks > 1)
                ParallelFor(numBlocks, body);
            else
                body(0, numBlocks);
        }

        // Radix-8 fused DIT butterfly. Inverse of ForwardRadix8Layer: applies
        // three DIT radix-2 layers (len/4 → len/2 → len) over 8 elements in a
        // single load/store.
        static inline void InverseRadix8Layer(ULong *a, Int n, Int outerLen, const ULong *roots)
        {
            Int halflen = outerLen >> 1;
            Int qlen4 = outerLen >> 2;
            Int qlen8 = outerLen >> 3;
            Int stride = n / outerLen;
            Int stride4 = stride << 2;
            Int numBlocks = n / outerLen;
            Int omega4_off = n / 4;
            auto body = [a, halflen, qlen4, qlen8, outerLen, stride, stride4, omega4_off, roots](Int bStart, Int bEnd) {
                for (Int b = bStart; b < bEnd; ++b)
                {
                    Int i = b * outerLen;
                    for (Int j = 0; j < qlen8; ++j)
                    {
                        ULong x0 = a[i + j];
                        ULong x1 = a[i + j + qlen8];
                        ULong x2 = a[i + j + qlen4];
                        ULong x3 = a[i + j + qlen4 + qlen8];
                        ULong x4 = a[i + j + halflen];
                        ULong x5 = a[i + j + halflen + qlen8];
                        ULong x6 = a[i + j + halflen + qlen4];
                        ULong x7 = a[i + j + halflen + qlen4 + qlen8];

                        // Layer A (innermost, len=L/4, twiddle = roots[j*stride4]).
                        ULong g_a = roots[j * stride4];
                        ULong v01 = ModularField::Mul(x1, g_a);
                        ULong y0 = ModularField::Add(x0, v01);
                        ULong y1 = ModularField::Sub(x0, v01);
                        ULong v23 = ModularField::Mul(x3, g_a);
                        ULong y2 = ModularField::Add(x2, v23);
                        ULong y3 = ModularField::Sub(x2, v23);
                        ULong v45 = ModularField::Mul(x5, g_a);
                        ULong y4 = ModularField::Add(x4, v45);
                        ULong y5 = ModularField::Sub(x4, v45);
                        ULong v67 = ModularField::Mul(x7, g_a);
                        ULong y6 = ModularField::Add(x6, v67);
                        ULong y7 = ModularField::Sub(x6, v67);

                        // Layer B (len=L/2). Pairs (0,2),(1,3),(4,6),(5,7).
                        ULong g_b1 = roots[2 * j * stride];
                        ULong g_b2 = roots[2 * (j + qlen8) * stride];
                        ULong w02 = ModularField::Mul(y2, g_b1);
                        ULong z0 = ModularField::Add(y0, w02);
                        ULong z2 = ModularField::Sub(y0, w02);
                        ULong w13 = ModularField::Mul(y3, g_b2);
                        ULong z1 = ModularField::Add(y1, w13);
                        ULong z3 = ModularField::Sub(y1, w13);
                        ULong w46 = ModularField::Mul(y6, g_b1);
                        ULong z4 = ModularField::Add(y4, w46);
                        ULong z6 = ModularField::Sub(y4, w46);
                        ULong w57 = ModularField::Mul(y7, g_b2);
                        ULong z5 = ModularField::Add(y5, w57);
                        ULong z7 = ModularField::Sub(y5, w57);

                        // Layer C (outermost, len=L). Pairs (0,4),(1,5),(2,6),(3,7).
                        ULong g_c0 = roots[j * stride];
                        ULong g_c1 = roots[(j + qlen8) * stride];
                        ULong g_c2 = roots[(j + qlen4) * stride];
                        ULong g_c3 = roots[(j + qlen4 + qlen8) * stride];
                        ULong w04 = ModularField::Mul(z4, g_c0);
                        ULong w15 = ModularField::Mul(z5, g_c1);
                        ULong w26 = ModularField::Mul(z6, g_c2);
                        ULong w37 = ModularField::Mul(z7, g_c3);

                        a[i + j]                          = ModularField::Add(z0, w04);
                        a[i + j + halflen]                = ModularField::Sub(z0, w04);
                        a[i + j + qlen8]                  = ModularField::Add(z1, w15);
                        a[i + j + halflen + qlen8]        = ModularField::Sub(z1, w15);
                        a[i + j + qlen4]                  = ModularField::Add(z2, w26);
                        a[i + j + halflen + qlen4]        = ModularField::Sub(z2, w26);
                        a[i + j + qlen4 + qlen8]          = ModularField::Add(z3, w37);
                        a[i + j + halflen + qlen4 + qlen8] = ModularField::Sub(z3, w37);
                    }
                }
            };
            if ((SizeT)(n / 2) >= ParallelMinSize() && numBlocks > 1)
                ParallelFor(numBlocks, body);
            else
                body(0, numBlocks);
        }

    public:
        static void Forward(std::vector<ULong> &a, const NTTPlan &plan)
        {
            Int n = (Int)a.size();
            if (n <= 1) return;
            ULong *aPtr = a.data();
            const ULong *roots = plan.forwardRoots.data();

            // Triple from the top: (n, n/2, n/4), (n/8, n/16, n/32), …
            // Straggler runs at len ∈ {4, 2} depending on log2(n) mod 3.
            Int len = n;
            while (len >= 8)
            {
                ForwardRadix8Layer(aPtr, n, len, roots);
                len >>= 3;
            }
            if (len == 4)
                ForwardRadix4Layer(aPtr, n, 4, roots);
            else if (len == 2)
                ForwardRadix2Layer(aPtr, n, 2, roots);
        }

        static void Inverse(std::vector<ULong> &a, const NTTPlan &plan)
        {
            Int n = (Int)a.size();
            ULong *aPtr = a.data();
            const ULong *roots = plan.inverseRoots.data();

            if (n >= 2)
            {
                // Triple from the bottom: (2,4,8), (16,32,64), …
                // If log2(n) % 3 ∈ {1,2}, a radix-2 or radix-4 head runs first
                // so subsequent triples start at outer_len ∈ {8, 16, 32}.
                Int logn = __builtin_ctzll((unsigned long long)n);
                Int rem = logn % 3;
                Int len;
                if (rem == 1)
                {
                    InverseRadix2Layer(aPtr, n, 2, roots);
                    len = 16;
                }
                else if (rem == 2)
                {
                    InverseRadix4Layer(aPtr, n, 4, roots);
                    len = 32;
                }
                else
                {
                    len = 8;
                }
                while (len <= n)
                {
                    InverseRadix8Layer(aPtr, n, len, roots);
                    len <<= 3;
                }
            }

            // Final scaling by inverse size — trivially parallel.
            ULong invSize = plan.inverseSize;
            auto scaleBody = [aPtr, invSize](Int s, Int e) {
                for (Int i = s; i < e; ++i)
                    aPtr[i] = ModularField::Mul(aPtr[i], invSize);
            };
            if ((SizeT)n >= ParallelMinSize())
                ParallelFor(n, scaleBody);
            else
                scaleBody(0, n);
        }
    };
}

#endif
