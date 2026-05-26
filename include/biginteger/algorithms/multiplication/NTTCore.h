#ifndef BIGMATH_NTT_CORE
#define BIGMATH_NTT_CORE

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "../../common/Util.h"

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

            std::vector<ULong> roots(n / 2);
            roots[0] = 1;
            for (Int k = 1; k < n / 2; ++k)
                roots[k] = ModularField::Mul(roots[k - 1], root);
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

        // Forward DIF and inverse DIT pair. The forward transform leaves data in
        // bit-reversed order; pointwise multiplication stays aligned and inverse DIT
        // returns natural-order coefficients.
        static void Forward(std::vector<ULong> &a, const NTTPlan &plan)
        {
            Int n = (Int)a.size();
            const std::vector<ULong> &roots = plan.forwardRoots;
            for (Int len = n; len > 1; len >>= 1)
            {
                Int halflen = len >> 1;
                Int stride = n / len;
                for (Int i = 0; i < n; i += len)
                {
                    for (Int j = 0; j < halflen; ++j)
                    {
                        ULong u = a[i + j];
                        ULong v = a[i + j + halflen];
                        a[i + j] = ModularField::Add(u, v);
                        a[i + j + halflen] = ModularField::Mul(ModularField::Sub(u, v), roots[j * stride]);
                    }
                }
            }
        }

        static void Inverse(std::vector<ULong> &a, const NTTPlan &plan)
        {
            Int n = (Int)a.size();
            const std::vector<ULong> &roots = plan.inverseRoots;
            for (Int len = 2; len <= n; len <<= 1)
            {
                Int halflen = len >> 1;
                Int stride = n / len;
                for (Int i = 0; i < n; i += len)
                {
                    for (Int j = 0; j < halflen; ++j)
                    {
                        ULong u = a[i + j];
                        ULong v = ModularField::Mul(a[i + j + halflen], roots[j * stride]);
                        a[i + j] = ModularField::Add(u, v);
                        a[i + j + halflen] = ModularField::Sub(u, v);
                    }
                }
            }

            for (Int i = 0; i < n; i++)
                a[i] = ModularField::Mul(a[i], plan.inverseSize);
        }
    };
}

#endif
