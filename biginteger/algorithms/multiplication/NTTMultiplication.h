#ifndef NTT_MULTIPLICATION
#define NTT_MULTIPLICATION

#include <vector>
#include <algorithm>
#include <bit>
#include <cstdint>
#include <unordered_map>

#include "../../common/Util.h"
#include "ClassicMultiplication.h"

using namespace std;

namespace BigMath
{
    // High-performance Number Theoretic Transform (NTT) Framework
    // Operating over the Fermat-like prime field modulo P = 2^64 - 2^32 + 1
    struct ModularField
    {
        static constexpr ULong P = 0xFFFFFFFF00000001ULL;
        // EPSILON = 2^64 - P = 2^32 - 1. Used in fast reduction:
        //   2^64 ≡  (2^32 - 1)  (mod P)
        //   2^96 ≡ -1           (mod P)
        static constexpr ULong EPSILON = 0xFFFFFFFFULL;

        // Fast reduction of a 128-bit value modulo P, exploiting
        //   x = x_lo + 2^64·x_hi_lo + 2^96·x_hi_hi
        //     ≡ x_lo + (2^32 - 1)·x_hi_lo - x_hi_hi   (mod P)
        // Eliminates the hardware divide hidden in `prod % P`.
        static inline ULong Reduce(ULong128 x)
        {
            ULong x_lo = (ULong)x;
            ULong x_hi = (ULong)(x >> 64);
            UInt x_hi_lo = (UInt)x_hi;
            UInt x_hi_hi = (UInt)(x_hi >> 32);

            // t0 = x_lo - x_hi_hi (mod P)
            ULong t0;
            if (__builtin_sub_overflow(x_lo, (ULong)x_hi_hi, &t0))
                t0 -= EPSILON;

            // delta = x_hi_lo * (2^32 - 1)
            ULong delta = ((ULong)x_hi_lo << 32) - (ULong)x_hi_lo;

            ULong result;
            if (__builtin_add_overflow(t0, delta, &result))
                result += EPSILON;

            if (result >= P)
                result -= P;
            return result;
        }

        // Exact modular multiplication of two 64-bit limbs producing a 64-bit result modulo P
        static inline ULong Mul(ULong a, ULong b)
        {
            return Reduce((ULong128)a * b);
        }

        // Modular addition modulo P
        static inline ULong Add(ULong a, ULong b)
        {
            ULong sum = a + b;
            if (sum >= P || sum < a)
            {
                sum -= P;
            }
            return sum;
        }

        // Modular subtraction modulo P
        static inline ULong Sub(ULong a, ULong b)
        {
            if (a >= b) return a - b;
            return a + P - b;
        }

        // Fast modular exponentiation: base^exp modulo P
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

        // Fermat's Little Theorem modular inverse: a^(P-2) modulo P
        static ULong Inv(ULong a)
        {
            return Power(a, P - 2);
        }
    };

    class NTTMultiplication
    {
    private:
        // Build twiddle table roots[k] = (primitive n-th root of unity)^k for k in [0, n/2).
        // For a Cooley-Tukey butterfly at level `len`, the inner twiddle index j maps
        // to roots[j * (n / len)] — fixed stride access, no per-butterfly multiplication.
        static vector<ULong> BuildRoots(Int n, bool invert)
        {
            ULong root = ModularField::Power(7, (ModularField::P - 1) / n);
            if (invert)
                root = ModularField::Inv(root);

            vector<ULong> roots(n / 2);
            roots[0] = 1;
            for (Int k = 1; k < n / 2; ++k)
                roots[k] = ModularField::Mul(roots[k - 1], root);
            return roots;
        }

        // Thread-local cache of twiddle tables keyed by NTT size. Subsequent transforms
        // of the same size skip the O(n) BuildRoots cost — meaningful when callers
        // multiply repeatedly at the same precision (Newton iterations, Karatsuba leaves).
        static const vector<ULong> &GetRoots(Int n, bool invert)
        {
            static thread_local unordered_map<Int, vector<ULong>> fwdCache;
            static thread_local unordered_map<Int, vector<ULong>> invCache;
            auto &cache = invert ? invCache : fwdCache;
            auto it = cache.find(n);
            if (it != cache.end())
                return it->second;
            return cache.emplace(n, BuildRoots(n, invert)).first->second;
        }

        // In-place iterative Cooley-Tukey NTT using a precomputed twiddle table.
        // Removes the serial `w = Mul(w, wlen)` dependency chain from the inner butterfly.
        static void ntt(vector<ULong> &a, const vector<ULong> &roots, bool invert)
        {
            Int n = (Int)a.size();
            for (Int i = 1, j = 0; i < n; ++i)
            {
                ULong bit = n >> 1;
                while (j & bit)
                {
                    j ^= bit;
                    bit >>= 1;
                }
                j ^= bit;
                if (i < j)
                    swap(a[i], a[j]);
            }

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
            if (invert)
            {
                ULong n_inv = ModularField::Inv(n);
                for (Int i = 0; i < n; i++)
                    a[i] = ModularField::Mul(a[i], n_inv);
            }
        }

        // Computes convolution of two coefficient vectors A and B.
        // Returns a vector C such that C[k] = sum_{i+j=k} A[i]*B[j] exactly.
        static vector<DataT> convolution(const vector<DataT> &A, const vector<DataT> &B)
        {
            ULong need = (ULong)(A.size() + B.size() - 1);
            DataT n = (DataT)std::bit_ceil(need);

            vector<ULong> fa(n, 0), fb(n, 0);
            for (SizeT i = 0; i < A.size(); i++)
                fa[i] = A[i];
            for (SizeT i = 0; i < B.size(); i++)
                fb[i] = B[i];

            // Cached twiddle tables — second multiply at same size pays only a hashmap lookup.
            const vector<ULong> &fwdRoots = GetRoots(n, false);
            const vector<ULong> &invRoots = GetRoots(n, true);

            ntt(fa, fwdRoots, false);
            ntt(fb, fwdRoots, false);
            for (Int i = 0; i < n; i++)
                fa[i] = ModularField::Mul(fa[i], fb[i]);
            ntt(fa, invRoots, true);

            return fa;
        }

    public:
        // Multiply two vectors of digits using NTT-based convolution.
        static vector<DataT> Multiply(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            if (IsZero(a) || IsZero(b)) // 0 times
                return vector<DataT>();

            // If b is a single digit, use the scalar multiplication
            if (b.size() == 1)
                return ClassicMultiplication::Multiply(a, b[0], base);
            // If a is a single digit, use the scalar multiplication
            if (a.size() == 1)
                return ClassicMultiplication::Multiply(b, a[0], base);

            // If the base is Base2_32, we split each 32-bit limb into two 16-bit halves
            // to avoid overflow in the 64-bit prime field.
            if (base == Base2_32)
            {
                vector<DataT> a16, b16;
                a16.reserve(a.size() * 2);
                b16.reserve(b.size() * 2);
                
                for (SizeT i = 0; i < a.size(); ++i)
                {
                    a16.push_back(a[i] & 0xFFFFULL);
                    a16.push_back(a[i] >> 16);
                }
                for (SizeT i = 0; i < b.size(); ++i)
                {
                    b16.push_back(b[i] & 0xFFFFULL);
                    b16.push_back(b[i] >> 16);
                }
                
                // Convolve in base 2^16
                vector<DataT> c16 = convolution(a16, b16);
                
                // Carry propagation in base 2^16
                ULong carry = 0;
                for (SizeT i = 0; i < c16.size(); i++)
                {
                    ULong total = c16[i] + carry;
                    c16[i] = (DataT)(total & 0xFFFFULL);
                    carry = total >> 16;
                }
                while (carry > 0)
                {
                    c16.push_back((DataT)(carry & 0xFFFFULL));
                    carry >>= 16;
                }
                
                // Reassemble 16-bit halves back to 32-bit limbs
                vector<DataT> c32;
                c32.reserve(c16.size() / 2 + 1);
                for (SizeT i = 0; i < c16.size(); i += 2)
                {
                    ULong low = c16[i];
                    ULong high = (i + 1 < c16.size()) ? c16[i + 1] : 0;
                    c32.push_back((DataT)(low | (high << 16)));
                }
                
                TrimZeros(c32);
                return c32;
            }
            else
            {
                // Fallback for smaller bases (like Base100) where coefficients can't overflow P
                vector<DataT> c = convolution(a, b);
                ULong carry = 0;
                for (SizeT i = 0; i < c.size(); i++)
                {
                    ULong total = c[i] + carry;
                    c[i] = (DataT)(total % base);
                    carry = total / base;
                }
                while (carry > 0)
                {
                    c.push_back((DataT)(carry % base));
                    carry /= base;
                }
                TrimZeros(c);
                return c;
            }
        }
    };
}

#endif // NTT_MULTIPLICATION
