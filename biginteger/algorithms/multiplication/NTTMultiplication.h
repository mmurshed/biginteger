#ifndef NTT_MULTIPLICATION
#define NTT_MULTIPLICATION

#include <vector>
#include <algorithm>
#include <cstdint>

#include "../../common/Util.h"
#include "ClassicMultiplication.h"

using namespace std;

namespace BigMath
{
    // High-performance Number Theoretic Transform (NTT) Framework
    // Operating over the Fermat-like prime field modulo P = 2^64 - 2^32 + 1
    struct ModularField
    {
        static constexpr uint64_t P = 0xFFFFFFFF00000001ULL;

        // Exact modular multiplication of two 64-bit limbs producing a 64-bit result modulo P
        static inline uint64_t Mul(uint64_t a, uint64_t b)
        {
            unsigned __int128 prod = (unsigned __int128)a * b;
            return (uint64_t)(prod % P);
        }

        // Modular addition modulo P
        static inline uint64_t Add(uint64_t a, uint64_t b)
        {
            uint64_t sum = a + b;
            if (sum >= P || sum < a)
            {
                sum -= P;
            }
            return sum;
        }

        // Modular subtraction modulo P
        static inline uint64_t Sub(uint64_t a, uint64_t b)
        {
            if (a >= b) return a - b;
            return a + P - b;
        }

        // Fast modular exponentiation: base^exp modulo P
        static uint64_t Power(uint64_t base, uint64_t exp)
        {
            uint64_t res = 1;
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
        static uint64_t Inv(uint64_t a)
        {
            return Power(a, P - 2);
        }
    };

    class NTTMultiplication
    {
    private:
        // In-place iterative Cooley-Tukey NTT.
        // If invert is true, computes the inverse NTT.
        static void ntt(vector<uint64_t> &a, bool invert)
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
                // Compute the primitive len-th root of unity (7 is a generator of the field Z_P^*)
                uint64_t wlen = ModularField::Power(7, (ModularField::P - 1) / len);
                if (invert)
                {
                    wlen = ModularField::Inv(wlen);
                }

                for (Int i = 0; i < n; i += len)
                {
                    uint64_t w = 1;
                    for (int j = 0; j < len / 2; j++)
                    {
                        uint64_t u = a[i + j];
                        uint64_t v = ModularField::Mul(a[i + j + len / 2], w);
                        a[i + j] = ModularField::Add(u, v);
                        a[i + j + len / 2] = ModularField::Sub(u, v);
                        w = ModularField::Mul(w, wlen);
                    }
                }
            }
            if (invert)
            {
                uint64_t n_inv = ModularField::Inv(n);
                for (Int i = 0; i < n; i++)
                    a[i] = ModularField::Mul(a[i], n_inv);
            }
        }

        // Computes convolution of two coefficient vectors A and B.
        // Returns a vector C such that C[k] = sum_{i+j=k} A[i]*B[j] exactly.
        static vector<DataT> convolution(const vector<DataT> &A, const vector<DataT> &B)
        {
            DataT n = 1;
            while (n < (Int)A.size() + (Int)B.size() - 1)
                n <<= 1;

            vector<uint64_t> fa(n, 0), fb(n, 0);
            for (SizeT i = 0; i < A.size(); i++)
                fa[i] = A[i];
            for (SizeT i = A.size(); i < (size_t)n; i++)
                fa[i] = 0;
            for (SizeT i = 0; i < B.size(); i++)
                fb[i] = B[i];
            for (SizeT i = B.size(); i < (size_t)n; i++)
                fb[i] = 0;

            ntt(fa, false);
            ntt(fb, false);
            for (Int i = 0; i < n; i++)
                fa[i] = ModularField::Mul(fa[i], fb[i]);
            ntt(fa, true);

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
                    uint64_t low = c16[i];
                    uint64_t high = (i + 1 < c16.size()) ? c16[i + 1] : 0;
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
