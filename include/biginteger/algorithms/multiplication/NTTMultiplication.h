#ifndef NTT_MULTIPLICATION
#define NTT_MULTIPLICATION

#include <vector>
#include <algorithm>
#include <bit>
#include <cstdint>

#include "../../common/Util.h"
#include "ClassicMultiplication.h"
#include "NTTCore.h"

using namespace std;

namespace BigMath
{
    class NTTMultiplication
    {
    private:
        static vector<DataT> FinalizeBase2_32(const vector<ULong> &coeffs, SizeT coeffCount)
        {
            vector<DataT> result;
            result.reserve(coeffCount / 2 + 2);

            ULong carry = 0;
            ULong low = 0;
            bool hasLow = false;

            for (SizeT i = 0; i < coeffCount; ++i)
            {
                ULong total = coeffs[i] + carry;
                ULong digit = total & 0xFFFFULL;
                carry = total >> 16;

                if (hasLow)
                {
                    result.push_back((DataT)(low | (digit << 16)));
                    hasLow = false;
                }
                else
                {
                    low = digit;
                    hasLow = true;
                }
            }

            while (carry > 0)
            {
                ULong digit = carry & 0xFFFFULL;
                carry >>= 16;

                if (hasLow)
                {
                    result.push_back((DataT)(low | (digit << 16)));
                    hasLow = false;
                }
                else
                {
                    low = digit;
                    hasLow = true;
                }
            }

            if (hasLow)
                result.push_back((DataT)low);

            TrimZeros(result);
            return result;
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
                ULong aCoeffSize = (ULong)a.size() * 2;
                ULong bCoeffSize = (ULong)b.size() * 2;
                ULong coeffCount = aCoeffSize + bCoeffSize - 1;
                DataT n = (DataT)std::bit_ceil(coeffCount);

                static thread_local vector<ULong> fa;
                static thread_local vector<ULong> fb;
                fa.assign(n, 0);
                fb.assign(n, 0);

                for (SizeT i = 0; i < a.size(); ++i)
                {
                    SizeT j = i * 2;
                    fa[j] = a[i] & 0xFFFFULL;
                    fa[j + 1] = a[i] >> 16;
                }
                for (SizeT i = 0; i < b.size(); ++i)
                {
                    SizeT j = i * 2;
                    fb[j] = b[i] & 0xFFFFULL;
                    fb[j + 1] = b[i] >> 16;
                }

                const NTTPlan &plan = NTTCore::GetPlan((Int)n);
                NTTCore::Forward(fa, plan);
                NTTCore::Forward(fb, plan);
                for (Int i = 0; i < (Int)n; i++)
                    fa[i] = ModularField::Mul(fa[i], fb[i]);
                NTTCore::Inverse(fa, plan);

                return FinalizeBase2_32(fa, (SizeT)coeffCount);
            }
            else
            {
                ULong coeffCount = (ULong)(a.size() + b.size() - 1);
                DataT n = (DataT)std::bit_ceil(coeffCount);

                static thread_local vector<ULong> fa;
                static thread_local vector<ULong> fb;
                fa.assign(n, 0);
                fb.assign(n, 0);

                for (SizeT i = 0; i < a.size(); i++)
                    fa[i] = a[i];
                for (SizeT i = 0; i < b.size(); i++)
                    fb[i] = b[i];

                const NTTPlan &plan = NTTCore::GetPlan((Int)n);
                NTTCore::Forward(fa, plan);
                NTTCore::Forward(fb, plan);
                for (Int i = 0; i < (Int)n; i++)
                    fa[i] = ModularField::Mul(fa[i], fb[i]);
                NTTCore::Inverse(fa, plan);

                vector<DataT> c;
                c.reserve((SizeT)coeffCount + 1);
                ULong carry = 0;
                for (SizeT i = 0; i < (SizeT)coeffCount; i++)
                {
                    ULong total = fa[i] + carry;
                    c.push_back((DataT)(total % base));
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
