#ifndef NTT_MULTIPLICATION
#define NTT_MULTIPLICATION

#include <vector>
#include <algorithm>
#include <bit>
#include <cstdint>

#include "../../common/Util.h"
#include "ClassicMultiplication.h"
#include "NTTCore.h"
#include "NTTFinalize.h"
#include "NTTMultiplicationCrt.h"

using namespace std;

namespace BigMath
{
    class NTTMultiplication
    {
    private:
    public:
        using PreparedOperand = NttCrt::PreparedOperand;

        // Precompute the CRT NTT spectra for `operand` so repeated
        // multiplications by operands up to `maxOtherLimbs` can skip one side's
        // packing and forward transforms. This targets repeated same-operand
        // workloads; one-off Multiply(a,b) remains unchanged.
        static PreparedOperand PrepareOperand(
            const vector<DataT> &operand,
            SizeT maxOtherLimbs,
            BaseT base)
        {
            return NttCrt::PrepareOperand(operand, maxOtherLimbs, base);
        }

        static vector<DataT> Multiply(
            const PreparedOperand &prepared,
            const vector<DataT> &other)
        {
            return NttCrt::Multiply(prepared, other);
        }

        // (a · b) mod (B^L − 1) via a cyclic NTT of half the full-product
        // length. See NttCrt::MultiplyMod2km1 for the caller contract.
        static vector<DataT> MultiplyMod2km1(
            const vector<DataT> &a,
            const vector<DataT> &b,
            SizeT L,
            BaseT base)
        {
            return NttCrt::MultiplyMod2km1(a, b, L, base);
        }

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

#if BIGMATH_NTT_CRT
            // Size-gated hybrid. The NEON Shoup butterflies (2026-06-12)
            // made CRT faster than single-prime Goldilocks at every measured
            // size, so the gate sits at 256 limbs sum — effectively always.
            // Override via -DBIGMATH_NTT_CRT_THRESHOLD=N.
#ifndef BIGMATH_NTT_CRT_THRESHOLD
#define BIGMATH_NTT_CRT_THRESHOLD 256
#endif
            if (a.size() + b.size() >= BIGMATH_NTT_CRT_THRESHOLD)
                return NttCrt::Multiply(a, b, base);
            // Fall through to Goldilocks below threshold.
#endif

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
                {
                    ULong *faPtr = fa.data();
                    ULong *fbPtr = fb.data();
                    auto body = [faPtr, fbPtr](Int s, Int e) {
                        for (Int i = s; i < e; ++i)
                            faPtr[i] = ModularField::Mul(faPtr[i], fbPtr[i]);
                    };
                    if ((SizeT)n >= ParallelMinSize())
                        ParallelFor((Int)n, body);
                    else
                        body(0, (Int)n);
                }
                NTTCore::Inverse(fa, plan);

                return NttFinalizeBase2_32(fa, (SizeT)coeffCount);
            }
            else if (base == Base2_64)
            {
                // Split each 64-bit limb into four 16-bit halves so coefficients
                // stay within the Goldilocks field's accumulation bound. Same
                // prime as Base2_32 (P = 2^64 - 2^32 + 1); convolution overflow
                // bound is still ~2^31 limbs (≈ 16 GB operands).
                ULong aCoeffSize = (ULong)a.size() * 4;
                ULong bCoeffSize = (ULong)b.size() * 4;
                ULong coeffCount = aCoeffSize + bCoeffSize - 1;
                DataT n = (DataT)std::bit_ceil(coeffCount);

                static thread_local vector<ULong> fa;
                static thread_local vector<ULong> fb;
                fa.assign(n, 0);
                fb.assign(n, 0);

                for (SizeT i = 0; i < a.size(); ++i)
                {
                    SizeT j = i * 4;
                    fa[j]     = a[i] & 0xFFFFULL;
                    fa[j + 1] = (a[i] >> 16) & 0xFFFFULL;
                    fa[j + 2] = (a[i] >> 32) & 0xFFFFULL;
                    fa[j + 3] = (a[i] >> 48) & 0xFFFFULL;
                }
                for (SizeT i = 0; i < b.size(); ++i)
                {
                    SizeT j = i * 4;
                    fb[j]     = b[i] & 0xFFFFULL;
                    fb[j + 1] = (b[i] >> 16) & 0xFFFFULL;
                    fb[j + 2] = (b[i] >> 32) & 0xFFFFULL;
                    fb[j + 3] = (b[i] >> 48) & 0xFFFFULL;
                }

                const NTTPlan &plan = NTTCore::GetPlan((Int)n);
                NTTCore::Forward(fa, plan);
                NTTCore::Forward(fb, plan);
                {
                    ULong *faPtr = fa.data();
                    ULong *fbPtr = fb.data();
                    auto body = [faPtr, fbPtr](Int s, Int e) {
                        for (Int i = s; i < e; ++i)
                            faPtr[i] = ModularField::Mul(faPtr[i], fbPtr[i]);
                    };
                    if ((SizeT)n >= ParallelMinSize())
                        ParallelFor((Int)n, body);
                    else
                        body(0, (Int)n);
                }
                NTTCore::Inverse(fa, plan);

                return NttFinalizeBase2_64(fa, (SizeT)coeffCount);
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
                {
                    ULong *faPtr = fa.data();
                    ULong *fbPtr = fb.data();
                    auto body = [faPtr, fbPtr](Int s, Int e) {
                        for (Int i = s; i < e; ++i)
                            faPtr[i] = ModularField::Mul(faPtr[i], fbPtr[i]);
                    };
                    if ((SizeT)n >= ParallelMinSize())
                        ParallelFor((Int)n, body);
                    else
                        body(0, (Int)n);
                }
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
