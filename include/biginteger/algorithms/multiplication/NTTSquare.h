/**
 * BigInteger Class
 * NTT-based squaring: single forward DIF + pointwise self-multiply + inverse DIT.
 * Half the FFT work of NTTMultiplication.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef NTT_SQUARE
#define NTT_SQUARE

#include <vector>
#include <bit>
#include <cstdint>

#include "../../common/Util.h"
#include "NTTMultiplication.h"
#include "NTTCore.h"
#include "NTTFinalize.h"
#include "ClassicSquare.h"

using namespace std;

namespace BigMath
{
  class NTTSquare
  {
  private:
  public:
    static vector<DataT> Square(vector<DataT> const &a, BaseT base)
    {
      if (IsZero(a))
        return vector<DataT>{0};
      if (a.size() == 1)
      {
        if (base == Base2_32)
        {
          ULong sq = (ULong)a[0] * (ULong)a[0];
          vector<DataT> r;
          r.push_back((DataT)(sq & 0xFFFFFFFFULL));
          if (sq >> 32) r.push_back((DataT)(sq >> 32));
          return r;
        }
        if (base == Base2_64)
        {
          ULong128 sq = (ULong128)a[0] * a[0];
          vector<DataT> r;
          r.push_back((DataT)(sq & 0xFFFFFFFFFFFFFFFFULL));
          if ((ULong)(sq >> 64)) r.push_back((DataT)(sq >> 64));
          return r;
        }
        ULong sq = (ULong)a[0] * (ULong)a[0];
        vector<DataT> r;
        r.push_back((DataT)(sq % base));
        if (sq / base) r.push_back((DataT)(sq / base));
        return r;
      }

      if (base == Base2_32)
      {
        ULong aCoeffSize = (ULong)a.size() * 2;
        ULong coeffCount = 2 * aCoeffSize - 1;
        DataT n = (DataT)std::bit_ceil(coeffCount);

        static thread_local vector<ULong> fa;
        fa.assign(n, 0);
        for (SizeT i = 0; i < a.size(); ++i)
        {
          SizeT j = i * 2;
          fa[j] = a[i] & 0xFFFFULL;
          fa[j + 1] = a[i] >> 16;
        }

        const NTTPlan &plan = NTTCore::GetPlan((Int)n);
        NTTCore::Forward(fa, plan);
        {
          ULong *faPtr = fa.data();
          auto body = [faPtr](Int s, Int e) {
            for (Int i = s; i < e; ++i)
              faPtr[i] = ModularField::Mul(faPtr[i], faPtr[i]);
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
        // Split each 64-bit limb into four 16-bit coefficients (mirrors
        // NTTMultiplication Base2_64 path). FinalizeBase2_64 reassembles 4-into-1.
        ULong aCoeffSize = (ULong)a.size() * 4;
        ULong coeffCount = 2 * aCoeffSize - 1;
        DataT n = (DataT)std::bit_ceil(coeffCount);

        static thread_local vector<ULong> fa;
        fa.assign(n, 0);
        for (SizeT i = 0; i < a.size(); ++i)
        {
          SizeT j = i * 4;
          fa[j]     = a[i] & 0xFFFFULL;
          fa[j + 1] = (a[i] >> 16) & 0xFFFFULL;
          fa[j + 2] = (a[i] >> 32) & 0xFFFFULL;
          fa[j + 3] = (a[i] >> 48) & 0xFFFFULL;
        }

        const NTTPlan &plan = NTTCore::GetPlan((Int)n);
        NTTCore::Forward(fa, plan);
        {
          ULong *faPtr = fa.data();
          auto body = [faPtr](Int s, Int e) {
            for (Int i = s; i < e; ++i)
              faPtr[i] = ModularField::Mul(faPtr[i], faPtr[i]);
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
        ULong coeffCount = (ULong)(2 * a.size() - 1);
        DataT n = (DataT)std::bit_ceil(coeffCount);

        static thread_local vector<ULong> fa;
        fa.assign(n, 0);
        for (SizeT i = 0; i < a.size(); i++)
          fa[i] = a[i];

        const NTTPlan &plan = NTTCore::GetPlan((Int)n);
        NTTCore::Forward(fa, plan);
        {
          ULong *faPtr = fa.data();
          auto body = [faPtr](Int s, Int e) {
            for (Int i = s; i < e; ++i)
              faPtr[i] = ModularField::Mul(faPtr[i], faPtr[i]);
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

#endif
