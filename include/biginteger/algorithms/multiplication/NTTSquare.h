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
#include "ClassicSquare.h"

using namespace std;

namespace BigMath
{
  class NTTSquare
  {
  private:
    friend class NTTMultiplication;

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
    static vector<DataT> Square(vector<DataT> const &a, BaseT base)
    {
      if (IsZero(a))
        return vector<DataT>{0};
      if (a.size() == 1)
      {
        ULong sq = (ULong)a[0] * (ULong)a[0];
        if (base == Base2_32)
        {
          vector<DataT> r;
          r.push_back((DataT)(sq & 0xFFFFFFFFULL));
          if (sq >> 32) r.push_back((DataT)(sq >> 32));
          return r;
        }
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
        for (Int i = 0; i < (Int)n; i++)
          fa[i] = ModularField::Mul(fa[i], fa[i]);
        NTTCore::Inverse(fa, plan);

        return FinalizeBase2_32(fa, (SizeT)coeffCount);
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
        for (Int i = 0; i < (Int)n; i++)
          fa[i] = ModularField::Mul(fa[i], fa[i]);
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
