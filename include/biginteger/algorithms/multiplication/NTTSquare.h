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
#include "ClassicSquare.h"

using namespace std;

namespace BigMath
{
  class NTTSquare
  {
  private:
    friend class NTTMultiplication;

    // Internals duplicated from NTTMultiplication since they're private there.
    // Keep in sync: ModularField + nttForward/nttInverse are stable.

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

    static void nttForward(vector<ULong> &a, const vector<ULong> &roots)
    {
      Int n = (Int)a.size();
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

    static void nttInverse(vector<ULong> &a, const vector<ULong> &roots)
    {
      Int n = (Int)a.size();
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
      ULong n_inv = ModularField::Inv(n);
      for (Int i = 0; i < n; i++)
        a[i] = ModularField::Mul(a[i], n_inv);
    }

    // Self-convolution: A * A. One FFT instead of two.
    static vector<DataT> convolutionSelfBase2_32(vector<DataT> const &A)
    {
      ULong aCoeffSize = (ULong)A.size() * 2;
      ULong need = 2 * aCoeffSize - 1;
      DataT n = (DataT)std::bit_ceil(need);

      vector<ULong> fa(n, 0);
      for (SizeT i = 0; i < A.size(); ++i)
      {
        SizeT j = i * 2;
        fa[j] = A[i] & 0xFFFFULL;
        fa[j + 1] = A[i] >> 16;
      }

      const vector<ULong> &fwdRoots = GetRoots(n, false);
      const vector<ULong> &invRoots = GetRoots(n, true);

      nttForward(fa, fwdRoots);
      for (Int i = 0; i < n; i++)
        fa[i] = ModularField::Mul(fa[i], fa[i]);
      nttInverse(fa, invRoots);

      return fa;
    }

    static vector<DataT> convolutionSelf(vector<DataT> const &A)
    {
      ULong need = (ULong)(2 * A.size() - 1);
      DataT n = (DataT)std::bit_ceil(need);

      vector<ULong> fa(n, 0);
      for (SizeT i = 0; i < A.size(); i++)
        fa[i] = A[i];

      const vector<ULong> &fwdRoots = GetRoots(n, false);
      const vector<ULong> &invRoots = GetRoots(n, true);

      nttForward(fa, fwdRoots);
      for (Int i = 0; i < n; i++)
        fa[i] = ModularField::Mul(fa[i], fa[i]);
      nttInverse(fa, invRoots);

      return fa;
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
        vector<DataT> c16 = convolutionSelfBase2_32(a);

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
        vector<DataT> c = convolutionSelf(a);
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

#endif
