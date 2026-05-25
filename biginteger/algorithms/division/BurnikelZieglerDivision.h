#ifndef BURNIKELZIEGLER_DIVISION
#define BURNIKELZIEGLER_DIVISION

#include <cstring>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;

#include "../../common/Comparator.h"
#include "../../common/Util.h"
#include "../Addition.h"
#include "../Multiplication.h"
#include "../Shift.h"
#include "FastDivision.h"

namespace BigMath
{
  class BurnikelZieglerDivision
  {
  private:
    static const SizeT BZ_THRESHOLD = 512;

    // out = (high << shift_limbs) + low. One allocation.
    static vector<DataT> AddShifted(
        vector<DataT> const &high, SizeT shift,
        const DataT *lowData, SizeT lowSize,
        BaseT base)
    {
      SizeT outSize = max(high.size() + shift, (size_t)lowSize) + 1;
      vector<DataT> out(outSize, 0);

      std::memcpy(out.data(), lowData, lowSize * sizeof(DataT));

      ULong carry = 0;
      SizeT i = 0;
      for (; i < high.size(); ++i)
      {
        ULong sum = (ULong)out[shift + i] + high[i] + carry;
        if (base == Base2_32)
        {
          out[shift + i] = (DataT)(sum & 0xFFFFFFFFULL);
          carry = sum >> 32;
        }
        else
        {
          out[shift + i] = (DataT)(sum % base);
          carry = sum / base;
        }
      }
      SizeT idx = shift + i;
      while (carry)
      {
        if (idx >= out.size())
          out.push_back(0);
        ULong sum = (ULong)out[idx] + carry;
        if (base == Base2_32)
        {
          out[idx] = (DataT)(sum & 0xFFFFFFFFULL);
          carry = sum >> 32;
        }
        else
        {
          out[idx] = (DataT)(sum % base);
          carry = sum / base;
        }
        ++idx;
      }
      return out;
    }

    // Span-based recursion. `a` and `b` are read-only views — the `high` slice
    // is passed as a subspan (zero allocation) instead of copied.
    static pair<vector<DataT>, vector<DataT>> DivideRecursive(
        span<const DataT> a,
        span<const DataT> b,
        BaseT base,
        bool computeRemainder)
    {
      if (a.size() <= BZ_THRESHOLD || b.size() <= BZ_THRESHOLD || a.size() <= 2 * b.size())
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      SizeT lowSize = 2 * (SizeT)b.size();
      // high = a[lowSize..end), low = a[0..lowSize). No allocation for either.
      span<const DataT> high = a.subspan(lowSize);

      auto highQr = DivideRecursive(high, b, base, true);

      // combined = (highQr.second << lowSize) + a[0..lowSize-1]
      vector<DataT> combined = AddShifted(highQr.second, lowSize, a.data(), lowSize, base);
      TrimZeros(combined);

      auto lowQr = FastDivision::DivideAndRemainder(combined, b, base, computeRemainder);

      // quotient = (highQr.first << lowSize) + lowQr.first
      vector<DataT> quotient = AddShifted(highQr.first, lowSize,
                                          lowQr.first.data(), (SizeT)lowQr.first.size(), base);
      TrimZeros(quotient);
      if (quotient.empty())
        quotient.push_back(0);

      return {quotient, computeRemainder ? lowQr.second : vector<DataT>()};
    }

  public:
    static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
        span<const DataT> a,
        span<const DataT> b,
        BaseT base,
        bool computeRemainder = true)
    {
      if (IsZero(b))
        throw invalid_argument("Division by zero");

      if (IsZero(a))
        return {vector<DataT>{0}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

      Int cmp = Compare(a, b);
      if (cmp < 0)
        return {vector<DataT>{0}, computeRemainder ? vector<DataT>(a.begin(), a.end()) : vector<DataT>()};
      if (cmp == 0)
        return {vector<DataT>{1}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

      return DivideRecursive(a, b, base, computeRemainder);
    }

    static vector<DataT> Divide(span<const DataT> a, span<const DataT> b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }

    // Vector overloads — backward compat.
    static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base,
        bool computeRemainder = true)
    {
      return DivideAndRemainder(span<const DataT>(a), span<const DataT>(b), base, computeRemainder);
    }

    static vector<DataT> Divide(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }
  };
}

#endif
