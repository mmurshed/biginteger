/**
 * BigMath: Multi-precision subtraction implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/Subtraction.h"

#include <algorithm>

namespace BigMath
{
  // Subtracts scalar b from a starting at limb aStart, borrowing upward.
  // Precondition: the value of a from aStart up is >= b (same contract as
  // before). For small bases (e.g. Base2_32) a 64-bit b can span multiple
  // limbs — the previous version assumed a single borrow and corrupted the
  // result for b >= base; this decomposes b limb by limb.
  void SubtractFrom(std::vector<DataT> &a, SizeT aStart, SizeT aEnd,
                    ULong b, BaseT base)
  {
    (void)aEnd; // kept for signature compatibility; subtraction borrows upward
    if (b == 0)
      return;

    ULong128 rem = b; // amount still to subtract at the current position
    SizeT pos = aStart;
    while (rem > 0)
    {
      if (pos >= a.size())
        a.resize(pos + 1, 0);

      ULong digit;
      if (base == Base2_64)
      {
        digit = (ULong)rem;
        rem >>= 64;
      }
      else
      {
        digit = (ULong)(rem % base);
        rem /= base;
      }

      if ((ULong)a[pos] >= digit)
        a[pos] = (DataT)((ULong)a[pos] - digit);
      else
      {
        // a[pos] + base - digit; for Base2_64 the wraparound of 64-bit
        // unsigned subtraction is exactly that.
        a[pos] = (base == Base2_64)
                     ? (DataT)((ULong)a[pos] - digit)
                     : (DataT)((ULong)a[pos] + (ULong)base - digit);
        rem += 1;
      }
      ++pos;
    }

    TrimZerosToOne(a);
  }

  void SubtractFrom(std::vector<DataT> &a, ULong b, BaseT base)
  {
    SubtractFrom(a, 0, a.size() - 1, b, base);
  }

  // Paper-pencil subtract; assumes a >= b. O(n) time, O(1) extra space.
  void Subtract(std::vector<DataT> const &a, SizeT aStart, SizeT aEnd,
                std::vector<DataT> const &b, SizeT bStart, SizeT bEnd,
                std::vector<DataT> &result, SizeT rStart,
                BaseT base)
  {
    aEnd = std::min(aEnd, (SizeT)(a.size() - 1));
    bEnd = std::min(bEnd, (SizeT)(b.size() - 1));

    Int size = std::max(Len(aStart, aEnd), Len(bStart, bEnd));

    if (base == Base2_64)
    {
      // 64-bit limb subtraction with explicit borrow tracking. Signed `Long`
      // can't hold a 64-bit limb, so do unsigned arithmetic and detect borrow
      // via comparison.
      ULong borrow = 0;
      for (Int i = 0; i < size; i++)
      {
        ULong ai = 0;
        Int aPos = aStart + i;
        if (aPos <= aEnd && aPos < (Int)a.size())
          ai = a[aPos];

        ULong bi = 0;
        Int bPos = bStart + i;
        if (bPos <= bEnd && bPos < (Int)b.size())
          bi = b[bPos];

        // Compute ai - bi - borrow with two-step borrow detection.
        ULong t1 = ai - borrow;
        ULong borrow1 = (ai < borrow) ? 1 : 0;
        ULong diff = t1 - bi;
        ULong borrow2 = (t1 < bi) ? 1 : 0;
        borrow = borrow1 + borrow2;

        Int rPos = rStart + i;
        if (rPos < (Int)result.size())
          result[rPos] = (DataT)diff;
      }
      return;
    }

    Long carry = 0;
    for (Int i = 0; i < size; i++)
    {
      Long digitOps = 0;

      Int aPos = aStart + i;
      if (aPos <= aEnd && aPos < (Int)a.size())
        digitOps = a[aPos];

      digitOps -= carry;

      Int bPos = bStart + i;
      if (bPos <= bEnd && bPos < (Int)b.size())
        digitOps -= b[bPos];

      carry = 0;
      if (digitOps < 0)
      {
        digitOps += base;
        carry = 1;
      }

      Int rPos = rStart + i;
      if (rPos < (Int)result.size())
        result[rPos] = (DataT)digitOps;
    }
  }

  void SubtractFrom(std::vector<DataT> &a, std::vector<DataT> const &b, BaseT base)
  {
    Subtract(a, 0, a.size() - 1,
             b, 0, b.size() - 1,
             a, 0,
             base);
  }

  void SubtractFrom(std::vector<DataT> &a, SizeT aStart, SizeT aEnd,
                    std::vector<DataT> const &b, SizeT bStart, SizeT bEnd,
                    BaseT base)
  {
    Subtract(a, aStart, aEnd,
             b, bStart, bEnd,
             a, aStart,
             base);
  }

  std::vector<DataT> Subtract(std::vector<DataT> const &a,
                              std::vector<DataT> const &b,
                              BaseT base)
  {
    SizeT size = (SizeT)std::max(a.size(), b.size()) + 1;
    std::vector<DataT> result(size);

    Subtract(a, 0, (SizeT)a.size() - 1,
             b, 0, (SizeT)b.size() - 1,
             result, 0, base);

    TrimZeros(result);
    return result;
  }

  void Subtract(std::vector<DataT> &result,
                std::vector<DataT> const &a,
                std::vector<DataT> const &b,
                BaseT base)
  {
    SizeT need = (SizeT)std::max(a.size(), b.size()) + 1;
    if (result.size() < need)
      result.assign(need, 0);
    else
      std::fill(result.begin(), result.end(), (DataT)0);

    Subtract(a, 0, (SizeT)a.size() - 1,
             b, 0, (SizeT)b.size() - 1,
             result, 0, base);

    TrimZeros(result);
  }
}
