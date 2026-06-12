/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef KNUTH_DIVISION
#define KNUTH_DIVISION

#include <cstddef>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>

#include "../../common/Util.h"
#include "../multiplication/ClassicMultiplication.h"
#include "../division/ClassicDivision.h"

namespace BigMath
{
  class KnuthDivision
  {
  public:
    // Divide a by b.
    // a and b are multi-precision numbers in little-endian order.
    // B is the base of each digit.
    // Returns a pair: (quotient, remainder)
    // See: D.E.Knuth 4.3.1
    // Runtime O(n^2), Space O(n)
    static std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainder(
        const std::vector<DataT> &a,
        const std::vector<DataT> &b,
        BaseT base,
        bool computeRemainder = true)
    {
      // Given nonnegative integers u = (u_m+n−1 . . . u_1 u_0)b and v = (v_n−1 . . . v_1 v_0)_b,
      // where v_n−1 != 0 and n > 1, we form the radix-b quotient ⌊u/v⌋ = (q_m q_m–1 . . . q_0)_b
      // and the remainder u mod v = (r_n−1 . . . r_1 r_0)b.

      // Check divisor nonzero
      if (IsZero(b))
        throw std::runtime_error("Division by zero.");

      // If |a| < |b| then quotient = 0 and remainder = a.
      if (Compare(a, b) < 0)
        return {std::vector<DataT>{0}, a};

      // Work on logically trimmed lengths: callers may pass vectors with
      // high zero limbs, and the algorithm needs b's true top digit for
      // normalization and a's true length for the digit count.
      SizeT n = (SizeT)b.size();
      while (n > 1 && b[n - 1] == 0)
        --n;
      SizeT aLen = (SizeT)a.size();
      while (aLen > 1 && a[aLen - 1] == 0)
        --aLen;
      SizeT m = (SizeT)(aLen - n); // m >= 0 since |a| >= |b|

      // Normalization: choose d so that b[n-1] >= B/2.
      // (There is always such a digit multiplier d in [1, B).)
      // Base2_64 is the sentinel 0, so compute B/(b_top+1) in 128 bits.
      DataT bTop = b[n - 1];
      DataT d;
      if (base == Base2_64)
        d = (bTop == ~(DataT)0)
                ? 1
                : (DataT)(LimbBase / ((ULong128)bTop + 1));
      else
        d = base / (bTop + 1);

      // Normalize u and v. Multiply trims trailing zero limbs, so restore
      // the fixed widths Knuth's loop indexes against: u as m+n+1 digits
      // (the extra top digit receives subtraction borrows), v as n digits.
      std::vector<DataT> u = ClassicMultiplication::Multiply(a, d, base);
      std::vector<DataT> v = ClassicMultiplication::Multiply(b, d, base);
      u.resize((size_t)m + n + 1, 0);
      v.resize(n, 0);

      // Initialize quotient.
      std::vector<DataT> q(m + 1, 0);

      // All per-digit arithmetic runs in 128 bits so both Base2_32 and
      // Base2_64 (sentinel 0, i.e. B = 2^64) limbs are handled exactly.
      const ULong128 B = (base == Base2_64) ? LimbBase : (ULong128)base;

      // Main loop: compute each quotient digit starting at index m downto 0.
      for (Int j = m; j >= 0; j--)
      {
        // D3: Compute estimated quotient digit qhat.
        // We form a 2-digit value from u[j+n] and u[j+n-1].
        ULong128 numerator = (ULong128)u[j + n] * B + u[j + n - 1];
        ULong128 qhat = numerator / v[n - 1];
        ULong128 rhat = numerator % v[n - 1];

        // IMPORTANT FIX: Correct qhat while it is too large.
        // In Knuth's algorithm, we must ensure that:
        //      qhat * v[n-2] <= (rhat * B + u[j+n-2])
        // (If n == 1, there is no v[n-2] so the loop does not run.)
        // qhat can initially be as large as B+1 (when u[j+n] == v[n-1]), so
        // keep decrementing while qhat >= B; only early-exit on rhat >= B
        // once qhat is a representable digit.
        while (qhat >= B || (n > 1 && qhat * v[n - 2] > (rhat * B + u[j + n - 2])))
        {
          qhat -= 1;
          rhat += v[n - 1];
          if (qhat < B && rhat >= B)
            break;
        }

        // D4: Multiply v by qhat and subtract from u[j...j+n].
        // If the subtraction borrows (i.e. the result goes negative),
        // then qhat was one too high; add back v and decrement qhat.

        // Propagate the borrow to the next digit.
        std::pair<bool, DataT> sub = subtract(u, v, (DataT)qhat, j, B);
        if (sub.first) // borrow produced
        {
          // Correction step: add back v.
          add(u, v, sub.second, j, B);
          qhat -= 1;
        }
        q[j] = (DataT)qhat;
      }

      TrimZerosToOne(q);

      std::vector<DataT> r;

      // Unnormalize the remainder.
      if (computeRemainder)
      {
        r = ClassicDivision::Divide(u, d, base);
        TrimZerosToOne(r);
      }

      return {q, r};
    }

    static std::vector<DataT> Divide(
        const std::vector<DataT> &a,
        const std::vector<DataT> &b,
        BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }

  private:
    // Subtracts qhat * v from u, starting at index j. B is the full base
    // value (2^64 for the Base2_64 sentinel), so all math is unsigned 128-bit
    // — the previous signed-Long version overflowed for qhat*v[i] >= 2^63.
    // Returns {true, borrow} if the subtraction went negative.
    static std::pair<bool, DataT> subtract(std::vector<DataT> &u, std::vector<DataT> const &v, DataT qhat, SizeT j, ULong128 B)
    {
      ULong128 borrow = 0;
      SizeT n = (SizeT)v.size();
      for (SizeT i = 0; i < n; i++)
      {
        ULong128 p = (ULong128)qhat * v[i] + borrow; // <= (B-1)^2 + (B-1) < 2^128
        DataT pd = (DataT)(p % B);
        borrow = p / B;
        DataT ud = u[i + j];
        if (ud < pd)
        {
          u[i + j] = (DataT)((ULong128)ud + B - pd);
          borrow += 1;
        }
        else
          u[i + j] = ud - pd;
      }

      // borrow < B at this point, so it fits a DataT.
      bool negative = (ULong128)u[j + n] < borrow;
      if (!negative)
        u[j + n] -= (DataT)borrow;

      return {negative, (DataT)borrow};
    }

    // Adds v to u starting at index j (used to undo an oversubtraction).
    static void add(std::vector<DataT> &u, const std::vector<DataT> &v, DataT borrow, SizeT j, ULong128 B)
    {
      DataT carry = 0;
      SizeT n = (SizeT)v.size();
      for (SizeT i = 0; i < n; i++)
      {
        ULong128 sum = (ULong128)u[i + j] + v[i] + carry;
        u[i + j] = (DataT)(sum % B);
        carry = (DataT)(sum / B);
      }
      // qhat was exactly one too high, so carry - borrow nets out mod B;
      // wrapping DataT arithmetic is intentional here.
      u[j + n] = u[j + n] + carry - borrow;
    }
  };
}

#endif
