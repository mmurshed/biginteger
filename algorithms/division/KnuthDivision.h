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
#include <algorithm>
using namespace std;

#include "../../BigInteger.h"
#include "../../BigIntegerUtil.h"
#include "../multiplication/ClassicMultiplication.h"

namespace BigMath
{
  class KnuthDivision
  {
  public:
    static void DivideTo(
        vector<DataT> &u,
        DataT d,
        BaseT base)
    {
      Divide(u, 0, u.size() - 1, d, u, 0, u.size() - 1, base);
    }

    static void DivideTo(
        vector<DataT> &u, SizeT uStart, SizeT uEnd,
        DataT d,
        BaseT base)
    {
      Divide(u, uStart, uEnd, d, u, uStart, uEnd, base);
    }

    static vector<DataT> Divide(
        vector<DataT> const &u,
        DataT d,
        BaseT base)
    {
      SizeT n = (SizeT)u.size();
      // Divide (u_n−1 . . . u_1 u_0)_b by d.
      vector<DataT> w(n);

      Divide(u, 0, u.size() - 1, d, w, 0, w.size() - 1, base);

      return w;
    }

    static void Divide(
        vector<DataT> const &u, SizeT uStart, SizeT uEnd,
        DataT d,
        vector<DataT> &w, SizeT wStart, SizeT wEnd,
        BaseT base)
    {
      // Number of limbs in [uStart..uEnd]
      SizeT n = (uEnd >= uStart) ? (uEnd - uStart + 1) : 0;
      if (n == 0)
      {
        // Nothing to divide
        w.clear();
        return;
      }

      // Make sure 'w' can hold at least n limbs for the quotient
      SizeT neededSize = wStart + n;
      if (w.size() < neededSize)
        w.resize(neededSize, 0);

      // 'r' is the running remainder.
      ULong r = 0;

      // Process from the highest index (uEnd) down to the lowest (uStart).
      //   For little-endian, u[uEnd] is the most significant limb.
      for (Int i = (Int)uEnd; i >= (Int)uStart; i--)
      {
        ULong val = r * base + u[i]; // Combine remainder + next limb
        DataT q = (DataT)(val / d);  // Quotient digit
        r = val % d;                 // Remainder

        // The quotient digit for this place goes at wPos (also in little-endian).
        SizeT wPos = wStart + (i - uStart);
        if (wPos <= wEnd)
          w[wPos] = q;
        else
          w.push_back(q);
      }

      // Remove any leading zero-limbs in the quotient.
      BigIntegerUtil::TrimZeros(w);
    }

    static pair< vector<DataT>, vector<DataT> > DivideAndRemainder(
        vector<DataT> const &u,
        DataT d,
        BaseT base)
    {
      SizeT n = (SizeT)u.size();
      // Divide (u_n−1 . . . u_1 u_0)_b by d.
      vector<DataT> w(n);
      vector<DataT> v(n);

      DivideAndRemainder(
          u, 0, u.size() - 1,
          d,
          w, 0, w.size() - 1,
          v, 0, v.size() - 1,
          base);

      return {w, v};
    }

    static void DivideAndRemainder(
        vector<DataT> const &u, SizeT uStart, SizeT uEnd,
        DataT d,
        vector<DataT> &w, SizeT wStart, SizeT wEnd,
        vector<DataT> &v, SizeT vStart, SizeT vEnd,
        BaseT base)
    {
      // Number of limbs in [uStart..uEnd]
      SizeT n = (uEnd >= uStart) ? (uEnd - uStart + 1) : 0;
      if (n == 0)
      {
        // Nothing to divide
        return;
      }

      // Make sure 'w' can hold at least n limbs for the quotient
      SizeT neededSize = wStart + n;
      if (w.size() < neededSize)
        w.resize(neededSize, 0);

      // 'r' is the running remainder.
      ULong r = 0;

      // Process from the highest index (uEnd) down to the lowest (uStart).
      //   For little-endian, u[uEnd] is the most significant limb.
      for (Int i = (Int)uEnd; i >= (Int)uStart; i--)
      {
        ULong val = r * base + u[i]; // Combine remainder + next limb
        DataT q = (DataT)(val / d);  // Quotient digit
        r = val % d;                 // Remainder

        // The quotient digit for this place goes at wPos (also in little-endian).
        SizeT wPos = wStart + (i - uStart);
        if (wPos <= wEnd)
          w[wPos] = q;
        else
          w.push_back(q);

        // The quotient digit for this place goes at wPos (also in little-endian).
        SizeT vPos = vStart + (i - uStart);
        if (vPos <= vEnd)
          v[vPos] = r;
        else
          v.push_back(r);
      }

      // Remove any leading zero-limbs in the quotient.
      BigIntegerUtil::TrimZeros(w);
      BigIntegerUtil::TrimZeros(v);
    }

    static void DivideTo(
        vector<DataT> &a,
        vector<DataT> const &b,
        BaseT base)
    {
      pair<vector<DataT>, vector<DataT>> results = DivideAndRemainder(a, b, base);
      a = results.first;
    }

    // Divide a by b.
    // a and b are multi-precision numbers in little-endian order.
    // B is the base of each digit.
    // Returns a pair: (quotient, remainder)
    // See: D.E.Knuth 4.3.1
    // Runtime O(n^2), Space O(n)
    static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
        const vector<DataT> &a,
        const vector<DataT> &b,
        BaseT base)
    {
      // Given nonnegative integers u = (u_m+n−1 . . . u_1 u_0)b and v = (v_n−1 . . . v_1 v_0)_b,
      // where v_n−1 != 0 and n > 1, we form the radix-b quotient ⌊u/v⌋ = (q_m q_m–1 . . . q_0)_b
      // and the remainder u mod v = (r_n−1 . . . r_1 r_0)b.

      // Check divisor nonzero
      if (b.empty() || (b.size() == 1 && b[0] == 0))
        throw runtime_error("Division by zero.");

      // If |a| < |b| then quotient = 0 and remainder = a.
      if (BigIntegerComparator::Compare(a, b) < 0)
        return {vector<DataT>{0}, a};

      SizeT n = (SizeT)b.size();
      SizeT m = (SizeT)(a.size() - n); // m >= 0

      // Normalization: choose d so that b[n-1] >= B/2.
      // (There is always such a digit multiplier d in [1, B).)
      DataT d = base / (b.back() + 1);

      // Normalize u and v.
      vector<DataT> u = ClassicMultiplication::Multiply(a, d, base);
      vector<DataT> v = ClassicMultiplication::Multiply(b, d, base);

      // Append an extra digit to u (u[m+n] becomes available).
      u.push_back(0);

      // Initialize quotient.
      vector<DataT> q(m + 1, 0);

      // Main loop: compute each quotient digit starting at index m downto 0.
      for (Int j = m; j >= 0; j--)
      {
        // D3: Compute estimated quotient digit qhat.
        // We form a 2-digit value from u[j+n] and u[j+n-1].
        ULong numerator = u[j + n] * base + u[j + n - 1];
        ULong qhat = numerator / v[n - 1];
        ULong rhat = numerator % v[n - 1];

        // IMPORTANT FIX: Correct qhat while it is too large.
        // In Knuth's algorithm, we must ensure that:
        //      qhat * v[n-2] <= (rhat * B + u[j+n-2])
        // (If n == 1, there is no v[n-2] so the loop does not run.)
        while (qhat == base || (n > 1 && qhat * v[n - 2] > (rhat * base + u[j + n - 2])))
        {
          qhat -= 1;
          rhat += v[n - 1];
          if (rhat >= base)
            break;
        }

        // D4: Multiply v by qhat and subtract from u[j...j+n].
        // If the subtraction borrows (i.e. the result goes negative),
        // then qhat was one too high; add back v and decrement qhat.

        // Propagate the borrow to the next digit.
        pair<bool, Long> sub = subtract(u, v, qhat, j, base);
        if (sub.first) // borrow produced
        {
          // Correction step: add back v.
          add(u, v, sub.second, j, base);
          qhat -= 1;
        }
        q[j] = qhat;
      }

      // Unnormalize the remainder.
      vector<DataT> r = Divide(u, d, base);

      // Remove potential leading zeros.
      BigIntegerUtil::TrimZeros(q);
      BigIntegerUtil::TrimZeros(r);
      return {q, r};
    }

  private:
    // Subtracts qhat * v from u, starting at index j.
    // Returns true if a borrow is produced (meaning the subtraction went negative).
    static pair<bool, Long> subtract(vector<DataT> &u, vector<DataT> &v, DataT qhat, SizeT j, BaseT base)
    {
      Long borrow = 0;
      SizeT n = v.size();
      for (SizeT i = 0; i < n; i++)
      {
        // Multiply v[i] by qhat and subtract along with the borrow.
        Long p = (Long)(qhat * v[i]);
        Long sub = (Long)u[i + j] - p - (Long)borrow;
        borrow = 0;
        if (sub < 0)
        {
          Long t = (-sub + base - 1) / base; // Compute how many times we need to add base.
          sub += t * base;                   // Adjust sub in one go.
          borrow += t;                       // Increase borrow by that number.
        }
        u[i + j] = sub;
      }

      Long sub1 = (Long)u[j + n] - (Long)borrow;
      if (sub1 >= 0)
        u[j + n] = sub1;

      return {sub1 < 0, borrow};
    }

    // Adds v to u starting at index j (used to undo an oversubtraction).
    static void add(vector<DataT> &u, const vector<DataT> &v, Long borrow, SizeT j, BaseT base)
    {
      DataT carry = 0;
      SizeT n = v.size();
      for (SizeT i = 0; i < n; i++)
      {
        ULong sum = (ULong)u[i + j] + (ULong)v[i] + (ULong)carry;
        u[i + j] = sum % base;
        carry = sum / base;
      }
      u[j + n] = u[j + n] + carry - borrow;
    }
  };
}

#endif
