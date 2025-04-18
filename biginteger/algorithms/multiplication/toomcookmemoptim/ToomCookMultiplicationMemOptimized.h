
/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef TOOMCOOK_MULTIPLICATION_MEM_OPTIMIZED
#define TOOMCOOK_MULTIPLICATION_MEM_OPTIMIZED

#include <vector>
#include <stack>
using namespace std;

#include "ToomCookDataMemOptimized.h"
#include "../../../BigInteger.h"
#include "../../../common/Util.h"
#include "../../Addition.h"
#include "../../Subtraction.h"
#include "../../multiplication/ClassicMultiplication.h"
#include "../../division/KnuthDivision.h"

namespace BigMath
{
  class ToomCookMultiplicationMemOptimized
  {
  private:
    static void ComputeAnswer(ToomCookDataMemOptim &tcd, SizeT k, BaseT base)
    {
      // Step 7. [Find a’s.]
      Long _2r = 2 * tcd.R(k);
      Long q = tcd.Q(k);
      Long p = tcd.P(k);

      Long wStart = tcd.W.ranges.size() - _2r - 1;

      // At this point stack W contains a sequence of numbers
      // ending with W(0), W(1), ..., W(2r) from bottom to
      // top, where each W(j) is a 2p-bit number.
      Long _2p = 2 * p;

      // Now for j = 1, 2, 3, ... , 2r, perform the following
      // Here j must increase and t must decrease.
      for (Long j = 1; j <= _2r; j++)
      {
        // loop: For t = 2r, 2r − 1, 2r – 2, ..., j,
        for (Long t = _2r; t >= j; t--)
        {
          Range wtr = tcd.W.ranges[wStart + t];      // W[t]
          Range wt1r = tcd.W.ranges[wStart + t - 1]; // W[t-1]

          // set W(t) ← ( W(t) – W(t − 1) ) / j.
          // The quantity ( W(t) - W(t-1) ) / j will always be a
          // nonnegative integer that fits in 2p bits.

          // W(t) –= W(t − 1)
          SubtractFrom(
              tcd.W.data, wtr.first, wtr.second,
              tcd.W.data, wt1r.first, wt1r.second,
              base);

          // W(t) /= j
          if (j > 1)
          {
            ClassicDivision::DivideTo(
                tcd.W.data, wtr.first, wtr.second,
                (DataT)j,
                base);
          }
        }
      }

      // Step 8. [Find W’s.]
      // For j = 2r − 1, 2r – 2, ..., 1, perform the following
      // Here j must decrease and t must increase.
      SizeT wEnd = _2p;
      for (Long j = _2r - 1; j >= 1; j--)
      {
        Range wtr = tcd.W.ranges[wStart + j]; // W[j]

        // For t = j, j + 1, ..., 2r − 1,
        for (Long t = j; t < _2r; t++)
        {
          Range wtr = tcd.W.ranges[wStart + t];      // W[t]
          Range wt1r = tcd.W.ranges[wStart + t + 1]; // W[t+1]

          // set W(t) ← W(t) – j * W(t + 1).
          // The result of this operation will again be
          // a nonnegative 2p-bit integer.

          // Wt1 = W(t+1) * j
          wEnd = ClassicMultiplication::Multiply(
              tcd.W.data, wt1r.first, wt1r.second,
              (ULong)j,
              tcd.WTemp, 0, wEnd,
              base);

          // set W(t) ← W(t) – j * W(t + 1).
          // W(t) ← W(t) – Wt1
          // W(t) -= Wt1
          SubtractFrom(
              tcd.W.data, wtr.first, wtr.second,
              tcd.WTemp, 0, wEnd,
              base);

          SetBit(tcd.WTemp, 0, tcd.WTemp.size() - 1, 0);
        }
      }

      // Step 9. [Set answer.]
      // Set w to the 2(q_k + q_k+1)-bit integer
      // ( ... (W(2r) * 2^q + W(2r-1)) * 2^q + ... + W(1)) * 2^q + W(0)

      // Long wsize = 2 * tcd.S(k);
      Long wsize = _2p + _2r * q;

      SizeT wAnsStart = wsize - _2p;
      SizeT wAnsEnd = wsize - 1;

      for (Long i = _2r; i >= 0; i--)
      {
        Range wtr = tcd.W.Pop(); // W[i]
        AddTo(
            tcd.WTemp, wAnsStart, wsize - 1,
            tcd.W.data, wtr.first, wtr.second,
            base);

        SetBit(tcd.W.data, wtr.first, wtr.second, 0);

        wAnsStart -= q;
        wAnsEnd = wAnsStart + _2p - 1;
      }

      // Remove W(2r), . . . , W(0) from stack W.
      tcd.W.Push(tcd.WTemp, wsize);
    }

    static void MultiplyRPlus1Parts(
        vector<DataT> const &u, Range ur,
        DataT j,
        vector<DataT> &result, Range rRes,
        Long q, Long r, BaseT base)
    {
      if (j == 0)
      {
        Copy(
            u, ur.first, ur.second,
            result, rRes.first, rRes.second);
      }

      SetBit(
          result, rRes.first, rRes.second,
          0);

      Long uEnd = ur.second;
      Long uStart = uEnd - q + 1;

      // Compute the p-bit numbers
      // Given U =  U_r ... U_1 U_0
      // ( ... (U_r * j + U_r-1) * j + ... + U_1) * j + U_0
      for (Long i = r; i >= 0; i--)
      {
        AddTo(
            result, rRes.first, rRes.second,
            u, uStart, uEnd,
            base);

        // Don't multiply if it's U_0
        // Or if result is already zero
        if (i > 0 && !IsZero(result, rRes.first, rRes.second))
        {
          // Case of j = zero is already eliminated
          ClassicMultiplication::MultiplyTo(
              result, rRes.first, rRes.second,
              j,
              base);
        }

        uEnd -= q;
        uStart = uEnd - q + 1;
      }
    }

    static void BreakIntoRPlus1Parts(ToomCookDataMemOptim &tcd, Long k, BaseT base)
    {
      Long p = tcd.P(k);
      Long q = tcd.Q(k);
      Long r = tcd.R(k);
      Long _2r = 2 * r;

      // Step 4. [Break into r + 1 parts.]
      // Let the number at the top of stack C be regarded as a list
      // of r + 1 numbers with q bits each, (U_r . . . U_1U_0)_2q .

      // The top of stack C now contains an (r + 1)q = (q_k + q_k+1)-bit number.
      // Remove U_r ... U_1 U_0 from stack C.
      Range ur = tcd.U.Pop();

      // Now the top of stack C contains another list of
      // r + 1 q-bit numbers, V_r ... V_1 V_0.
      // Remove V_r ... V_1 V_0 from stack C.
      Range vr = tcd.V.Pop();

      Range rRange = {0, p - 1};

      // For j = 2r, 2r-1, ... , 1, 0
      for (Long j = _2r; j >= 0; j--)
      {
        // code
        tcd.C.push(j == _2r ? ToomCookDataMemOptim::CODE2 : ToomCookDataMemOptim::CODE3);

        // Compute the p-bit numbers
        // ( ... (V_r * j + V_r-1) * j + ... + V_1) * j + V_0
        MultiplyRPlus1Parts(
            tcd.V.data, vr,
            j,
            tcd.VTemp, rRange,
            q, r, base);

        tcd.U.Push(ur.first + rRange.first, ur.first + rRange.second);

        // Compute the p-bit numbers
        // ( ... (U_r * j + U_r-1) * j + ... + U_1) * j + U_0
        // and successively put these values onto stack U.
        MultiplyRPlus1Parts(
            tcd.U.data, ur,
            j,
            tcd.UTemp, rRange,
            q, r, base);

        tcd.V.Push(vr.first + rRange.first, vr.first + rRange.second);

        rRange.first += p;
        rRange.second = rRange.first + p - 1;
      }

      copy(
          tcd.UTemp.begin(),
          tcd.UTemp.begin() + rRange.first,
          tcd.U.data.begin() + ur.first);
      copy(
          tcd.VTemp.begin(),
          tcd.VTemp.begin() + rRange.first,
          tcd.V.data.begin() + vr.first);

      // Stack C contains
      // code-2, V(2r), U(2r),
      // code-3, V(2r-1), U(2r-1),
      // ... ,
      // code-3, V(1), U(1),
      // code-3, V(0), U(0)
    }

    static void BreakIntoParts(ToomCookDataMemOptim &tcd, SizeT k, BaseT base)
    {
      while (--k > 0)
      {
        BreakIntoRPlus1Parts(tcd, k, base);
      }

      // At this point k is 0
      // the top of stack C now contains
      // two 32-bit numbers, u and v
      // remove them
      Range ur = tcd.U.Pop();
      Range vr = tcd.V.Pop();

      // set w ← uv using a built-in routine for multiplying
      // 32-bit numbers, and go to step 10.
      SizeT wStart = tcd.W.ranges.empty() ? 0 : tcd.W.Top().second + 1;
      SizeT wEnd = ClassicMultiplication::Multiply(
          tcd.U.data, ur.first, ur.second,
          tcd.V.data, vr.first, vr.second,
          tcd.W.data, wStart,
          base);

      tcd.W.Push(wStart, wEnd);
    }

    static vector<DataT> Multiply(ToomCookDataMemOptim &tcd, BaseT base)
    {
      SizeT k = tcd.K;
      Int code = ToomCookDataMemOptim::CODE3;

      // Step 10. [Return.]
      while (code != ToomCookDataMemOptim::CODE1)
      {
        // If it is code-3, go to step T6.
        if (code == ToomCookDataMemOptim::CODE3)
        {
          // Step 6. [Save one product.]
          // Go back to step T3.
          BreakIntoParts(tcd, k, base);
          k = 0;
        }
        // If it is code-2, put w onto stack W and go to step 7.
        else if (code == ToomCookDataMemOptim::CODE2)
        {
          ComputeAnswer(tcd, k, base);
        }

        k++; // k ← k + 1
        code = tcd.C.top();
        tcd.C.pop();
      }

      // And if it is code-1, terminate the algorithm (w is the answer).
      return tcd.W.data;
    }

  public:
    /*
     * Given a positive integer n and two nonnegative n-bit integers
     * u and v, this algorithm forms their 2n-bit product, w.
     *
     * Four auxiliary stacks are used to hold the long numbers that
     * are manipulated during the procedure:
     *
     * Stack U, V: Temporary storage of U(j) and V(j) in step 4.
     * Stack C: Numbers to be multiplied, and control codes.
     * Stack W: Storage W(j).
     *
     * These stacks may contain either binary numbers or special
     * control symbols called code-1, code-2, and code-3. The
     * algorithm also constructs an auxiliary table of numbers qk, rk;
     * this table is maintained in such a manner that it may be stored
     * as a linear list, where a single pointer that traverses the list
     * (moving back and forth) can be used to access the current table
     * entry of interest.
     */
    static vector<DataT> Multiply(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base,
        bool trim = true)
    {
      SizeT n = (SizeT)max(a.size(), b.size());
      ToomCookDataMemOptim tcd(n);

      // Step 2. [Put u, v on stack.]

      // u and v as numbers of exactly q_k-1 + q_k bits each.
      ULong p = tcd.P(tcd.K);
      tcd.U.Push(a, p);
      tcd.V.Push(b, p);

      // Put code-1 on stack C
      tcd.C.push(ToomCookDataMemOptim::CODE1);

      vector<DataT> w = Multiply(tcd, base);

      if (trim)
      {
        TrimZeros(w);
      }

      return w;
    }
  };
}

#endif
