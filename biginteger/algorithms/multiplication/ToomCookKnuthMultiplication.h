
/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef TOOMCOOK_KNUTH_MULTIPLICATION
#define TOOMCOOK_KNUTH_MULTIPLICATION

#include <vector>
#include <stack>
using namespace std;

#include "../../BigInteger.h"
#include "../../common/Util.h"
#include "../Addition.h"
#include "../Subtraction.h"
#include "../multiplication/ClassicMultiplication.h"
#include "../division/KnuthDivision.h"
#include "ToomCookData.h"

namespace BigMath
{
  class ToomCookKnuthMultiplication
  {
    enum
    {
      CODE1 = -1,
      CODE2 = -2,
      CODE3 = -3
    };

  private:
    static vector<DataT> ComputeAnswer(ToomCookData &tcd, SizeT k, BaseT base)
    {
      // Step 7. [Find a’s.]
      // Set r ← r_k
      Long r = tcd.r_table[k];
      Long _2r = 2 * r;
      // q ← q_k
      Long q = tcd.q_table[k];
      // p ← q_k–1 + q_k
      Long p = tcd.q_table[k - 1] + tcd.q_table[k];

      Long wStart = tcd.W.size() - _2r - 1;

      // At this point stack W contains a sequence of numbers
      // ending with W(0), W(1), ..., W(2r) from bottom to
      // top, where each W(j) is a 2p-bit number.
      Long _2p = 2 * p;
      for (Long j = 0; j <= _2r; j++)
      {
        Long wPos = wStart + j;
        TrimZeros(tcd.W[wPos], _2p);
        Resize(tcd.W[wPos], _2p);
      }

      // Now for j = 1, 2, 3, ... , 2r, perform the following
      // Here j must increase and t must decrease.
      for (Long j = 1; j <= _2r; j++)
      {
        // loop: For t = 2r, 2r − 1, 2r – 2, ..., j,
        for (Long t = _2r; t >= j; t--)
        {
          Long wPos = wStart + t;

          // set W(t) ← ( W(t) – W(t − 1) ) / j.
          // The quantity ( W(t) - W(t-1) ) / j will always be a
          // nonnegative integer that fits in 2p bits.

          // W(t) –= W(t − 1)
          SubtractFrom(
              tcd.W[wPos], 0, (SizeT)tcd.W[wPos].size() - 1,
              tcd.W[wPos - 1], 0, (SizeT)tcd.W[wPos - 1].size() - 1,
              base);

          // W(t) /= j
          if (j > 1)
          {
            // vector<DataT> bigJ = Convert(j);
            ClassicDivision::DivideTo(tcd.W[wPos], j, base);
          }
        }
      }

      // Step 8. [Find W’s.]
      // For j = 2r − 1, 2r – 2, ..., 1, perform the following
      // Here j must decrease and t must increase.
      for (Long j = _2r - 1; j >= 1; j--)
      {
        // For t = j, j + 1, ..., 2r − 1,
        for (Long t = j; t < _2r; t++)
        {
          Long wPos = wStart + t;

          // set W(t) ← W(t) – j * W(t + 1).
          // The result of this operation will again be
          // a nonnegative 2p-bit integer.

          // vector<DataT> bigJ = Convert(j);
          // Wt1 = W(t+1) * j
          vector<DataT> Wt1 = ClassicMultiplication::Multiply(
              tcd.W[wPos + 1],
              j,
              base);

          // set W(t) ← W(t) – j * W(t + 1).
          // W(t) ← W(t) – Wt1
          // W(t) -= Wt1
          SubtractFrom(
              tcd.W[wPos], 0, (SizeT)tcd.W[wPos].size() - 1,
              Wt1, 0, (SizeT)Wt1.size() - 1,
              base);
        }
      }

      // Step 9. [Set answer.]
      // Set w to the 2(q_k + q_k+1)-bit integer
      // ( ... (W(2r) * 2^q + W(2r-1)) * 2^q + ... + W(1)) * 2^q + W(0)
      // Long wsize = 2 * (q_table[k] + q_table[k+1]);
      Long wsize = _2p + _2r * q;
      vector<DataT> w(wsize, 0);

      SizeT wAnsStart = wsize - _2p;
      SizeT wAnsEnd = wsize - 1;

      for (Long i = _2r; i >= 0; i--)
      {
        Long wPos = wStart + i;
        AddTo(
            w, wAnsStart, wAnsEnd,
            tcd.W[wPos], 0, tcd.W[wPos].size() - 1,
            base);

        wAnsStart -= q;
        wAnsEnd = wAnsStart + _2p - 1;
      }

      // Remove W(2r), . . . , W(0) from stack W.
      for (Long i = _2r; i >= 0; i--)
      {
        tcd.W.pop_back();
      }

      return w;
    }

    static vector<DataT> MultiplyRPlus1Part(vector<DataT> U, DataT j, Long p, Long q, Long r, BaseT base)
    {
      vector<DataT> Uj(p, 0);

      if (j == 0)
      {
        Copy(U, 0, q - 1, Uj, 0, p - 1);
        return Uj;
      }

      Long uEnd = U.size() - 1;
      Long uStart = uEnd - q + 1;

      // Compute the p-bit numbers
      // Given U =  U_r ... U_1 U_0
      // ( ... (U_r * j + U_r-1) * j + ... + U_1) * j + U_0
      for (Long i = r; i >= 0; i--)
      {
        AddTo(
            Uj, 0, p - 1,
            U, uStart, uEnd,
            base);

        // Don't multiply if it's U_0
        // Or if Uj is already zero
        if (i > 0 && !IsZero(Uj))
        {
          // Case of j = zero is already elminiated
          ClassicMultiplication::MultiplyTo(Uj, j, base);
        }

        uEnd -= q;
        uStart = uEnd - q + 1;

        TrimZeros(Uj, p);
        Resize(Uj, p);
      }

      return Uj;
    }

    static void BreakIntoRPlus1Parts(ToomCookData &tcd, Long p, Long q, Long r, BaseT base)
    {
      // Step 4. [Break into r + 1 parts.]
      // Let the number at the top of stack C be regarded as a list
      // of r + 1 numbers with q bits each, (U_r . . . U_1U_0)_2q .

      // The top of stack C now contains an (r + 1)q = (q_k + q_k+1)-bit number.
      // Remove U_r ... U_1 U_0 from stack C.
      vector<DataT> U = tcd.Cval[tcd.C.top()];
      tcd.C.pop();
      tcd.Cval.pop_back();

      // Now the top of stack C contains another list of
      // r + 1 q-bit numbers, V_r ... V_1 V_0.
      // Remove V_r ... V_1 V_0 from stack C.
      vector<DataT> V = tcd.Cval[tcd.C.top()];
      tcd.C.pop();
      tcd.Cval.pop_back();

      Long _2r = 2 * r;

      // For j = 2r, 2r-1, ... , 1, 0
      for (Long j = _2r; j >= 0; j--)
      {
        // code
        tcd.C.push(j == _2r ? CODE2 : CODE3);

        // vector<DataT> jbig = Convert(j);

        // Compute the p-bit numbers
        // ( ... (V_r * j + V_r-1) * j + ... + V_1) * j + V_0
        vector<DataT> Vj = MultiplyRPlus1Part(V, j, p, q, r, base);
        tcd.Cval.push_back(Vj);
        tcd.C.push((SizeT)tcd.Cval.size() - 1);

        // Compute the p-bit numbers
        // ( ... (U_r * j + U_r-1) * j + ... + U_1) * j + U_0
        // and successively put these values onto stack U.
        vector<DataT> Uj = MultiplyRPlus1Part(U, j, p, q, r, base);
        tcd.Cval.push_back(Uj);
        tcd.C.push((SizeT)tcd.Cval.size() - 1);
      }
      // Stack C contains
      // code-2, V(2r), U(2r),
      // code-3, V(2r-1), U(2r-1),
      // ... ,
      // code-3, V(1), U(1),
      // code-3, V(0), U(0)
    }

    static vector<DataT> DivideRecursive(ToomCookData &tcd, SizeT k, BaseT base)
    {
      // Step 3. [Check recursion level.]
      // Decrease k by 1.
      --k;

      // If k = 0, the top of stack C now contains
      // two 32-bit numbers, u and v
      if (k == 0)
      {
        // remove them
        vector<DataT> u = tcd.Cval[tcd.C.top()];
        tcd.C.pop();
        tcd.Cval.pop_back();

        vector<DataT> v = tcd.Cval[tcd.C.top()];
        tcd.C.pop();
        tcd.Cval.pop_back();
        // set w ← uv using a built-in routine for multiplying
        // 32-bit numbers, and go to step 10.
        vector<DataT> w = ClassicMultiplication::Multiply(u, v, base);
        TrimZeros(w);
        return w;
      }

      // If k > 0
      // set r ← r_k
      Long r = tcd.r_table[k];
      // q ← q_k
      Long q = tcd.q_table[k];
      // p ← q_k–1 + q_k
      Long p = tcd.q_table[k - 1] + tcd.q_table[k];

      // go on to step 4.
      BreakIntoRPlus1Parts(tcd, p, q, r, base);

      // Step 5. [Recurse.]
      // Go back to step 3.
      return DivideRecursive(tcd, k, base);
    }

    static vector<DataT> Divide(ToomCookData &tcd, SizeT k, BaseT base)
    {
      // While k > 0
      while (--k > 0)
      {
        Long r = tcd.r_table[k];                      // r ← r_k
        Long q = tcd.q_table[k];                      // q ← q_k
        Long p = tcd.q_table[k - 1] + tcd.q_table[k]; // p ← q_k–1 + q_k

        BreakIntoRPlus1Parts(tcd, p, q, r, base);
      }

      // At this point k is 0
      // the top of stack C now contains
      // two 32-bit numbers, u and v
      // remove them
      vector<DataT> u = tcd.Cval[tcd.C.top()];
      tcd.C.pop();
      tcd.Cval.pop_back();

      vector<DataT> v = tcd.Cval[tcd.C.top()];
      tcd.C.pop();
      tcd.Cval.pop_back();
      // set w ← uv using a built-in routine for multiplying
      // 32-bit numbers, and go to step 10.
      vector<DataT> w = ClassicMultiplication::Multiply(u, v, base);
      TrimZeros(w);
      return w;
    }

    static vector<DataT> Multiply(ToomCookData &tcd, SizeT k, BaseT base)
    {
      Int code = CODE3;

      // Step 10. [Return.]
      while (code != CODE1)
      {
        // If it is code-3, go to step T6.
        if (code == CODE3)
        {
          // Step 6. [Save one product.]
          // Go back to step T3.
          vector<DataT> w = Divide(tcd, k, base);
          tcd.W.push_back(w);
          k = 0;
        }
        // If it is code-2, put w onto stack W and go to step 7.
        else if (code == CODE2)
        {
          vector<DataT> w = ComputeAnswer(tcd, k, base);
          tcd.W.push_back(w);
        }

        // Set k ← k + 1.
        k++;
        code = tcd.C.top();
        tcd.C.pop();
      }

      // And if it is code-1, terminate the algorithm (w is the answer).
      return tcd.W.at(0);
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
        BaseT base)
    {
      if (IsZero(a) || IsZero(b)) // 0 times
        return vector<DataT>();

      // If b is a single digit, use the scalar multiplication
      if (b.size() == 1)
        return ClassicMultiplication::Multiply(a, b[0], base);
      // If a is a single digit, use the scalar multiplication
      if (a.size() == 1)
        return ClassicMultiplication::Multiply(b, a[0], base);

      SizeT n = (SizeT)max(a.size(), b.size());

      ToomCookData tcd(n);

      // Set stacks U, V, C, and W empty.
      // Stack C contains control code.
      // code-1 as -1, code-2 as -2, and code-3 as -3.
      // It also contains the locations of Cval as positive integers.

      // Step 2. [Put u, v on stack.]
      Long p = tcd.q_table[tcd.k - 1] + tcd.q_table[tcd.k];

      // u and v as numbers of exactly q_k-1 + q_k bits each.
      vector<DataT> u(p, 0);
      vector<DataT> v(p, 0);

      Copy(a, u);
      Copy(b, v);

      // Put code-1 on stack C
      tcd.C.push(CODE1);

      // Place v onto stack C
      tcd.Cval.push_back(v);
      tcd.C.push(0);

      // Place u onto stack C
      tcd.Cval.push_back(u);
      tcd.C.push(1);

      vector<DataT> w = Multiply(tcd, tcd.k, base);

      TrimZeros(w);

      return w;
    }
  };
}

#endif
