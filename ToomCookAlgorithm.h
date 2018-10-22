
/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef TOOM_COOK_ALGORITHM_H
#define TOOM_COOK_ALGORITHM_H

#include <vector>
#include <stack>
using namespace std;

#include "BigInteger.h"
#include "BigIntegerUtil.h"
#include "ClassicalAlgorithms.h"

namespace BigMath
{
  class ToomCookAlgorithm
  {
    public:  
    ToomCookAlgorithm()
    {
      CODE1 = -1;
      CODE2 = -2;
      CODE3 = -3;
    }
    private:
    vector< vector<DataT> > W;
    stack<Int> C; 
    vector< vector<DataT> > Cval;

    vector<Long> q_table;
    vector<Long> r_table;
    
    Int CODE1;
    Int CODE2;
    Int CODE3;

    public:
    vector<DataT> MultiplyRPart(vector<DataT> U, SizeT r, SizeT j, SizeT s, ULong base)
    {
      vector<DataT> Uj(s + 1);
      BigIntegerUtil::SetBit(Uj, 0, U.size() - 1, 0);

      Int p = U.size() / r;
      Int q = U.size() % r;
      Int k = 0;
      SizeT uEnd = 1;

      // Compute the p-bit numbers
      // Given U =  U_r ... U_1 U_0
      // ( ... (U_r * j + U_r-1) * j + ... + U_1) * j + U_0 
      for(Int i = U.size() - 1; i > 0; i = k - 1)
      {
        k = i - p + 1;
        // If the digits couldn't be equally divided
        // include one extra digit until all the extra
        //  digits (q), are accounted for
        if(q > 0)
        {
          k--;
          q--;
        }
        
        ClassicalAlgorithms::AddToUnsigned(
          Uj, 0, uEnd,
          U, k, i,
          base);
        
        // Don't multiply if it's U_0
        if(k > 0)
        {
          ClassicalAlgorithms::MultiplyToUnsigned(Uj, j, base);
        }
        uEnd = BigIntegerUtil::FindNonZeroByte(Uj);
      }

      BigIntegerUtil::TrimZeros(Uj);

      return Uj;
    }

    vector<DataT> ComputeCode2(SizeT k, ULong base)
    {
        // Step 7. [Find a’s.]
        // Set r ← r_k
        Long r = r_table[k];
        Long _2r = 2 * r;
        // q ← q_k
        Long q = q_table[k];
        // p ← q_k–1 + q_k
        Long p = q_table[k-1] + q_table[k];

        // At this point stack W contains a sequence of numbers
        // ending with W(0), W(1), ..., W(2r) from bottom to 
        // top, where each W(j) is a 2p-bit number.

        // Now for j = 1, 2, 3, ... , 2r, perform the following 
        // Here j must increase and t must decrease.
        for(Long j = 1; j <= _2r; j++)
        {
          // loop: For t = 2r, 2r − 1, 2r – 2, ..., j,  
          for(Long t = _2r; t >= j; t--)
          {
            // set W(t) ← ( W(t) – W(t − 1) ) / j.
            // The quantity ( W(t) - W(t-1) ) / j will always be a 
            // nonnegative integer that fits in 2p bits.

            // W(t) –= W(t − 1)
            ClassicalAlgorithms::SubtractFromUnsigned(
              W[t], 0, W[t].size() - 1,
              W[t-1], 0, W[t-1].size() - 1,
              base);

            // W(t) /= j
            ClassicalAlgorithms::DivideToUnsigned(W[t], j, base);
          }
        }

        // Step 8. [Find W’s.]        
        // For j = 2r − 1, 2r – 2, ..., 1, perform the following 
        // Here j must decrease and t must increase. 
        for(Long j = _2r - 1; j >= 1; j--)
        {
          // loop: For t = j, j + 1, ..., 2r − 1, 
          for(Long t = j; t < _2r; t++)
          {
            // set W(t) ← W(t) – j W(t + 1).
            // The result of this operation will again be 
            // a nonnegative 2p-bit integer.

            // Wt1 = j * W(t+1)
            vector<DataT> Wt1 = ClassicalAlgorithms::MultiplyUnsigned(W[t+1], j, base);

            // W(t) -= Wt1
            ClassicalAlgorithms::SubtractFromUnsigned(
              W[t], 0, W[t].size() - 1,
              Wt1, 0, Wt1.size() - 1,
              base);
          }
        }

        // Step 9. [Set answer.]
        // Set w to the 2(q_k + q_k+1)-bit integer
        // ( ... (W(2r) * 2^q + W(2r-1)) * 2^q + ... + W(1)) * 2^q + W(0)
        Long qp = twopow(q);

        vector<DataT> w(2 * (q_table[k] + q_table[k+1]) );
        BigIntegerUtil::SetBit(w, 0, w.size() - 1, 0);

        for(Int i = _2r; i >= 0; i--)
        {
          ClassicalAlgorithms::AddToUnsigned(
            w, 0, w.size() - 1,
            W[i], 0, W[i].size() - 1,
            base);
          // Don't multiply if it's W_0
          if(i > 0)
          {
            ClassicalAlgorithms::MultiplyToUnsigned(w, qp, base);
          }
        }

        // Remove W(2r), . . . , W(0) from stack W.
        W.clear();

        return w;
    }

    void BreakIntoR1Parts(Long r, Long p, ULong base)
    {
      // Step 4. [Break into r + 1 parts.]
      // Let the number at the top of stack C be regarded as a list 
      // of r + 1 numbers with q bits each, (U_r . . . U_1U_0)_2q . 
      
      // The top of stack C now contains an (r + 1)q = (q_k + q_k+1)-bit number.
      // Remove U_r ... U_1 U_0 from stack C.
      vector<DataT> U = Cval[C.top()];
      C.pop();
      
      // Now the top of stack C contains another list of
      // r + 1 q-bit numbers, V_r ... V_1 V_0.
      // Remove V_r ... V_1 V_0 from stack C.
      vector<DataT> V = Cval[C.top()];
      C.pop();

      Int _2r = 2 * r;
      
      // For j = 0, 1, ... , 2r
      for(Int j = _2r; j >= 0; j--)
      {
        // code
        C.push( j == _2r ? CODE2 : CODE3 );

        // Computer the p-bit numbers
        // ( ... (V_r * j + V_r-1) * j + ... + V_1) * j + V_0 
        Cval.push_back(MultiplyRPart(V, r + 1, j, p, base));
        C.push(Cval.size() - 1);

        // Compute the p-bit numbers
        // ( ... (U_r * j + U_r-1) * j + ... + U_1) * j + U_0 
        // and successively put these values onto stack U. 
        Cval.push_back(MultiplyRPart(U, r + 1, j, p, base));
        C.push(Cval.size() - 1);
        
      }
      // Stack C contains
      // code-2, V(2r), U(2r),
      // code-3, V(2r-1), U(2r-1),
      // ... ,
      // code-3, V(1), U(1),
      // code-3, V(0), U(0)
    }

    vector<DataT> MultiplyUnsignedRecursive(SizeT k, ULong base)
    {
      // Step 3. [Check recursion level.]
      // Decrease k by 1.
      --k;

      // If k = 0, the top of stack C now contains
      // two 32-bit numbers, u and v
      if(k == 0)
      {
        // remove them
        vector<DataT> u = Cval[C.top()];
        C.pop();
        vector<DataT> v = Cval[C.top()];
        C.pop();
        // set w ← uv using a built-in routine for multiplying 
        // 32-bit numbers, and go to step 10.
        return ClassicalAlgorithms::MultiplyUnsigned(u, v, base);
      }
      
      // If k > 0
      // set r ← r_k
      Long r = r_table[k];
      // q ← q_k
      Long q = q_table[k];
      // p ← q_k–1 + q_k
      Long p = q_table[k-1] + q_table[k];

      // go on to step 4.
      BreakIntoR1Parts(r, p, base);     

      // Step 5. [Recurse.]
      // Go back to step 3.
      vector<DataT> w = MultiplyUnsignedRecursive(k, base);
      
      // Step 10. [Return.]
      Int code = C.top();
      C.pop();

      // Set k ← k + 1.
      k++;
  
      // If it is code-3, go to step T6.
      if(code == CODE3)
      {
        // Step 6. [Save one product.]
        // At this point the multiplication algorithm has set w 
        // to one of the products W(j) = U(j)V(j).
        
        // Put w onto stack W.
        // This number w contains 2(q_k + q_k−1) bits.
        W.push_back(w);

        // Go back to step T3.
        w = MultiplyUnsignedRecursive(k, base);
      }
        // If it is code-2, put w onto stack W and go to step 7.
      else if(code == CODE2)
      {
        W.push_back(w);
        w = ComputeCode2(k, base);
      }
        
      // And if it is code-1, terminate the algorithm (w is the answer).
      return w;
    }

    // Setp 1. [Compute q, r tables.]
    SizeT ComputeQRTable(SizeT n)
    {
      // Set k ← 1
      SizeT k = 1;

      // q_0 ← q_1 ← 16
      q_table.push_back(16);
      q_table.push_back(16);
      
      // r_0 ← r_1 ← 4
      r_table.push_back(4);
      r_table.push_back(4);

      // Q ← 4
      Long Q = 4;
      // R ← 2
      Long R = 2;

      // Now if q_k−1 + q_k < n
      // and repeat this operation until q_k−1 + q_k ≥ n. 
      while(q_table[k-1] + q_table[k] < n)
      {
        // Set k ← k + 1
        k++;
        // Q ← Q + R
        Q += R;
        
        // R ← ⌊sqrt(Q)⌋
        // (Note: The calculation of R ← ⌊sqrt(Q)⌋ does not 
        // require a square root to be taken, since we may 
        // simply set R ← R + 1 if (R + 1)^2 ≤ Q and leave 
        // R unchanged if (R + 1)^2 > Q.
        if( sqr(R+1) <= Q )
          R++;

        // q_k ← 2^Q
        q_table.push_back( twopow(Q) );
        // r_k ← 2^R
        r_table.push_back( twopow(R) );
        
        // In this step we build the sequences
        //
        // k   = 0     1    2    3    4     5     6
        // q_k = 2^4  2^4  2^6  2^8  2^10  2^13  2^16
        // r_k = 2^2  2^2  2^2  2^2  2^3   2^3   2^4
        //
        // The multiplication of 70000-bit numbers would cause 
        // this step to terminate with k = 6, since 70000 < 2^13 + 2^16.)
      }

      return k;
    }

    static inline Long sqr(Long n) { return n * n; }
    static inline Long twopow(Int pow) {  return 1L << pow; }

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
    vector<DataT> MultiplyUnsigned(vector<DataT> const& a, vector<DataT> const& b, ULong base)
    {
      SizeT n = max(a.size(), b.size());
      vector<DataT> u(a);
      vector<DataT> v(b);

      BigIntegerUtil::MakeSameSize(u, v);

      SizeT k = ComputeQRTable(n);

      // Set stacks U, V, C, and W empty.
      // Stack C contains control code.
      // code-1 as -1, code-2 as -2, and code-3 as -3.
      // It also contains the locations of Cval as positive integers.

      // Step 2. [Put u, v on stack.]
      // Put code-1 on stack C, then place u and v onto stack C as 
      // numbers of exactly q_k−1 + q_k bits each.

      C.push(CODE1);
      Cval.push_back(u);
      C.push(0);
      Cval.push_back(v);
      C.push(1);

      return MultiplyUnsignedRecursive(k, base);
    }
   };
}

#endif