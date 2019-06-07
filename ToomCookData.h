
/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef TOOM_COOK_DATA_H
#define TOOM_COOK_DATA_H

#include <vector>
#include <stack>
using namespace std;

#include "BigIntegerUtil.h"
#include "StackPool.h"

namespace BigMath
{
  class ToomCookData
  {
    public:
    StackPool<DataT> U;
    vector<DataT> UTemp;
    StackPool<DataT> V;
    vector<DataT> VTemp;
    StackPool<DataT> W;
    vector<DataT> WTemp;

    stack<Int> C;
  
    vector<Long> q_table;
    vector<Long> r_table;

    enum {CODE1 = -1, CODE2 = -2, CODE3 = -3};   

    static inline Long sqr(Long n) { return n * n; }
    static inline Long twopow(Int pow) {  return 1L << pow; }

    inline Long P(SizeT k) const { return q_table[k] + q_table[k-1]; } // p ← q_k–1 + q_k
    inline Long Q(SizeT k) const { return q_table[k]; } // q ← q_k
    inline Long R(SizeT k) const { return r_table[k]; } // r ← r_k
    SizeT K;

    ToomCookData(SizeT n)
    {
      ComputeTablesAndAllocate(n);
    }

    // Setp 1. [Compute q, r tables.]
    // In this step we build the sequences
    //
    // k   = 0     1    2    3    4     5     6
    // q_k = 2^4  2^4  2^6  2^8  2^10  2^13  2^16
    // r_k = 2^2  2^2  2^2  2^2  2^3   2^3   2^4
    //
    // The multiplication of 70000-bit numbers would cause 
    // this step to terminate with k = 6, since 70000 < 2^13 + 2^16.)
    void ComputeTablesAndAllocate(SizeT n)
    {
      K = 1;                // k ← 1

      // q_0 ← q_1 ← 16
      q_table.push_back(4);
      q_table.push_back(4);
      
      // r_0 ← r_1 ← 4
      r_table.push_back(2);
      r_table.push_back(2);

      Int Q = 2;                 // Q ← 4
      Int R = 1;                 // R ← 2

      ULong p = P(K);
      ULong rSize = 2 * r_table[K] + 1;
      ULong uSize = p * rSize;
      ULong wSize = 2 * uSize;
      ULong uSizeMax = uSize;
      ULong wSizeMax = wSize;

      // Now if q_k−1 + q_k < n
      // and repeat this operation until q_k−1 + q_k ≥ n. 
      while(p < n)
      {
        K++;                     // k ← k + 1
        Q += R;                  // Q ← Q + R
        
        // R ← ⌊sqrt(Q)⌋
        // Instead of sqrt calculation
        // set R ← R + 1 if (R + 1)^2 ≤ Q 
        // leave R unchanged if (R + 1)^2 > Q.
        if( sqr(R+1) <= Q )
          R++;

        q_table.push_back( twopow(Q) ); // q_k ← 2^Q
        r_table.push_back( twopow(R) ); // r_k ← 2^R

        p = P(K);
        if(p < n)
        {
          rSize = 2 * r_table[K];
          ULong uSizeTemp = rSize * p;
          uSize += uSizeTemp;
          uSizeMax = uSizeTemp + p;
          wSizeMax  = 2 * (uSizeTemp + p);
          wSize += wSizeMax;
        }
      }

      U.Resize(uSize);
      V.Resize(uSize);
      W.Resize(wSize);
      UTemp.resize(uSizeMax);
      VTemp.resize(uSizeMax);
      WTemp.resize(wSizeMax);
    }
   };
}

#endif
