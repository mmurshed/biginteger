
/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef TOOMCOOK_DATA
#define TOOMCOOK_DATA

#include <vector>
#include <stack>
using namespace std;

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
    inline Long S(SizeT k) const { return r_table[k] + q_table[k+1]; } // s ← q_k + q_k+1
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
    void ComputeTablesAndAllocate(SizeT n, Int R = 2) // R ← 2
    {
      K = 1;                // k ← 1

      Int Q = R * 2;                 // Q ← 4

      // q_0 ← q_1 ← 16
      q_table.push_back( twopow(Q) );
      q_table.push_back( q_table[0] );
      
      // r_0 ← r_1 ← 4
      r_table.push_back( twopow(R) );
      r_table.push_back( r_table[0] );

      ULong p = P(K);
      ULong rSize = 2 * r_table[K] + 1;
      ULong uSize = p * rSize;
      ULong uSizeMax = uSize;
      ULong wSize = 2 * uSize;
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
        Long q = q_table[K];
        Long r = r_table[K];
        if(p < n)
        {
          ULong uSizeTemp = 2 * r * p;
          uSize += uSizeTemp; // 2r * p
          uSizeMax = uSizeTemp + p; // p*(2r+1)
          
          wSizeMax  = 2 * uSizeMax; // 2p*(2r+1)
          wSize += 3*p*(2*r+1); // 3p*(2r+1)
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
