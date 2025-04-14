/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef TOOMCOOK_DATA
#define TOOMCOOK_DATA

#include <vector>
#include <stack>
using namespace std;

#include "../../common/Util.h"

namespace BigMath
{

    class ToomCookData
    {
    public:
        vector<vector<DataT>> W;
        stack<Int> C;
        vector<vector<DataT>> Cval;

        vector<Long> q_table;
        vector<Long> r_table;

        SizeT k;

        ToomCookData(SizeT n)
        {
            k = ComputeTables(n);
        }

    private:
        // Setp 1. [Compute q, r tables.]
        // In this step we build the sequences
        //
        // k   = 0     1    2    3    4     5     6
        // q_k = 2^4  2^4  2^6  2^8  2^10  2^13  2^16
        // r_k = 2^2  2^2  2^2  2^2  2^3   2^3   2^4
        //
        // The multiplication of 70000-bit numbers would cause
        // this step to terminate with k = 6, since 70000 < 2^13 + 2^16.)
        SizeT ComputeTables(SizeT n)
        {
            SizeT k = 1; // k ← 1

            // q_0 ← q_1 ← 16
            q_table.push_back(16);
            q_table.push_back(16);

            // r_0 ← r_1 ← 4
            r_table.push_back(4);
            r_table.push_back(4);

            Int Q = 4; // Q ← 4
            Int R = 2; // R ← 2

            // Now if q_k−1 + q_k < n
            // and repeat this operation until q_k−1 + q_k ≥ n.
            while (q_table[k - 1] + q_table[k] < n)
            {
                k++;    // k ← k + 1
                Q += R; // Q ← Q + R

                // R ← ⌊sqrt(Q)⌋
                // Instead of sqrt calculation
                // set R ← R + 1 if (R + 1)^2 ≤ Q
                // leave R unchanged if (R + 1)^2 > Q.
                if (sqr(R + 1) <= Q)
                    R++;

                q_table.push_back(twopow(Q)); // q_k ← 2^Q
                r_table.push_back(twopow(R)); // r_k ← 2^R
            }

            return k;
        }

        static inline Long sqr(Long n) { return n * n; }
        static inline Long twopow(Int pow) { return 1L << pow; }
    };
}

#endif