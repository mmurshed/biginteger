#ifndef MONTGOMERY_DIVISION
#define MONTGOMERY_DIVISION

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "../../BigIntegerBuilder.h"
#include "../../BigIntegerUtil.h"
#include "../../ops/BigIntegerIO.h"
#include "../../ops/MultiplicationStrategy.h"
#include "../CommonAlgorithms.h"
#include "../addition/ClassicAddition.h"
#include "../subtraction/ClassicSubtraction.h"
#include "../divison/KnuthDivision.h"

using namespace std;

namespace BigMath
{
    class MontgomeryDivision
    {
        // Helper function to compute the greatest common divisor (GCD)
        static DataT gcd(DataT a, DataT b)
        {
            while (b != 0)
            {
                DataT temp = b;
                b = a % b;
                a = temp;
            }
            return a;
        }

        public:
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            if (BigIntegerUtil::IsZero(b))
            {
                throw invalid_argument("Division by zero");
            }

            int m = b.size();
            if (m == 0)
            {
                return {vector<DataT>(), a};
            }

            DataT b0 = b[0];
            if (gcd(b0, base) != 1)
            {
                // Fall back to classic division if b0 and base are not coprime
                return KnuthDivision::DivideAndRemainder(a, b, base);
            }

            // Compute mu = -b0^{-1} mod base
            DataT mu = 0;
            bool found = false;
            for (SizeT i = 1; i < base; ++i)
            {
                if ((b0 * i) % base == 1)
                {
                    mu = (base - i) % base;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                return KnuthDivision::DivideAndRemainder(a, b, base);
            }

            vector<DataT> r = a;
            for (int i = 0; i < m; ++i)
            {
                DataT r_i = i < r.size() ? r[i] : 0;
                DataT q_i = (r_i * mu) % base;

                vector<DataT> q_i_vec = {q_i};
                vector<DataT> qb = MultiplicationStrategy::MultiplyUnsigned(b, q_i_vec, base);

                // Shift left by i digits (equivalent to multiplying by base^i)
                vector<DataT> shifted_qb(i, 0);
                shifted_qb.insert(shifted_qb.end(), qb.begin(), qb.end());

                r = ClassicAddition::Add(r, shifted_qb, base);

                // Divide t by base (remove the least significant digit)
                if (!r.empty())
                {
                    r.erase(r.begin());
                }
                else
                {
                    r = {0};
                }
            }

            // After reduction, ensure t < b
            while (BigIntegerComparator::Compare(r, b) >= 0)
            {
                r = ClassicSubtraction::Subtract(r, b, base);
            }

            BigIntegerUtil::TrimZeros(r);

            // Compute quotient: (a - remainder) / b
            vector<DataT> a_minus_r = ClassicSubtraction::Subtract(a, r, base);
            auto qr = KnuthDivision::DivideAndRemainder(a_minus_r, b, base);
            BigIntegerUtil::TrimZeros(qr.first);

            return {qr.first, r};
        }
    };
}
#endif