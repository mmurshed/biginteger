/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef CLASSIC_DIVISION
#define CLASSIC_DIVISION

#include <cstddef>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;

#include "../../BigInteger.h"
#include "../../common/Util.h"
#include "../../algorithms/Multiplication.h"
#include "../Shift.h"
#include "../Addition.h"
#include "../Subtraction.h"

namespace BigMath
{
    class ClassicDivision
    {
    public:
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            if (IsZero(b))
            {
                throw invalid_argument("Division by zero");
            }

            int cmp = Compare(a, b);
            if (cmp < 0)
            {
                return {vector<DataT>(), a};
            }
            else if (cmp == 0)
            {
                return {vector<DataT>{1}, vector<DataT>()};
            }

            vector<DataT> remainder = a;
            vector<DataT> quotient;

            int max_d = a.size() - b.size();
            for (int d = max_d; d >= 0; --d)
            {
                vector<DataT> shifted_divisor;
                shifted_divisor.insert(shifted_divisor.end(), d, 0);
                shifted_divisor.insert(shifted_divisor.end(), b.begin(), b.end());

                if (Compare(shifted_divisor, remainder) > 0)
                {
                    if (quotient.size() > static_cast<size_t>(d))
                    {
                        quotient[d] = 0;
                    }
                    continue;
                }

                DataT q_low = 1, q_high = base - 1;
                DataT q = 0;
                while (q_low <= q_high)
                {
                    DataT candidate = (q_low + q_high) / 2;
                    vector<DataT> candidate_vec = {candidate};
                    vector<DataT> product = Multiply(shifted_divisor, candidate_vec, base);

                    int cmp = Compare(product, remainder);
                    if (cmp <= 0)
                    {
                        q = candidate;
                        q_low = candidate + 1;
                    }
                    else
                    {
                        q_high = candidate - 1;
                    }
                }

                if (q == 0)
                {
                    continue;
                }

                vector<DataT> q_vec = {q};
                vector<DataT> product = Multiply(shifted_divisor, q_vec, base);
                remainder = Subtract(remainder, product, base);

                if (quotient.size() <= static_cast<size_t>(d))
                {
                    quotient.resize(d + 1, 0);
                }
                quotient[d] = q;
            }

            TrimZeros(quotient);
            if (quotient.empty())
            {
                quotient.push_back(0);
            }

            TrimZeros(remainder);

            return {quotient, remainder};
        }
    };
}

#endif
