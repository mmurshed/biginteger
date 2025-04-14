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
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT> &a, const vector<DataT> &b, DataT base)
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

            while (Compare(remainder, b) >= 0)
            {
                int d = remainder.size() - b.size();
                vector<DataT> shifted_divisor;

                while (d >= 0)
                {
                    shifted_divisor.clear();
                    shifted_divisor.insert(shifted_divisor.end(), d, 0);
                    shifted_divisor.insert(shifted_divisor.end(), b.begin(), b.end());
                    if (Compare(shifted_divisor, remainder) <= 0)
                    {
                        break;
                    }
                    d--;
                }

                if (d < 0)
                {
                    break;
                }

                DataT q = 0;
                for (DataT candidate = base - 1; candidate >= 1; --candidate)
                {
                    vector<DataT> candidate_vec = {candidate};
                    vector<DataT> product = Multiply(shifted_divisor, candidate_vec, base);
                    if (Compare(product, remainder) <= 0)
                    {
                        q = candidate;
                        break;
                    }
                }

                if (q == 0)
                {
                    throw runtime_error("Failed to find quotient digit");
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

            // Convert quotient to little-endian with proper trimming
            vector<DataT> little_endian_quotient;
            for (size_t i = 0; i < quotient.size(); ++i)
            {
                if (i < quotient.size() && quotient[i] != 0)
                {
                    little_endian_quotient.clear();
                    little_endian_quotient.insert(little_endian_quotient.end(), quotient.begin() + i, quotient.end());
                    break;
                }
            }
            if (little_endian_quotient.empty())
            {
                little_endian_quotient.push_back(0);
            }

            reverse(little_endian_quotient.begin(), little_endian_quotient.end());
            TrimZeros(little_endian_quotient);

            TrimZeros(remainder);

            return {little_endian_quotient, remainder};
        }
    };
}

#endif
