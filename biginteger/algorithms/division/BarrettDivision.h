#ifndef BARRETT_DIVISION
#define BARRETT_DIVISION

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "../../common/Builder.h"
#include "../../common/Util.h"
#include "../../ops/IO.h"
#include "../../algorithms/Multiplication.h"
#include "../Shift.h"
#include "../Addition.h"
#include "../Subtraction.h"
#include "../division/KnuthDivision.h"

using namespace std;

namespace BigMath
{
    class BarrettDivision
    {
    public:
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT> &a, const vector<DataT> &b, DataT base)
        {
            if (IsZero(b))
            {
                throw invalid_argument("Division by zero");
            }

            int m = b.size();

            // Compute base^(2m)
            vector<DataT> base_power_2m(2 * m, 0);
            base_power_2m.push_back(1);

            // Compute mu = floor(base^(2m) / b)
            auto mu_div = KnuthDivision::DivideAndRemainder(base_power_2m, b, base);
            vector<DataT> mu = mu_div.first;

            // Compute a * mu
            vector<DataT> a_times_mu = Multiply(a, mu, base);

            // Compute q_hat = floor(a_times_mu / base^(2m))
            size_t shift_amount = 2 * m;
            vector<DataT> q_hat;

            if (a_times_mu.size() >= shift_amount)
            {
                q_hat.assign(a_times_mu.begin() + shift_amount, a_times_mu.end());
            }
            else
            {
                q_hat.clear();
            }

            TrimZeros(q_hat);
            if (q_hat.empty())
            {
                q_hat.push_back(0);
            }

            // Adjust q_hat if it's larger than the actual quotient
            vector<DataT> product = Multiply(q_hat, b, base);
            int cmp = Compare(product, a);

            if (cmp > 0)
            {
                q_hat = Subtract(q_hat, {1}, base);
                product = Multiply(q_hat, b, base);
            }

            vector<DataT> remainder = Subtract(a, product, base);

            // Compute how many times b fits into remainder and adjust
            auto delta_qr = KnuthDivision::DivideAndRemainder(remainder, b, base);
            vector<DataT> delta_q = delta_qr.first;
            remainder = delta_qr.second;
            q_hat = Add(q_hat, delta_q, base);

            // Final trimming
            TrimZeros(q_hat);
            TrimZeros(remainder);

            return {q_hat, remainder};
        }
    };
}
#endif