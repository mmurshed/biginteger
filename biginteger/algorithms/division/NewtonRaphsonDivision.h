#ifndef NEWTONRAPHSON_DIVISION
#define NEWTONRAPHSON_DIVISION

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "../../common/Builder.h"
#include "../../common/Util.h"
#include "../../ops/IO.h"
#include "../Multiplication.h"
#include "../Shift.h"
#include "../Addition.h"
#include "../Subtraction.h"
#include "KnuthDivision.h"

using namespace std;

namespace BigMath
{
    class NewtonRaphsonDivision
    {
    public:
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            if (IsZero(b))
            {
                throw invalid_argument("Division by zero");
            }

            int m = b.size();

            int extra_digits = m;
            int scaling_digits = 2 * m + extra_digits;
            vector<DataT> ten_power_F(scaling_digits, 0);
            ten_power_F.push_back(1); // Represents 10^scaling_digits

            vector<DataT> X;
            const int MAX_FLT_DIGITS = 15; // Max digits for accurate 64-bit float

            // --- Optimized Initial Approximation for Small Divisors ---
            if (m <= MAX_FLT_DIGITS)
            {
                // Extract B_high as the full trimmed divisor (little-endian to big-endian)
                vector<DataT> b_reversed(b.rbegin(), b.rend());
                TrimZeros(b_reversed); // Remove leading zeros (now in big-endian)

                // Convert b_reversed to a 64-bit integer
                uint64_t B_high = 0;
                for (auto d : b_reversed)
                {
                    B_high = B_high * 10 + d;
                }

                // Compute X0_high = 10^(scaling_digits) / B_high using floating-point
                double scale = pow(10.0, scaling_digits);
                double X0_high_dbl = scale / B_high;
                uint64_t X0_high = static_cast<uint64_t>(X0_high_dbl);

                // Convert X0_high to little-endian vector<DataT>
                vector<DataT> X0_vec;
                uint64_t temp = X0_high;
                do
                {
                    X0_vec.push_back(temp % 10);
                    temp /= 10;
                } while (temp > 0);
                X = X0_vec;
            }
            else
            {
                // Fallback to KnuthDivision for large divisors
                auto div_result = KnuthDivision::DivideAndRemainder(ten_power_F, b, base);
                X = div_result.first;
            }

            // --- Newton-Raphson Iterations (Optimized) ---
            int iterations = 3; // Reduced due to better initial guess
            for (int iter = 0; iter < iterations; ++iter)
            {
                vector<DataT> product = Multiply(b, X, base);
                vector<DataT> two_times_ten_power = Multiply(ten_power_F, 2, base);

                if (Compare(two_times_ten_power, product) < 0)
                {
                    throw runtime_error("Subtraction underflow");
                }
                vector<DataT> temp = Subtract(two_times_ten_power, product, base);
                vector<DataT> delta = Multiply(X, temp, base);

                // Truncate instead of dividing by ten_power_F
                if (delta.size() >= scaling_digits)
                {
                    X = vector<DataT>(delta.begin() + scaling_digits, delta.end());
                }
                else
                {
                    X = {0};
                }
                TrimZeros(X);
            }

            // --- Compute Quotient and Remainder ---
            vector<DataT> a_times_X = Multiply(a, X, base);
            vector<DataT> quotient;
            if (a_times_X.size() >= scaling_digits)
            {
                quotient = vector<DataT>(a_times_X.begin() + scaling_digits, a_times_X.end());
            }
            else
            {
                quotient = {0};
            }
            TrimZeros(quotient);

            // --- Final Correction ---
            vector<DataT> q_times_b = Multiply(quotient, b, base);
            vector<DataT> diff = Subtract(a, q_times_b, base);
            if (Compare(diff, b) >= 0)
            {
                auto extra_div = KnuthDivision::DivideAndRemainder(diff, b, base);
                AddTo(quotient, extra_div.first, base);
                diff = extra_div.second;
            }

            TrimZeros(quotient);
            TrimZeros(diff);
            return {quotient, diff};
        }
    };
}
#endif