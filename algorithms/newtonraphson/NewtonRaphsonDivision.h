#ifndef NEWTONRAPHSON_DIVISION
#define NEWTONRAPHSON_DIVISION

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "../../BigIntegerBuilder.h"
#include "../../BigIntegerUtil.h"
#include "../../ops/BigIntegerIO.h"
#include "../../ops/MultiplicationStrategy.h"
#include "../classic/CommonAlgorithms.h"
#include "../classic/ClassicAddition.h"
#include "../classic/ClassicSubtraction.h"
#include "../classic/ClassicDivision.h"

using namespace std;

namespace BigMath
{
    class NewtonRaphsonDivision
    {
    public:
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            if (BigIntegerUtil::IsZero(b))
            {
                throw invalid_argument("Division by zero");
            }

            int m = b.size();

            // Use extra guard digits for precision. Here we choose extra_digits = m.
            int extra_digits = m;
            // Scaling factor exponent = 2*m + extra_digits.
            int scaling_digits = 2 * m + extra_digits;

            // Build scaling factor: ten_power_F = 10^(scaling_digits) in little-endian.
            vector<DataT> ten_power_F(scaling_digits, 0);
            ten_power_F.push_back(1); // The most significant digit.

            // Compute initial approximation X0 = ten_power_F / b.
            auto div_result = ClassicDivision::DivideAndRemainder(ten_power_F, b, base);
            vector<DataT> X = div_result.first;

            // Perform a fixed number of Newton-Raphson iterations to refine X.
            int iterations = 5; // Adjust as needed.
            for (int iter = 0; iter < iterations; iter++)
            {
                // product = b * X.
                vector<DataT> product = MultiplicationStrategy::MultiplyUnsigned(b, X, base);

                // two_times_ten_power = 2 * ten_power_F.
                vector<DataT> two = {2};
                vector<DataT> two_times_ten_power = MultiplicationStrategy::MultiplyUnsigned(ten_power_F, two, base);

                // Ensure subtraction is valid.
                if (BigIntegerComparator::Compare(two_times_ten_power, product) < 0)
                {
                    throw runtime_error("Subtraction underflow in Newton-Raphson iteration");
                }
                vector<DataT> temp = ClassicSubtraction::Subtract(two_times_ten_power, product, base);

                // delta = X * temp.
                vector<DataT> delta = MultiplicationStrategy::MultiplyUnsigned(X, temp, base);

                // Update X: X = delta / ten_power_F.
                auto delta_div = ClassicDivision::DivideAndRemainder(delta, ten_power_F, base);
                X = delta_div.first;
            }

            // Multiply a by refined X.
            vector<DataT> a_times_X = MultiplicationStrategy::MultiplyUnsigned(a, X, base);

            // Preliminary quotient = (a * X) / ten_power_F.
            auto div_q = ClassicDivision::DivideAndRemainder(a_times_X, ten_power_F, base);
            vector<DataT> quotient = div_q.first;

            // --- Robust Final Quotient Adjustment ---
            // First, if quotient is too high, decrement until (quotient * b) <= a.
            vector<DataT> q_times_b = MultiplicationStrategy::MultiplyUnsigned(quotient, b, base);
            while (BigIntegerComparator::Compare(q_times_b, a) > 0)
            {
                quotient = ClassicSubtraction::Subtract(quotient, vector<DataT>{1}, base);
                q_times_b = MultiplicationStrategy::MultiplyUnsigned(quotient, b, base);
            }

            // Now, adjust upward. Instead of incrementing one by one, compute a correction factor.
            while (true)
            {
                vector<DataT> quotient_plus_one = ClassicAddition::Add(quotient, vector<DataT>{1}, base);
                vector<DataT> q_plus_one_times_b = MultiplicationStrategy::MultiplyUnsigned(quotient_plus_one, b, base);

                // If (quotient + 1) * b > a then our quotient is finalized.
                if (BigIntegerComparator::Compare(q_plus_one_times_b, a) > 0)
                {
                    break;
                }

                // Compute the difference diff = a - (quotient * b).
                vector<DataT> diff = ClassicSubtraction::Subtract(a, q_times_b, base);
                // Compute how many multiples of b fit in diff.
                auto extra_div = ClassicDivision::DivideAndRemainder(diff, b, base);
                vector<DataT> extra = extra_div.first;
                // Ensure that at least one unit is added.
                if (BigIntegerComparator::Compare(extra, vector<DataT>{1}) < 0)
                {
                    extra = {1};
                }
                // Update quotient: quotient = quotient + extra.
                quotient = ClassicAddition::Add(quotient, extra, base);
                q_times_b = MultiplicationStrategy::MultiplyUnsigned(quotient, b, base);
            }

            // Final remainder = a - (quotient * b).
            vector<DataT> remainder = ClassicSubtraction::Subtract(a, q_times_b, base);

            BigIntegerUtil::TrimZeros(quotient);
            BigIntegerUtil::TrimZeros(remainder);
            return {quotient, remainder};
        }
    };
}
#endif