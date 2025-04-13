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

            // Use extra guard digits for precision; here we choose extra_digits = m.
            int extra_digits = m;
            int scaling_digits = 2 * m + extra_digits;

            // Build the scaling factor: ten_power_F = 10^(scaling_digits) in little-endian representation.
            vector<DataT> ten_power_F(scaling_digits, 0);
            ten_power_F.push_back(1);

            // Compute the initial approximation X0 = ten_power_F / b.
            auto div_result = ClassicDivision::DivideAndRemainder(ten_power_F, b, base);
            vector<DataT> X = div_result.first;

            // Perform a fixed number of Newton-Raphson iterations to refine X.
            int iterations = 5; // Adjust for required precision.
            for (int iter = 0; iter < iterations; iter++)
            {
                vector<DataT> product = MultiplicationStrategy::MultiplyUnsigned(b, X, base);

                vector<DataT> two = {2};
                vector<DataT> two_times_ten_power = MultiplicationStrategy::MultiplyUnsigned(ten_power_F, two, base);

                if (BigIntegerComparator::Compare(two_times_ten_power, product) < 0)
                {
                    throw runtime_error("Subtraction underflow in Newton-Raphson iteration");
                }
                vector<DataT> temp = ClassicSubtraction::Subtract(two_times_ten_power, product, base);

                vector<DataT> delta = MultiplicationStrategy::MultiplyUnsigned(X, temp, base);

                auto delta_div = ClassicDivision::DivideAndRemainder(delta, ten_power_F, base);
                X = delta_div.first;
            }

            // Multiply dividend a by the refined reciprocal X.
            vector<DataT> a_times_X = MultiplicationStrategy::MultiplyUnsigned(a, X, base);

            // Compute the preliminary quotient: (a * X) / ten_power_F.
            auto div_q = ClassicDivision::DivideAndRemainder(a_times_X, ten_power_F, base);
            vector<DataT> quotient = div_q.first;

            // --- Robust Final Correction ---
            // Instead of incrementing one by one (or using an exponential/binary search),
            // compute diff = a - (quotient * b) and use full-precision division to determine
            // how many additional multiples of b to add.
            vector<DataT> one = {1};
            while (true)
            {
                vector<DataT> q_times_b = MultiplicationStrategy::MultiplyUnsigned(quotient, b, base);
                // Log the current product.
                vector<DataT> diff = ClassicSubtraction::Subtract(a, q_times_b, base);
                // If diff is less than b then the quotient is final.
                if (BigIntegerComparator::Compare(diff, b) < 0)
                    break;
                // Compute extra = floor(diff / b).
                auto extra_div = ClassicDivision::DivideAndRemainder(diff, b, base);
                vector<DataT> extra = extra_div.first;
                if (BigIntegerComparator::Compare(extra, one) < 0)
                    extra = one;
                quotient = ClassicAddition::Add(quotient, extra, base);
            }

            // Final remainder = a - (quotient * b).
            vector<DataT> final_product = MultiplicationStrategy::MultiplyUnsigned(quotient, b, base);
            vector<DataT> remainder = ClassicSubtraction::Subtract(a, final_product, base);

            BigIntegerUtil::TrimZeros(quotient);
            BigIntegerUtil::TrimZeros(remainder);
            return {quotient, remainder};
        }
    };
}
#endif