#ifndef BURNIKEL_ZIEGLER_DIVISION
#define BURNIKEL_ZIEGLER_DIVISION

#include <stdexcept>
#include <utility>
using namespace std;

#include "../../BigInteger.h"
#include "../BigIntegerUtil.h"
#include "../classic/CommonAlgorithms.h"
#include "../classic/ClassicDivision.h"
#include "../../ops/BigIntegerAddition.h"
#include "../../ops/BigIntegerShift.h"
#include "../../ops/BigIntegerComparison.h"

namespace BigMath
{
    class BurnikelZieglerDivision
    {
        static const SizeT DIVISION_THRESHOLD = 4;

    private:
        static int Normalize(vector<DataT> &a, vector<DataT> &b)
        {
            const Int limbBits = sizeof(DataT) * CHAR_BIT;
            vector<DataT> dvec = b;
            DataT msd = dvec.back();

            Int s = 0;
            while (((UInt)msd << s) < (1U << (limbBits - 1)))
                s++;

            a = CommonAlgorithms::ShiftLeftBits(a, s);
            b = CommonAlgorithms::ShiftLeftBits(b, s);
            return s;
        }

        static vector<DataT> Denormalize(const vector<DataT> &r, int s)
        {
            return CommonAlgorithms::ShiftRightBits(r, s);
        }

        /// Combines two normalized BigIntegers into a single number where:
        /// combined = (low + high * base^num_low_digits)
        static vector<DataT> Combine(const vector<DataT> &high, const vector<DataT> &low)
        {
            const auto &low_digits = low;
            const auto &high_digits = high;

            vector<DataT> combined_digits;
            combined_digits.reserve(low_digits.size() + high_digits.size());

            // Append low-order digits first (little-endian: least significant first)
            combined_digits.insert(combined_digits.end(),
                                   low_digits.begin(),
                                   low_digits.end());

            // Append high-order digits (more significant) after low digits
            combined_digits.insert(combined_digits.end(),
                                   high_digits.begin(),
                                   high_digits.end());

            // Remove trailing zeros (which are leading zeros in big-endian representation)
            while (!combined_digits.empty() && combined_digits.back() == 0)
            {
                combined_digits.pop_back();
            }

            // Handle all-zero case
            if (combined_digits.empty())
            {
                return vector<DataT>{0};
            }

            return std::move(combined_digits);
        }

        static pair<vector<DataT>, vector<DataT>> SplitAt(const vector<DataT> &num, SizeT n)
        {
            vector<DataT> digits = num;
            if (n >= digits.size())
                return {vector<DataT>{0}, num};

            return {
                vector<DataT>(digits.begin() + n, digits.end()),
                vector<DataT>(digits.begin(), digits.begin() + n)};
        }

        static pair<vector<DataT>, vector<DataT>> Divide2n1n(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            // If a < b, then quotient is 0
            if (BigIntegerComparator::Compare(a, b) < 0)
                return {vector<DataT>{0}, a};

            if (b.size() <= DIVISION_THRESHOLD)
                return ClassicDivision::DivideAndRemainder(a, b, base);

            SizeT n = b.size() / 2;
            // If n is 0 or n is not a proper split or a is too small, fallback
            if (n == 0 || n >= b.size() || a.size() < n)
                return ClassicDivision::DivideAndRemainder(a, b, base);

            auto [b1, b0] = SplitAt(b, n);

            auto [a12, a3] = SplitAt(a, n);

            // If a12 is trivial (i.e. zero), use classical division.
            if (BigIntegerUtil::IsZero(a12))
                return ClassicDivision::DivideAndRemainder(a, b, base);

            auto [a1, a2] = SplitAt(a12, n);

            auto [q1, r1] = Divide3n2n(Combine(a1, a2), b1, base);
            auto r1a3 = Combine(r1, a3);
            auto [q0, r] = Divide3n2n(r1a3, b, base);

            return {Combine(q1, q0), r};
        }

        static pair<vector<DataT>, vector<DataT>> Divide3n2n(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            // If a < b, quotient is 0.
            if (BigIntegerComparator::Compare(a, b) < 0)
                return {vector<DataT>{0}, a};

            if (b.size() <= DIVISION_THRESHOLD)
                return ClassicDivision::DivideAndRemainder(a, b, base);

            SizeT n = b.size() / 2;
            if (n == 0 || n >= b.size() || a.size() < 2 * n)
                return ClassicDivision::DivideAndRemainder(a, b, base);

            auto [b1, b0] = SplitAt(b, n);
            auto [a12, a3] = SplitAt(a, 2 * n);
            if (BigIntegerUtil::IsZero(a12))
                return ClassicDivision::DivideAndRemainder(a, b, base);

            auto [q1, r1] = Divide2n1n(a12, b, base);
            auto r1a3 = Combine(r1, a3);
            auto [q0, r] = Divide2n1n(r1a3, b, base);

            vector<DataT> q = Combine(q1, q0);
            return {q, r};
        }

    public:
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            if (BigIntegerUtil::IsZero(b))
                throw runtime_error("Division by zero");
            if (BigIntegerComparator::Compare(a, b) < 0)
                return {vector<DataT>{0}, a};

            vector<DataT> a_norm = a;
            vector<DataT> b_norm = b;
            int shift = Normalize(a_norm, b_norm);

            auto [q, r] = Divide3n2n(a_norm, b_norm, base);
            return {q, Denormalize(r, shift)};
        }
    };
} // namespace BigMath

#endif // BURNIKEL_ZIEGLER_DIVISION