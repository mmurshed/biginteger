#ifndef BURNIKEL_ZIEGLER_DIVISION
#define BURNIKEL_ZIEGLER_DIVISION

#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <utility>
using namespace std;

#include "../../BigInteger.h"
#include "../../BigIntegerUtil.h"
#include "../classic/CommonAlgorithms.h"
#include "../classic/ClassicDivision.h"
#include "../../ops/BigIntegerAddition.h"
#include "../../ops/BigIntegerShift.h"
#include "../../ops/BigIntegerComparison.h"
#include "../../ops/BigIntegerIO.h"

namespace BigMath
{
    class BurnikelZieglerDivision
    {
        static const SizeT DIVISION_THRESHOLD = 4;

    private:
        static int Normalize1(vector<DataT> &a, vector<DataT> &b)
        {
            const Int limbBits = sizeof(DataT) * CHAR_BIT;
            vector<DataT> dvec = b;
            DataT msd = dvec.back();

            Int s = 0;
            while (((ULong)msd << s) < (static_cast<ULong>(1) << (limbBits - 1)))
                s++;

            a = CommonAlgorithms::ShiftLeftBits(a, s);
            b = CommonAlgorithms::ShiftLeftBits(b, s);
            return s;
        }

        // Corrected normalization with proper limb shifting
        static int Normalize(vector<DataT> &a, vector<DataT> &b)
        {
            const int limbBits = sizeof(DataT) * CHAR_BIT;
            DataT &msd = b.back();
            int s = 0;

            // Calculate required shift to set MSB of MSD
            while ((msd & (1 << (limbBits - 1))) == 0)
            {
                DataT carry = 0;

                // Shift left entire divisor b
                for (auto &limb : b)
                {
                    DataT new_carry = limb >> (limbBits - 1);
                    limb = (limb << 1) | carry;
                    carry = new_carry;
                }
                if (carry != 0)
                    b.push_back(carry);

                // Shift left dividend a similarly
                carry = 0;
                for (auto &limb : a)
                {
                    DataT new_carry = limb >> (limbBits - 1);
                    limb = (limb << 1) | carry;
                    carry = new_carry;
                }
                if (carry != 0)
                    a.push_back(carry);

                s++;
                msd = b.back(); // Update MSD after shift
            }
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
            vector<DataT> combined_digits;
            combined_digits.reserve(low.size() + high.size());

            // Append low-order digits first (low block must remain at its fixed length).
            combined_digits.insert(combined_digits.end(), low.begin(), low.end());
            // Append high-order digits (more significant) after low digits.
            combined_digits.insert(combined_digits.end(), high.begin(), high.end());

            // DO NOT trim trailing zeros—these zeros are significant for the fixed block size.
            // If you need to debug or print a "clean" version for display, do that separately.

            return combined_digits;
        }

        static pair<vector<DataT>, vector<DataT>> SplitAt(const vector<DataT> &num, SizeT n)
        {
            SizeT size = num.size();
            std::cerr << "[SplitAt] num: " << num << ", split index n: " << n << "\n";

            if (n >= size)
            {
                std::cerr << "[SplitAt] Split index >= size → returning {0, num}\n";
                return {vector<DataT>{0}, num};
            }

            vector<DataT> high(num.begin() + n, num.end());
            vector<DataT> low(num.begin(), num.begin() + n);

            std::cerr << "[SplitAt] high: " << high << ", low: " << low << "\n";
            return {high, low};
        }

        static pair<vector<DataT>, vector<DataT>> Divide2n1n(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            std::cerr << "[Divide2n1n] a: " << a << ", b: " << b << "\n";

            if (BigIntegerComparator::Compare(a, b) < 0)
            {
                std::cerr << "[Divide2n1n] Base case: " << a << " < " << b << endl;
                return {vector<DataT>{0}, a};
            }

            if (b.size() <= DIVISION_THRESHOLD)
            {
                std::cerr << "[Divide2n1n] Base case: b.size() <= threshold\n";
                auto result = ClassicDivision::DivideAndRemainder(a, b, base);
                std::cerr << "[Divide2n1n] Base case result: " << result.first << ", " << result.second << "\n";
                return result;
            }

            SizeT n = b.size() / 2;
            if (n == 0 || n >= b.size() || a.size() < n)
            {
                std::cerr << "[Divide2n1n] Base case: invalid n or small a\n";
                auto result = ClassicDivision::DivideAndRemainder(a, b, base);
                std::cerr << "[Divide2n1n] Base case result: " << result.first << ", " << result.second << "\n";
                return result;
            }

            auto [b1, b0] = SplitAt(b, n);
            auto [a12, a3] = SplitAt(a, n);

            if (BigIntegerUtil::IsZero(a12))
            {
                std::cerr << "[Divide2n1n] Base case: a12 is zero\n";
                auto result = ClassicDivision::DivideAndRemainder(a, b, base);
                std::cerr << "[Divide2n1n] Base case result: " << result.first << ", " << result.second << "\n";
                return result;
            }

            auto [a1, a2] = SplitAt(a12, n);

            std::cerr << "[Divide2n1n] Recursing on q1 with Combine(a1,a2), b1\n";
            auto [q1, r1] = Divide3n2n(Combine(a1, a2), b1, base);
            auto r1a3 = Combine(r1, a3);
            std::cerr << "[Divide2n1n] Recursing on q0 with r1a3, b\n";
            auto [q0, r] = Divide3n2n(r1a3, b, base);

            auto q = Combine(q1, q0);

            return {q, r};
        }

        static pair<vector<DataT>, vector<DataT>> Divide3n2n(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            std::cerr << "[Divide3n2n] a: " << a << ", b: " << b << "\n";

            if (BigIntegerComparator::Compare(a, b) < 0)
            {
                std::cerr << "[Divide3n2n] Base case: " << a << " < " << b << endl;
                return {vector<DataT>{0}, a};
            }

            if (b.size() <= DIVISION_THRESHOLD)
            {
                std::cerr << "[Divide3n2n] Base case: b.size() <= threshold\n";
                auto result = ClassicDivision::DivideAndRemainder(a, b, base);
                std::cerr << "[Divide3n2n] Base case result: " << result.first << ", " << result.second << "\n";
                return result;
            }

            SizeT n = b.size() / 2;
            if (n == 0 || n >= b.size() || a.size() < 2 * n)
            {
                std::cerr << "[Divide3n2n] Base case: invalid n or small a\n";
                auto result = ClassicDivision::DivideAndRemainder(a, b, base);
                std::cerr << "[Divide3n2n] Base case result: " << result.first << ", " << result.second << "\n";
                return result;
            }

            auto [b1, b0] = SplitAt(b, n);
            auto [a12, a3] = SplitAt(a, 2 * n);

            if (BigIntegerUtil::IsZero(a12))
            {
                std::cerr << "[Divide3n2n] Base case: a12 is zero\n";
                auto result = ClassicDivision::DivideAndRemainder(a, b, base);
                std::cerr << "[Divide2n1n] Base case result: " << result.first << ", " << result.second << "\n";
                return result;
            }

            std::cerr << "[Divide3n2n] Recursing on q1 with a12, b\n";
            auto [q1, r1] = Divide2n1n(a12, b, base);
            auto r1a3 = Combine(r1, a3);
            std::cerr << "[Divide3n2n] Recursing on q0 with r1a3, b\n";
            auto [q0, r] = Divide2n1n(r1a3, b, base);

            vector<DataT> q = Combine(q1, q0);

            return {q, r};
        }

    public:
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT> &a, const vector<DataT> &b, BaseT base)
        {
            std::cerr << "[DivideAndRemainder] ENTRY: a = " << a << ", b = " << b << "\n";

            if (BigIntegerUtil::IsZero(b))
                throw runtime_error("Division by zero");

            if (BigIntegerComparator::Compare(a, b) < 0)
            {
                std::cerr << "[DivideAndRemainder] Base case: " << a << " < " << b << endl;
                return {vector<DataT>{0}, a};
            }

            vector<DataT> a_norm = a;
            vector<DataT> b_norm = b;
            int shift = Normalize(a_norm, b_norm);

            std::cerr << "[DivideAndRemainder] Normalized: a = " << a_norm << ", b = " << b_norm << "\n";

            auto [q, r] = Divide3n2n(a_norm, b_norm, base);
            auto r_denorm = Denormalize(r, shift);

            BigIntegerUtil::TrimZeros(q);
            BigIntegerUtil::TrimZeros(r_denorm);
            std::cerr << "[DivideAndRemainder] DONE. q = " << q << ", r = " << r_denorm << "\n";
            return {q, r_denorm};
        }
    };
} // namespace BigMath

#endif // BURNIKEL_ZIEGLER_DIVISION