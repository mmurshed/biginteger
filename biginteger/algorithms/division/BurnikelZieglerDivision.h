#ifndef BURNIKELZIEGLER_DIVISION
#define BURNIKELZIEGLER_DIVISION

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
    class BurnikelZieglerDivision
    {
        static vector<DataT> MultiplyByBasePower(const vector<DataT>& num, SizeT power, BaseT base) {
            vector<DataT> result = num;
            result.insert(result.begin(), power, 0);
            return result;
        }

        public:
        
        static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(const vector<DataT>& a, const vector<DataT>& b, BaseT base) {      
            if (IsZero(b)) {
                throw invalid_argument("Division by zero");
            }
        
            int cmp = Compare(a, b);
            if (cmp < 0) {
                return {vector<DataT>(), a};
            } else if (cmp == 0) {
                return {vector<DataT>{1}, vector<DataT>()};
            }
        
            SizeT m = b.size();
            SizeT n = a.size();
        
            // Base case: use Knuth's algorithm for small divisors or when dividend is small
            if (m == 1 || n <= 2 * m) {
                return KnuthDivision::DivideAndRemainder(a, b, base);
            }
        
            // Split the dividend into high and low parts
            SizeT splitPos = 2 * m;
            if (splitPos > n) splitPos = n;
        
            vector<DataT> high(a.end() - splitPos, a.end());
            vector<DataT> low(a.begin(), a.end() - splitPos);
        
            // Recursively divide the high part by the divisor
            auto [qHigh, rHigh] = DivideAndRemainder(high, b, base);
        
            // Combine remainder with low part and divide
            vector<DataT> combined = MultiplyByBasePower(rHigh, low.size(), base);
            combined = Add(combined, low, base);
            auto [qLow, rLow] = DivideAndRemainder(combined, b, base);
        
            // Combine the quotients
            vector<DataT> qHighShifted = MultiplyByBasePower(qHigh, low.size(), base);
            vector<DataT> quotient = Add(qHighShifted, qLow, base);
        
            // Trim leading zeros (trailing in little-endian)
            TrimZeros(quotient);
            TrimZeros(rLow);
        
            return {quotient, rLow};
        }        
    };
}
#endif