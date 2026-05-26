#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "biginteger/BigInteger.h"
#include "biginteger/ops/IO.h"
#include "biginteger/ops/Comparison.h"
#include "biginteger/common/Parser.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"

using namespace std;
using namespace BigMath;

string GenerateRandomDigits(int length)
{
    static mt19937 gen(42);
    uniform_int_distribution<> dis(0, 9);
    uniform_int_distribution<> first_dis(1, 9);

    string s = "";
    s += to_string(first_dis(gen));
    for (int i = 1; i < length; ++i)
    {
        s += to_string(dis(gen));
    }
    return s;
}

int main()
{
    vector<int> sizes = {100, 10000, 100000, 500000};
    for (int digits : sizes)
    {
        cout << "Testing size " << digits << "..." << endl;
        string s_a = GenerateRandomDigits(digits);
        string s_b = GenerateRandomDigits(digits);

        BigInteger a = Parse(s_a.c_str());
        BigInteger b = Parse(s_b.c_str());

        const auto &vecA = a.GetInteger();
        const auto &vecB = b.GetInteger();

        vector<DataT> rClassic = ClassicMultiplication::Multiply(vecA, vecB, BigInteger::Base());
        vector<DataT> rKara = KaratsubaMultiplication::Multiply(vecA, vecB, BigInteger::Base());

        if (Compare(rClassic, rKara) != 0)
        {
            cout << "FAIL at size " << digits << "!" << endl;
            cout << "rClassic size: " << rClassic.size() << ", rKara size: " << rKara.size() << endl;

            // Find the first index where they differ
            size_t limit = min(rClassic.size(), rKara.size());
            for (size_t i = 0; i < limit; ++i)
            {
                if (rClassic[i] != rKara[i])
                {
                    cout << "First difference at limb index " << i << ":" << endl;
                    cout << "  rClassic[" << i << "] = " << rClassic[i] << endl;
                    cout << "  rKara[" << i << "]    = " << rKara[i] << endl;
                    cout << "  Difference = " << (long long)rClassic[i] - (long long)rKara[i] << endl;
                    
                    // Print surrounding limbs
                    for (int k = max(0, (int)i - 2); k <= min((int)limit - 1, (int)i + 2); ++k)
                    {
                        cout << "  Limb " << k << " -> Classic: " << rClassic[k] << ", Kara: " << rKara[k] << endl;
                    }
                    break;
                }
            }
            return 1;
        }
    }

    cout << "SUCCESS: All sizes matched!" << endl;
    return 0;
}
