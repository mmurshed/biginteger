#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <random>

#include "biginteger/BigInteger.h"
#include "biginteger/ops/IO.h"
#include "biginteger/ops/Comparison.h"
#include "biginteger/common/Parser.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
#include "biginteger/algorithms/multiplication/NTTMultiplication.h"

using namespace std;
using namespace BigMath;

// Helper to generate a random number string of given digit length
string GenerateRandomDigits(int length)
{
    static mt19937 gen(42); // Seeded for reproducibility
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

void RunBenchmarkForSize(int digits)
{
    cout << "---------------------------------------------" << endl;
    cout << "Benchmarking with two " << digits << "-digit random numbers" << endl;
    cout << "---------------------------------------------" << endl;

    string s_a = GenerateRandomDigits(digits);
    string s_b = GenerateRandomDigits(digits);

    BigInteger a = Parse(s_a.c_str());
    BigInteger b = Parse(s_b.c_str());

    const auto &vecA = a.GetInteger();
    const auto &vecB = b.GetInteger();

    // 1. Classical Multiplication
    auto start = chrono::high_resolution_clock::now();
    vector<DataT> rClassic = ClassicMultiplication::Multiply(vecA, vecB, BigInteger::Base());
    auto end = chrono::high_resolution_clock::now();
    double timeClassic = chrono::duration<double, milli>(end - start).count();

    // 2. Karatsuba Multiplication
    start = chrono::high_resolution_clock::now();
    vector<DataT> rKara = KaratsubaMultiplication::Multiply(vecA, vecB, BigInteger::Base());
    end = chrono::high_resolution_clock::now();
    double timeKara = chrono::duration<double, milli>(end - start).count();

    // 3. NTT Multiplication
    start = chrono::high_resolution_clock::now();
    vector<DataT> rNTT = NTTMultiplication::Multiply(vecA, vecB, BigInteger::Base());
    end = chrono::high_resolution_clock::now();
    double timeNTT = chrono::duration<double, milli>(end - start).count();

    // Validate correctness
    bool cmp_kara_ok = (Compare(rClassic, rKara) == 0);
    bool cmp_ntt_ok = (Compare(rClassic, rNTT) == 0);

    if (!cmp_kara_ok || !cmp_ntt_ok)
    {
        cerr << "ERROR: Algorithmic results differ!" << endl;
        cerr << "Classic limbs: " << rClassic.size() << endl;
        cerr << "Karatsuba: " << (cmp_kara_ok ? "Match" : "FAILED")
             << " limbs=" << rKara.size() << endl;
        cerr << "NTT:       " << (cmp_ntt_ok ? "Match" : "FAILED")
             << " limbs=" << rNTT.size() << endl;
        exit(1);
    }

    cout << "Validation: SUCCESS (All products match exactly)" << endl;
    printf("Classical:  %8.3f ms\n", timeClassic);
    printf("Karatsuba:  %8.3f ms (Speedup vs Classic: %4.1fx)\n", timeKara, timeClassic / timeKara);
    printf("NTT:        %8.3f ms (Speedup vs Classic: %4.1fx)\n", timeNTT, timeClassic / timeNTT);
}

int main()
{
    cout << "BigInteger Multiplication Performance Benchmark" << endl;
    cout << "Internal representation: Base 2^32 (" << BigInteger::Base() << ")" << endl;
    cout << "---------------------------------------------" << endl;

    RunBenchmarkForSize(100);
    RunBenchmarkForSize(10000);
    RunBenchmarkForSize(100000);
    RunBenchmarkForSize(500000);

    return 0;
}
