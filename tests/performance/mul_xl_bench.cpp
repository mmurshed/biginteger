// XL multiplication bench — measures NTT/CRT path at 1M-50M digits.
//
// Build with the cmake `bigmath_add_test` macro or:
//   c++ -std=c++20 -O3 -march=native -Iinclude \
//       tests/performance/mul_xl_bench.cpp src/**/*.cpp bigdecimal/**/*.cpp \
//       -o mul_xl_bench

#include <chrono>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/ops/Operations.h"
#include "biginteger/algorithms/multiplication/NTTMultiplication.h"

using namespace std;
using namespace BigMath;
using Clock = chrono::high_resolution_clock;

static vector<DataT> RandomLimbs(size_t n, uint64_t seed)
{
    mt19937_64 gen(seed);
    uniform_int_distribution<uint64_t> dis;
    vector<DataT> v(n);
    for (size_t i = 0; i < n; ++i)
        v[i] = (DataT)dis(gen);
    if (v.back() == 0) v.back() = 1;
    return v;
}

template <typename F>
static double TimeMs(F &&fn, int iters)
{
    double best = 1e300;
    for (int i = 0; i < iters; ++i)
    {
        auto t0 = Clock::now();
        fn();
        auto t1 = Clock::now();
        double ms = chrono::duration<double, milli>(t1 - t0).count();
        if (ms < best) best = ms;
    }
    return best;
}

int main(int argc, char **argv)
{
    BaseT base = BigInteger::Base();
    printf("NTT XL bench  base=%s\n", base == Base2_64 ? "Base2_64" : "Base2_32");
    printf("%-12s %14s %12s\n", "limbs", "ms (best)", "iters");
    printf("%-12s %14s %12s\n", "-----", "---------", "-----");

    // Limb counts roughly mapping to digit decimals via *9.6
    // 1M, 2M, 5M, 10M, 20M, 50M limbs
    vector<size_t> sizes;
    if (argc > 1)
    {
        for (int k = 1; k < argc; ++k)
            sizes.push_back((size_t)atoll(argv[k]));
    }
    else
    {
        sizes = {100'000, 500'000, 1'000'000, 2'000'000, 5'000'000};
    }

    for (size_t n : sizes)
    {
        auto a = RandomLimbs(n, 0xA1B2C3D4ULL ^ n);
        auto b = RandomLimbs(n, 0xE5F6A7B8ULL ^ n);

        int iters = n <= 100'000 ? 5 : (n <= 1'000'000 ? 3 : 1);
        double ms = TimeMs([&]() {
            auto c = NTTMultiplication::Multiply(a, b, base);
            if (c.empty()) abort();
        }, iters);
        printf("%-12zu %14.3f %12d\n", n, ms, iters);
        fflush(stdout);
    }
    return 0;
}
