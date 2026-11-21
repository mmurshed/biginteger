// Toom threshold sweep for the Karatsuba recursion (64-bit limbs).
//
// Header-only: it does not link libbigmath, so each build can set its own
// thresholds without mixing definitions with the library's. Build one binary
// per value and compare the rows:
//
//   for T in 1000000 128 168 200; do
//     c++ -std=c++20 -O3 -march=native -DNDEBUG -Iinclude \
//         -DBIGMATH_MUL_TOOM33_THRESHOLD=$T -DBIGMATH_MUL_TOOM44_THRESHOLD=1000000 \
//         tests/performance/toom_threshold_sweep.cpp -o toom_$T && ./toom_$T
//   done
//
// "check" as the first argument verifies products of 2,664 shapes against
// the schoolbook kernel instead of timing. 1000000 disables a split.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>
#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
using namespace BigMath;
static volatile uint64_t sink;
int main(int argc, char **argv)
{
    std::string mode = argc > 1 ? argv[1] : "time";
    std::mt19937_64 gen(9);
    auto rnd = [&](size_t n, int kind) { std::vector<DataT> v(n); for (auto &x : v) x = kind == 1 ? ~(DataT)0 : (kind == 2 ? (gen() & 1 ? ~(DataT)0 : 0) : gen()); v.back() |= 1; return v; };
    if (mode == "check")
    {
        long bad = 0, cases = 0;
        for (size_t la = 25; la <= 700; la += (la < 120 ? 1 : 37))
            for (size_t lb : {la, la - la / 7, la - la / 3, (la * 2) / 3 + 1, la / 2 + 1, la + la / 5, la * 2 + 3, la * 5})
                for (int kind = 0; kind < 3; ++kind)
                {
                    if (lb < 2) continue;
                    auto a = rnd(la, kind), b = rnd(lb, (kind + 1) % 3);
                    auto c = KaratsubaMultiplication::Multiply(a, b, Base2_64);
                    std::vector<DataT> ref(la + lb);
                    LimbKernel::MulBasecase(ref.data(), a.data(), la, b.data(), lb);
                    while (ref.size() > 1 && ref.back() == 0) ref.pop_back();
                    ++cases;
                    if (c != ref) { if (++bad < 5) printf("MISMATCH %zu x %zu kind %d\n", la, lb, kind); }
                }
        printf("check T=%d: %ld cases, %ld mismatches\n", BIGMATH_MUL_TOOM33_THRESHOLD, cases, bad);
        return bad != 0;
    }
    printf("T3=%7d T4=%7d", BIGMATH_MUL_TOOM33_THRESHOLD, BIGMATH_MUL_TOOM44_THRESHOLD);
    for (size_t n : {128, 160, 192, 224, 256, 320, 384, 448, 512, 640, 768, 1024})
    {
        auto a = rnd(n, 0), b = rnd(n, 0);
        long reps = std::max(20L, 40000000L / (long)(n * n));
        double best = 1e18;
        for (int r = 0; r < 9; ++r) { auto t0 = std::chrono::steady_clock::now(); for (long i = 0; i < reps; ++i) { auto c = KaratsubaMultiplication::Multiply(a, b, Base2_64); sink += c[0]; } best = std::min(best, std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() / reps); }
        printf(" %8.2f", best);
    }
    printf("\n");
}
