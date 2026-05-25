// Benchmark BigMath vs GMP (mpz) on multiplication, division, decimal parse, ToString.
//
// Build:
//   c++ -std=c++20 -O2 -I/opt/homebrew/include -L/opt/homebrew/lib \
//       tests/performance/bench_vs_gmp.cpp -o /private/tmp/bench_vs_gmp -lgmp \
//       && /private/tmp/bench_vs_gmp
//
// Linux: drop the -I/-L for /opt/homebrew; ensure libgmp-dev installed.

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <gmp.h>

#include "../../biginteger/BigInteger.h"
#include "../../biginteger/ops/Operations.h"
#include "../../biginteger/common/Builder.h"
#include "../../biginteger/common/Parser.h"

using namespace std;
using namespace BigMath;

using Clock = chrono::high_resolution_clock;

static string GenerateRandomDigits(int length, uint64_t seed)
{
    mt19937_64 gen(seed);
    uniform_int_distribution<int> dis(0, 9);
    uniform_int_distribution<int> first_dis(1, 9);

    string s;
    s.reserve(length);
    s += char('0' + first_dis(gen));
    for (int i = 1; i < length; ++i)
        s += char('0' + dis(gen));
    return s;
}

// Run `fn` once if `single_shot`, otherwise repeat `iters` times and report the minimum.
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

static int IterationsForDigits(int digits)
{
    if (digits <= 2000) return 7;
    if (digits <= 20000) return 5;
    if (digits <= 100000) return 3;
    return 1;
}

struct Row
{
    string op;
    string size_label;
    double bm_ms;
    double gmp_ms;
};

static void PrintHeader()
{
    printf("%-10s %-22s %14s %14s %10s\n", "op", "size", "BigMath ms", "GMP ms", "BM/GMP");
    printf("%-10s %-22s %14s %14s %10s\n", "----", "----", "----------", "------", "------");
}

static void PrintRow(const Row &r)
{
    double ratio = (r.gmp_ms > 0.0) ? r.bm_ms / r.gmp_ms : 0.0;
    printf("%-10s %-22s %14.3f %14.3f %10.2fx\n",
           r.op.c_str(), r.size_label.c_str(), r.bm_ms, r.gmp_ms, ratio);
    fflush(stdout);
}

static void BenchMul(int digits_a, int digits_b, vector<Row> &rows)
{
    string sa = GenerateRandomDigits(digits_a, 0xA1ull ^ (uint64_t)digits_a);
    string sb = GenerateRandomDigits(digits_b, 0xB2ull ^ (uint64_t)digits_b);

    BigInteger a = Parse(sa.c_str());
    BigInteger b = Parse(sb.c_str());

    mpz_t ga, gb, gc;
    mpz_init_set_str(ga, sa.c_str(), 10);
    mpz_init_set_str(gb, sb.c_str(), 10);
    mpz_init(gc);

    int iters = IterationsForDigits(max(digits_a, digits_b));

    double bm = TimeMs([&]() {
        BigInteger c = a * b;
        if (c.size() == 0) abort();
    }, iters);

    double gp = TimeMs([&]() {
        mpz_mul(gc, ga, gb);
    }, iters);

    char buf[64];
    snprintf(buf, sizeof(buf), "%dx%d digits", digits_a, digits_b);
    Row row{"mul", buf, bm, gp};
    rows.push_back(row);
    PrintRow(row);

    mpz_clear(ga);
    mpz_clear(gb);
    mpz_clear(gc);
}

static void BenchDiv(int digits_a, int digits_b, vector<Row> &rows)
{
    string sa = GenerateRandomDigits(digits_a, 0xC3ull ^ (uint64_t)digits_a ^ ((uint64_t)digits_b << 16));
    string sb = GenerateRandomDigits(digits_b, 0xD4ull ^ (uint64_t)digits_b);

    BigInteger a = Parse(sa.c_str());
    BigInteger b = Parse(sb.c_str());

    mpz_t ga, gb, gq, gr;
    mpz_init_set_str(ga, sa.c_str(), 10);
    mpz_init_set_str(gb, sb.c_str(), 10);
    mpz_init(gq);
    mpz_init(gr);

    int iters = IterationsForDigits(max(digits_a, digits_b));

    double bm = TimeMs([&]() {
        BigInteger q = a / b;
        if (q.size() == 0) abort();
    }, iters);

    double gp = TimeMs([&]() {
        mpz_tdiv_qr(gq, gr, ga, gb);
    }, iters);

    char buf[64];
    snprintf(buf, sizeof(buf), "%dx%d digits", digits_a, digits_b);
    Row row{"div", buf, bm, gp};
    rows.push_back(row);
    PrintRow(row);

    mpz_clear(ga);
    mpz_clear(gb);
    mpz_clear(gq);
    mpz_clear(gr);
}

static void BenchParse(int digits, vector<Row> &rows)
{
    string s = GenerateRandomDigits(digits, 0xE5ull ^ (uint64_t)digits);
    const char *cs = s.c_str();

    int iters = IterationsForDigits(digits);

    double bm = TimeMs([&]() {
        BigInteger b = Parse(cs);
        if (b.size() == 0) abort();
    }, iters);

    double gp = TimeMs([&]() {
        mpz_t g;
        mpz_init(g);
        mpz_set_str(g, cs, 10);
        mpz_clear(g);
    }, iters);

    char buf[64];
    snprintf(buf, sizeof(buf), "%d digits", digits);
    Row row{"parse", buf, bm, gp};
    rows.push_back(row);
    PrintRow(row);
}

static void BenchToString(int digits, vector<Row> &rows)
{
    string s = GenerateRandomDigits(digits, 0xF6ull ^ (uint64_t)digits);
    BigInteger b = Parse(s.c_str());

    mpz_t g;
    mpz_init_set_str(g, s.c_str(), 10);

    int iters = IterationsForDigits(digits);
    // ToString is slow; cap iters harder for very large inputs.
    if (digits >= 100000) iters = 1;

    double bm = TimeMs([&]() {
        string out = ToString(b);
        if (out.empty()) abort();
    }, iters);

    double gp = TimeMs([&]() {
        char *out = mpz_get_str(nullptr, 10, g);
        if (!out) abort();
        free(out);
    }, iters);

    char buf[64];
    snprintf(buf, sizeof(buf), "%d digits", digits);
    Row row{"tostr", buf, bm, gp};
    rows.push_back(row);
    PrintRow(row);

    mpz_clear(g);
}

int main()
{
    cout << "BigMath vs GMP " << gmp_version << " (Apple Silicon / brew gmp)\n\n";
    PrintHeader();

    vector<Row> rows;

    // ===== Multiplication =====
    cout << "\n[multiplication, balanced]\n";
    PrintHeader();
    for (int d : {1000, 5000, 10000, 50000, 100000, 500000, 1000000})
        BenchMul(d, d, rows);

    cout << "\n[multiplication, skewed]\n";
    PrintHeader();
    BenchMul(100000, 10000, rows);
    BenchMul(500000, 50000, rows);
    BenchMul(1000000, 100000, rows);

    // ===== Division =====
    cout << "\n[division, balanced]\n";
    PrintHeader();
    for (int d : {1000, 5000, 10000, 50000, 100000, 500000})
        BenchDiv(d, d, rows);

    cout << "\n[division, skewed (a >> b) — Newton/BZ band]\n";
    PrintHeader();
    BenchDiv(40000, 10000, rows);
    BenchDiv(100000, 10000, rows);
    BenchDiv(200000, 50000, rows);
    BenchDiv(500000, 100000, rows);

    // ===== Decimal Parse =====
    cout << "\n[decimal parse]\n";
    PrintHeader();
    for (int d : {1000, 10000, 50000, 100000, 500000, 1000000})
        BenchParse(d, rows);

    // ===== ToString =====
    cout << "\n[ToString (known BigMath bottleneck)]\n";
    PrintHeader();
    for (int d : {1000, 10000, 50000, 100000})
        BenchToString(d, rows);

    cout << "\nDone. " << rows.size() << " measurements.\n";
    return 0;
}
