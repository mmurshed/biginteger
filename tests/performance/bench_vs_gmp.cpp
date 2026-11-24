// Benchmark BigMath vs GMP (mpz) on multiplication, division, decimal parse, ToString.
//
// Build with cmake --build build --target bench_vs_gmp (requires GMP).

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include <gmp.h>

#include "biginteger/BigInteger.h"
#include "biginteger/ops/Operations.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "bigdecimal/BigDecimal.h"

using namespace std;
using namespace BigMath;

using Clock = chrono::steady_clock;

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

// Report the fastest of the requested iterations.
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
    int reps;
    const char *phase = "best";
};

static void PrintHeader()
{
    printf("%-10s %-22s %14s %14s %10s %6s %8s\n", "op", "size", "BigMath ms", "GMP ms", "BM/GMP", "reps", "phase");
    printf("%-10s %-22s %14s %14s %10s %6s %8s\n", "----", "----", "----------", "------", "------", "----", "-----");
}

static void PrintRow(const Row &r)
{
    double ratio = (r.gmp_ms > 0.0) ? r.bm_ms / r.gmp_ms : numeric_limits<double>::quiet_NaN();
    printf("%-10s %-22s %14.3f %14.3f %10.2fx %6d %8s\n",
           r.op.c_str(), r.size_label.c_str(), r.bm_ms, r.gmp_ms, ratio, r.reps, r.phase);
    fflush(stdout);
}

static void BenchMul(int digits_a, int digits_b, vector<Row> &rows)
{
    fprintf(stderr, "  [mul %dx%d] setup...\n", digits_a, digits_b); fflush(stderr);
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
    Row row{"mul", buf, bm, gp, iters};
    rows.push_back(row);
    PrintRow(row);

    mpz_clear(ga);
    mpz_clear(gb);
    mpz_clear(gc);
}

static void BenchDiv(int digits_a, int digits_b, vector<Row> &rows)
{
    fprintf(stderr, "  [div %dx%d] setup...\n", digits_a, digits_b); fflush(stderr);
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

    char buf[64];
    snprintf(buf, sizeof(buf), "%dx%d digits", digits_a, digits_b);
    auto check = [](const BigInteger &value, mpz_srcptr reference) {
        mpz_t z;
        mpz_init(z);
        mpz_import(z, value.size(), -1, sizeof(DataT), 0, 64 - LimbBits, value.Limbs().data());
        if (value.IsNegative()) mpz_neg(z, z);
        if (mpz_cmp(z, reference)) throw runtime_error("division result differs from GMP");
        mpz_clear(z);
    };
    for (bool withRemainder : {false, true}) {
        BigInteger q, r;
        for (bool warm : {false, true}) {
            int reps = warm ? iters : 1;
            double bm = TimeMs([&]() {
                if (withRemainder) tie(q, r) = DivideAndRemainder(a, b);
                else q = a / b;
            }, reps);
            double gp = TimeMs([&]() {
                if (withRemainder) mpz_tdiv_qr(gq, gr, ga, gb);
                else mpz_tdiv_q(gq, ga, gb);
            }, reps);
            check(q, gq);
            if (withRemainder) check(r, gr);
            rows.push_back({withRemainder ? "div_qr" : "div_q", buf, bm, gp, reps,
                            warm ? "warm" : "first"});
            PrintRow(rows.back());
        }
    }

    mpz_clear(ga);
    mpz_clear(gb);
    mpz_clear(gq);
    mpz_clear(gr);
}

static void BenchParse(int digits, vector<Row> &rows)
{
    fprintf(stderr, "  [parse %d] setup...\n", digits); fflush(stderr);
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
    Row row{"parse", buf, bm, gp, iters};
    rows.push_back(row);
    PrintRow(row);
}

static void BenchToString(int digits, vector<Row> &rows)
{
    fprintf(stderr, "  [tostr %d] setup...\n", digits); fflush(stderr);
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
    Row row{"tostr", buf, bm, gp, iters};
    rows.push_back(row);
    PrintRow(row);

    mpz_clear(g);
}

// Run one binary-origin fixture in a fresh process so the first call
// includes formatter setup and starts with empty decimal caches.
static void BenchBinaryFormat(int digits, uint64_t seed)
{
    mt19937_64 gen(seed ^ uint64_t(digits));
    vector<uint8_t> bytes((uint64_t(digits) * 3322 + 7999) / 8000);
    for (auto &byte : bytes)
        byte = uint8_t(gen());
    bytes.front() |= 128;
    const bool negative = (seed & 1) == 0;
    const auto value = BigInteger::FromByteArray(bytes, negative);
    mpz_t reference;
    mpz_init(reference);
    mpz_import(reference, bytes.size(), 1, 1, 0, 0, bytes.data());
    if (negative)
        mpz_neg(reference, reference);

    const auto once = [&](bool gmpFirst) {
        string actual;
        char *expected = nullptr;
        const auto big = [&] { actual = ToString(value); };
        const auto gmp = [&] { expected = mpz_get_str(nullptr, 10, reference); };
        double bigMs, gmpMs;
        if (gmpFirst) {
            gmpMs = TimeMs(gmp, 1);
            bigMs = TimeMs(big, 1);
        } else {
            bigMs = TimeMs(big, 1);
            gmpMs = TimeMs(gmp, 1);
        }
        if (expected == nullptr || actual != expected)
            throw runtime_error("ToString result differs from GMP");
        free(expected);
        return pair{bigMs, gmpMs};
    };

    const auto first = once(seed & 1);
    const int repetitions = digits <= 10000 ? 17 : 7;
    vector<double> big, gmp;
    for (int i = 0; i < repetitions; ++i) {
        const auto [bigMs, gmpMs] = once((i + seed) & 1);
        big.push_back(bigMs);
        gmp.push_back(gmpMs);
    }
    sort(big.begin(), big.end());
    sort(gmp.begin(), gmp.end());
    printf("approx_digits,seed,bigmath_first_ms,gmp_first_ms,bigmath_warm_ms,gmp_warm_ms\n");
    printf("%d,%llu,%.9f,%.9f,%.9f,%.9f\n", digits, static_cast<unsigned long long>(seed),
           first.first, first.second, big[repetitions / 2], gmp[repetitions / 2]);
    mpz_clear(reference);
}

static void BenchBigDecimalAdd(int int_digits, int frac_digits, vector<Row> &rows)
{
    string sa_int = GenerateRandomDigits(int_digits, 0x11ull ^ (uint64_t)int_digits);
    string sa_frac = GenerateRandomDigits(frac_digits, 0x22ull ^ (uint64_t)frac_digits);
    string sb_int = GenerateRandomDigits(int_digits, 0x33ull ^ (uint64_t)int_digits);
    string sb_frac = GenerateRandomDigits(frac_digits, 0x44ull ^ (uint64_t)frac_digits);

    string sa = sa_int + "." + sa_frac;
    string sb = sb_int + "." + sb_frac;

    BigDecimal a = BigDecimal::FromString(sa);
    BigDecimal b = BigDecimal::FromString(sb);

    unsigned long prec_bits = (unsigned long)((int_digits + frac_digits) * 3.322) + 64;
    mpf_t ga, gb, gc;
    mpf_init2(ga, prec_bits);
    mpf_init2(gb, prec_bits);
    mpf_init2(gc, prec_bits);
    mpf_set_str(ga, sa.c_str(), 10);
    mpf_set_str(gb, sb.c_str(), 10);

    int iters = IterationsForDigits(int_digits + frac_digits);

    double bm = TimeMs([&]() {
        BigDecimal c = a + b;
        if (c.Scale() < -1000000) abort();
    }, iters);

    double gp = TimeMs([&]() {
        mpf_add(gc, ga, gb);
    }, iters);

    char buf[64];
    snprintf(buf, sizeof(buf), "%d.%d digits", int_digits, frac_digits);
    Row row{"dec_add", buf, bm, gp, iters};
    rows.push_back(row);
    PrintRow(row);

    mpf_clear(ga);
    mpf_clear(gb);
    mpf_clear(gc);
}

static void BenchBigDecimalMul(int int_digits, int frac_digits, vector<Row> &rows)
{
    string sa_int = GenerateRandomDigits(int_digits, 0x55ull ^ (uint64_t)int_digits);
    string sa_frac = GenerateRandomDigits(frac_digits, 0x66ull ^ (uint64_t)frac_digits);
    string sb_int = GenerateRandomDigits(int_digits, 0x77ull ^ (uint64_t)int_digits);
    string sb_frac = GenerateRandomDigits(frac_digits, 0x88ull ^ (uint64_t)frac_digits);

    string sa = sa_int + "." + sa_frac;
    string sb = sb_int + "." + sb_frac;

    BigDecimal a = BigDecimal::FromString(sa);
    BigDecimal b = BigDecimal::FromString(sb);

    unsigned long prec_bits = (unsigned long)((int_digits + frac_digits) * 3.322) + 64;
    mpf_t ga, gb, gc;
    mpf_init2(ga, prec_bits);
    mpf_init2(gb, prec_bits);
    mpf_init2(gc, prec_bits);
    mpf_set_str(ga, sa.c_str(), 10);
    mpf_set_str(gb, sb.c_str(), 10);

    int iters = IterationsForDigits(int_digits + frac_digits);

    double bm = TimeMs([&]() {
        BigDecimal c = a * b;
        if (c.Scale() < -1000000) abort();
    }, iters);

    double gp = TimeMs([&]() {
        mpf_mul(gc, ga, gb);
    }, iters);

    char buf[64];
    snprintf(buf, sizeof(buf), "%d.%d digits", int_digits, frac_digits);
    Row row{"dec_mul", buf, bm, gp, iters};
    rows.push_back(row);
    PrintRow(row);

    mpf_clear(ga);
    mpf_clear(gb);
    mpf_clear(gc);
}

static void BenchBigDecimalDiv(int int_digits, int frac_digits, int scale, vector<Row> &rows)
{
    string sa_int = GenerateRandomDigits(int_digits, 0x99ull ^ (uint64_t)int_digits);
    string sa_frac = GenerateRandomDigits(frac_digits, 0xAAull ^ (uint64_t)frac_digits);
    string sb_int = GenerateRandomDigits(int_digits, 0xBBull ^ (uint64_t)int_digits);
    string sb_frac = GenerateRandomDigits(frac_digits, 0xCCull ^ (uint64_t)frac_digits);

    string sa = sa_int + "." + sa_frac;
    string sb = sb_int + "." + sb_frac;

    BigDecimal a = BigDecimal::FromString(sa);
    BigDecimal b = BigDecimal::FromString(sb);

    unsigned long prec_bits = (unsigned long)((int_digits + frac_digits + scale) * 3.322) + 64;
    mpf_t ga, gb, gc;
    mpf_init2(ga, prec_bits);
    mpf_init2(gb, prec_bits);
    mpf_init2(gc, prec_bits);
    mpf_set_str(ga, sa.c_str(), 10);
    mpf_set_str(gb, sb.c_str(), 10);

    int iters = IterationsForDigits(int_digits + frac_digits);

    double bm = TimeMs([&]() {
        BigDecimal c = a.Divide(b, scale, RoundingMode::HALF_EVEN);
        if (c.Scale() < -1000000) abort();
    }, iters);

    double gp = TimeMs([&]() {
        mpf_div(gc, ga, gb);
    }, iters);

    char buf[64];
    snprintf(buf, sizeof(buf), "%d.%d / %d dp", int_digits, frac_digits, scale);
    Row row{"dec_div", buf, bm, gp, iters};
    rows.push_back(row);
    PrintRow(row);

    mpf_clear(ga);
    mpf_clear(gb);
    mpf_clear(gc);
}

static void BenchBigDecimalParse(int digits, vector<Row> &rows)
{
    string s = GenerateRandomDigits(digits, 0xDDull ^ (uint64_t)digits);
    s.insert(s.size() / 2, ".");
    const char *cs = s.c_str();

    unsigned long prec_bits = (unsigned long)(digits * 3.322) + 64;
    int iters = IterationsForDigits(digits);

    double bm = TimeMs([&]() {
        BigDecimal b = BigDecimal::FromString(cs);
        if (b.Scale() < -1000000) abort();
    }, iters);

    double gp = TimeMs([&]() {
        mpf_t g;
        mpf_init2(g, prec_bits);
        mpf_set_str(g, cs, 10);
        mpf_clear(g);
    }, iters);

    char buf[64];
    snprintf(buf, sizeof(buf), "%d digits", digits);
    Row row{"dec_parse", buf, bm, gp, iters};
    rows.push_back(row);
    PrintRow(row);
}

static void BenchBigDecimalToString(int digits, vector<Row> &rows)
{
    string s = GenerateRandomDigits(digits, 0xEEull ^ (uint64_t)digits);
    s.insert(s.size() / 2, ".");
    BigDecimal b = BigDecimal::FromString(s);

    unsigned long prec_bits = (unsigned long)(digits * 3.322) + 64;
    mpf_t g;
    mpf_init2(g, prec_bits);
    mpf_set_str(g, s.c_str(), 10);

    int iters = IterationsForDigits(digits);
    if (digits >= 100000) iters = 1;

    double bm = TimeMs([&]() {
        string out = b.ToPlainString();
        if (out.empty()) abort();
    }, iters);

    double gp = TimeMs([&]() {
        mp_exp_t exp;
        char *out = mpf_get_str(nullptr, &exp, 10, 0, g);
        if (!out) abort();
        free(out);
    }, iters);

    char buf[64];
    snprintf(buf, sizeof(buf), "%d digits", digits);
    Row row{"dec_tostr", buf, bm, gp, iters};
    rows.push_back(row);
    PrintRow(row);

    mpf_clear(g);
}

int main(int argc, char **argv)
{
    if (argc == 4 && string_view(argv[1]) == "--format") {
        const auto read = [](string_view text, auto &value) {
            const auto end = text.data() + text.size();
            const auto result = from_chars(text.data(), end, value);
            return result.ec == errc{} && result.ptr == end;
        };
        int digits;
        uint64_t seed;
        if (read(argv[2], digits) && digits > 0 && read(argv[3], seed)) {
            BenchBinaryFormat(digits, seed);
            return 0;
        }
        cerr << "usage: bench_vs_gmp --format digits seed\n";
        return 2;
    }

    cout << "BigMath vs GMP " << gmp_version << "\n\n";
    PrintHeader();

    vector<Row> rows;

    if (argc != 1) {
        auto readDigits = [](string_view text, int &digits) {
            auto end = text.data() + text.size();
            auto result = from_chars(text.data(), end, digits);
            return result.ec == errc{} && result.ptr == end && digits > 0;
        };
        int a, b;
        if (argc == 4 && string_view(argv[1]) == "--div" &&
            readDigits(argv[2], a) && readDigits(argv[3], b)) {
            BenchDiv(a, b, rows);
            return 0;
        }
        cerr << "usage: bench_vs_gmp [--div dividend_digits divisor_digits | --format digits seed]\n";
        return 2;
    }

    // ===== Multiplication =====
    cout << "\n[multiplication, balanced]\n";
    PrintHeader();
    for (int d : {1000, 5000, 10000, 50000, 100000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000, 100000000, 200000000})
        BenchMul(d, d, rows);

    cout << "\n[multiplication, skewed]\n";
    PrintHeader();
    BenchMul(100000, 10000, rows);
    BenchMul(500000, 50000, rows);
    BenchMul(1000000, 100000, rows);
    BenchMul(2000000, 200000, rows);
    BenchMul(5000000, 500000, rows);
    BenchMul(10000000, 1000000, rows);
    BenchMul(20000000, 2000000, rows);
    BenchMul(50000000, 5000000, rows);
    BenchMul(100000000, 10000000, rows);
    BenchMul(200000000, 20000000, rows);

    // ===== Division =====
    cout << "\n[division, balanced]\n";
    PrintHeader();
    for (int d : {1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000})
        BenchDiv(d, d, rows);

    cout << "\n[division, skewed (a >> b) — Newton/BZ band]\n";
    PrintHeader();
    BenchDiv(40000, 10000, rows);
    BenchDiv(100000, 10000, rows);
    BenchDiv(200000, 50000, rows);
    BenchDiv(500000, 100000, rows);
    BenchDiv(1000000, 200000, rows);
    BenchDiv(2000000, 500000, rows);
    BenchDiv(5000000, 1000000, rows);
    BenchDiv(10000000, 2000000, rows);
    BenchDiv(20000000, 4000000, rows);
    BenchDiv(50000000, 10000000, rows);
    BenchDiv(100000000, 20000000, rows);
    BenchDiv(200000000, 40000000, rows);

    // ===== Decimal Parse =====
    cout << "\n[decimal parse]\n";
    PrintHeader();
    for (int d : {1000, 10000, 50000, 100000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000})
        BenchParse(d, rows);

    // ===== ToString =====
    cout << "\n[ToString (known BigMath bottleneck)]\n";
    PrintHeader();
    for (int d : {1000, 10000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000})
        BenchToString(d, rows);

    // ===== BigDecimal Arithmetic =====
    cout << "\n[BigDecimal arithmetic]\n";
    cout << "\n--- Addition ---\n";
    PrintHeader();
    for (int d : {100, 1000, 5000, 20000})
        BenchBigDecimalAdd(d, d/10, rows);

    cout << "\n--- Multiplication ---\n";
    PrintHeader();
    for (int d : {100, 1000, 5000, 20000})
        BenchBigDecimalMul(d, d/10, rows);

    cout << "\n--- Division ---\n";
    PrintHeader();
    for (int d : {100, 1000, 5000, 20000})
        BenchBigDecimalDiv(d, d/10, d/10, rows);

    cout << "\n--- Division at varying target scales (operand=2000 int + 200 frac) ---\n";
    PrintHeader();
    for (int scale : {0, 10, 100, 1000, 5000})
        BenchBigDecimalDiv(2000, 200, scale, rows);

    // ===== BigDecimal Parse =====
    cout << "\n[BigDecimal decimal parse]\n";
    PrintHeader();
    for (int d : {100, 1000, 10000, 50000})
        BenchBigDecimalParse(d, rows);

    // ===== BigDecimal ToString =====
    cout << "\n[BigDecimal ToString]\n";
    PrintHeader();
    for (int d : {100, 1000, 10000, 50000})
        BenchBigDecimalToString(d, rows);

    cout << "\nDone. " << rows.size() << " measurements.\n";
    return 0;
}
