// Verify BigMath vs GMP for correctness at sizes that exercise CRT NTT.

#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>

#include <gmp.h>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/Operations.h"

using namespace std;
using namespace BigMath;

static string RandDigits(int n, uint64_t seed) {
    mt19937_64 g(seed);
    uniform_int_distribution<int> dis(0, 9), first(1, 9);
    string s;
    s.reserve(n);
    s += char('0' + first(g));
    for (int i = 1; i < n; ++i) s += char('0' + dis(g));
    return s;
}

static bool Check(int da, int db) {
    string sa = RandDigits(da, 0xA1ULL ^ da);
    string sb = RandDigits(db, 0xB2ULL ^ db ^ ((uint64_t)da << 24));
    BigInteger a = Parse(sa.c_str());
    BigInteger b = Parse(sb.c_str());
    BigInteger c = a * b;

    mpz_t ga, gb, gc;
    mpz_init_set_str(ga, sa.c_str(), 10);
    mpz_init_set_str(gb, sb.c_str(), 10);
    mpz_init(gc);
    mpz_mul(gc, ga, gb);

    char *gmpStr = mpz_get_str(nullptr, 10, gc);
    string bigStr = ToString(c);
    bool ok = (bigStr == gmpStr);
    printf("mul %d x %d : %s  (bm_digits=%zu, gmp_digits=%zu)\n",
        da, db, ok ? "OK" : "FAIL", bigStr.size(), strlen(gmpStr));
    if (!ok) {
        printf("  bm_first50: %.50s\n", bigStr.c_str());
        printf("  gmp_first50: %.50s\n", gmpStr);
        printf("  bm_last50: %s\n", bigStr.c_str() + (bigStr.size() > 50 ? bigStr.size() - 50 : 0));
        printf("  gmp_last50: %s\n", gmpStr + (strlen(gmpStr) > 50 ? strlen(gmpStr) - 50 : 0));
    }
    free(gmpStr);
    mpz_clear(ga); mpz_clear(gb); mpz_clear(gc);
    return ok;
}

int main(int argc, char **argv) {
    bool ok = true;
    if (argc > 1) {
        for (int k = 1; k + 1 < argc; k += 2)
            ok = Check(atoi(argv[k]), atoi(argv[k+1])) && ok;
    } else {
        // Sizes spanning the CRT-prime length threshold:
        // p2 supports NTT up to 2^22 ≈ 4M coeffs → up to ~20M digits operands.
        // Beyond that we expect the existing CRT path to break.
        for (int d : {1'000'000, 5'000'000, 10'000'000, 20'000'000, 30'000'000})
            ok = Check(d, d) && ok;
    }
    printf("\n%s\n", ok ? "ALL OK" : "FAILURES");
    return ok ? 0 : 1;
}
