// Raw-limb GMP probe for the MFA multiplication band (mem_pass_fusion.md
// correctness protocol, step 3).
//
// Drives NttCrt::Multiply directly on random Base2_64 limb vectors and
// compares against GMP via mpz_import/mpz_export. Deliberately NEVER uses
// Parse/ToString: those route through the very multiplication band under
// test, so a corrupt MFA layer could cancel itself out and the comparison
// would prove nothing (the PR #103 lesson).
//
// Default tiers put the transform length at 2^24 and 2^25 — at and above the
// MFA gate — plus one 8:1 skewed shape for the lowered high-skew cutover.
// Runtime is minutes and memory is several GB; this is a pre-merge probe,
// not a ctest.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#include <gmp.h>

#include "biginteger/common/Constants.h"
#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"

using namespace BigMath;

static std::vector<DataT> RandLimbs(SizeT n, uint64_t seed)
{
    std::mt19937_64 g(seed);
    std::vector<DataT> v(n);
    for (SizeT i = 0; i < n; ++i) v[i] = g();
    if (v.back() == 0) v.back() = 1; // keep the top limb significant
    return v;
}

static bool Check(SizeT la, SizeT lb, uint64_t seed)
{
    std::vector<DataT> a = RandLimbs(la, seed);
    std::vector<DataT> b = RandLimbs(lb, seed ^ 0x9E3779B97F4A7C15ULL);

    std::vector<DataT> c = NttCrt::Multiply(a, b, Base2_64);

    mpz_t ga, gb, gc;
    mpz_init(ga); mpz_init(gb); mpz_init(gc);
    mpz_import(ga, la, -1, sizeof(DataT), 0, 0, a.data());
    mpz_import(gb, lb, -1, sizeof(DataT), 0, 0, b.data());
    mpz_mul(gc, ga, gb);

    size_t gn = 0;
    std::vector<DataT> g(la + lb, 0);
    mpz_export(g.data(), &gn, -1, sizeof(DataT), 0, 0, gc);
    mpz_clear(ga); mpz_clear(gb); mpz_clear(gc);

    bool ok = (c.size() == gn) && (memcmp(c.data(), g.data(), gn * sizeof(DataT)) == 0);
    printf("mul %zu x %zu limbs : %s  (bm_limbs=%zu, gmp_limbs=%zu)\n",
           (size_t)la, (size_t)lb, ok ? "OK" : "FAIL", c.size(), gn);
    if (!ok)
    {
        size_t lim = std::min(c.size(), gn);
        for (size_t i = 0; i < lim; ++i)
            if (c[i] != g[i])
            {
                printf("  first diff at limb %zu: bm=%016llx gmp=%016llx\n",
                       i, (unsigned long long)c[i], (unsigned long long)g[i]);
                break;
            }
    }
    return ok;
}

int main(int argc, char **argv)
{
    bool ok = true;
    if (argc > 1)
    {
        for (int k = 1; k + 1 < argc; k += 2)
            ok &= Check((SizeT)atoll(argv[k]), (SizeT)atoll(argv[k + 1]), 0xC0FFEE ^ k);
    }
    else
    {
        // 2 coeffs/limb: la+lb limbs → ~2(la+lb) coeffs → n = bit_ceil.
        ok &= Check(1u << 22, 1u << 22, 0xC0FFEE01); // n = 2^24 (MFA gate)
        ok &= Check(1u << 23, 1u << 23, 0xC0FFEE02); // n = 2^25
        ok &= Check(1u << 23, 1u << 20, 0xC0FFEE03); // 8:1 skew, n = 2^25
    }
    printf("%s\n", ok ? "ALL OK" : "FAILURES");
    return ok ? 0 : 1;
}
