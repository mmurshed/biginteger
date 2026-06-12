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
// MFA gate — plus one 8:1 skewed shape for the lowered high-skew cutover,
// cyclic (MultiplyMod2km1) rows spanning the plain/MFA boundary and the old
// 2^22 cap, and two full-pipeline division tiers.
// Runtime is minutes and memory is several GB; this is a pre-merge probe,
// not a ctest.
//
// Transient-break runs: kill the probe after the cyc rows. The div rows do
// NOT fail fast on a corrupt multiply layer — Newton's fixup fallback spins
// for tens of minutes instead (the PR #103 signature), so a broken build
// hangs there rather than printing FAIL.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#include <gmp.h>

#include "biginteger/common/Constants.h"
#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"
#include "biginteger/algorithms/Division.h"

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

// Cyclic product check: NttCrt::MultiplyMod2km1 vs GMP's a*b mod (2^(64L)-1),
// folded via repeated shift-and-add. Exercises the cyclic MFA routing at
// n = 2*L (Base2_64).
static bool CheckCyclic(SizeT la, SizeT lb, SizeT L, uint64_t seed)
{
    std::vector<DataT> a = RandLimbs(la, seed);
    std::vector<DataT> b = RandLimbs(lb, seed ^ 0xDEADBEEFCAFEF00DULL);

    std::vector<DataT> c = NttCrt::MultiplyMod2km1(a, b, L, Base2_64);

    mpz_t ga, gb, gc, gm;
    mpz_init(ga); mpz_init(gb); mpz_init(gc); mpz_init(gm);
    mpz_import(ga, la, -1, sizeof(DataT), 0, 0, a.data());
    mpz_import(gb, lb, -1, sizeof(DataT), 0, 0, b.data());
    mpz_mul(gc, ga, gb);
    // gm = 2^(64L) - 1; fold gc into [0, gm) (gc < gm^2, two folds + final).
    mpz_set_ui(gm, 1); mpz_mul_2exp(gm, gm, 64 * L); mpz_sub_ui(gm, gm, 1);
    while (mpz_cmp(gc, gm) >= 0)
    {
        mpz_t hi, lo;
        mpz_init(hi); mpz_init(lo);
        mpz_fdiv_q_2exp(hi, gc, 64 * L);
        mpz_fdiv_r_2exp(lo, gc, 64 * L);
        mpz_add(gc, hi, lo);
        mpz_clear(hi); mpz_clear(lo);
    }

    size_t gn = 0;
    std::vector<DataT> g(L + 1, 0);
    mpz_export(g.data(), &gn, -1, sizeof(DataT), 0, 0, gc);
    if (gn == 0) { g.assign(1, 0); gn = 1; } // canonical zero
    mpz_clear(ga); mpz_clear(gb); mpz_clear(gc); mpz_clear(gm);

    bool ok = (c.size() == gn) && (memcmp(c.data(), g.data(), gn * sizeof(DataT)) == 0);
    printf("cyc %zu x %zu mod B^%zu-1 : %s  (bm_limbs=%zu, gmp_limbs=%zu)\n",
           (size_t)la, (size_t)lb, (size_t)L, ok ? "OK" : "FAIL", c.size(), gn);
    return ok;
}

// Division check: quotient limb-for-limb vs GMP. Exercises the full Newton
// pipeline including the (now MFA-routed) cyclic wrapped products.
static bool CheckDiv(SizeT la, SizeT lb, uint64_t seed)
{
    std::vector<DataT> a = RandLimbs(la, seed);
    std::vector<DataT> b = RandLimbs(lb, seed ^ 0x517CC1B727220A95ULL);

    std::vector<DataT> q = Divide(a, b, Base2_64);

    mpz_t ga, gb, gq;
    mpz_init(ga); mpz_init(gb); mpz_init(gq);
    mpz_import(ga, la, -1, sizeof(DataT), 0, 0, a.data());
    mpz_import(gb, lb, -1, sizeof(DataT), 0, 0, b.data());
    mpz_tdiv_q(gq, ga, gb);

    size_t gn = 0;
    std::vector<DataT> g(la - lb + 2, 0);
    mpz_export(g.data(), &gn, -1, sizeof(DataT), 0, 0, gq);
    mpz_clear(ga); mpz_clear(gb); mpz_clear(gq);

    SizeT qn = q.size();
    while (qn > 1 && q[qn - 1] == 0) --qn; // BigMath may keep a zero top limb
    bool ok = (qn == gn) && (memcmp(q.data(), g.data(), gn * sizeof(DataT)) == 0);
    printf("div %zu / %zu limbs : %s  (bm_limbs=%zu, gmp_limbs=%zu)\n",
           (size_t)la, (size_t)lb, ok ? "OK" : "FAIL", (size_t)qn, gn);
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
        ok &= Check(1u << 22, 1u << 22, 0xC0FFEE01); // n = 2^24
        ok &= Check(1u << 23, 1u << 23, 0xC0FFEE02); // n = 2^25
        ok &= Check(1u << 23, 1u << 20, 0xC0FFEE03); // 8:1 skew, n = 2^25
        // Cyclic: n = 2L. Below the MFA gate (plain), at it, above it, and
        // above the old 2^22 cap (newly admissible after the 2^26 raise).
        ok &= CheckCyclic(1u << 18, 1u << 18, 1u << 18, 0xC0FFEE04); // n = 2^19, plain
        ok &= CheckCyclic(1u << 19, 1u << 19, 1u << 19, 0xC0FFEE05); // n = 2^20, MFA gate
        ok &= CheckCyclic(1u << 21, 1u << 21, 1u << 21, 0xC0FFEE06); // n = 2^22, MFA
        ok &= CheckCyclic(1u << 22, 1u << 22, 1u << 22, 0xC0FFEE07); // n = 2^23, beyond old cap
        ok &= CheckCyclic((1u << 22) + 12345, (1u << 22) - 777, 1u << 23, 0xC0FFEE08); // n = 2^24, uneven
        // Division through the full Newton pipeline at MFA-band sizes.
        ok &= CheckDiv(5200000, 1040000, 0xC0FFEE09);   // ≈100M ÷ 20M digits
        ok &= CheckDiv(10400000, 2080000, 0xC0FFEE0A);  // ≈200M ÷ 40M digits
    }
    printf("%s\n", ok ? "ALL OK" : "FAILURES");
    return ok ? 0 : 1;
}
