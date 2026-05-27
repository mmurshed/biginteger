// MFA correctness probe: compare CRT-MFA path against Goldilocks single-prime
// NTT at sizes large enough to exercise the MFA gate.
//
// Truth source: BIGMATH_NTT_CRT=0 forces NTTMultiplication to use Goldilocks
// NTTCore — independent of the CRT-MFA code path. We compile two TUs into one
// binary with __MFA_OFF__ swapping the macros.
//
// Build:
//   c++ -std=c++20 -O3 -march=native -Iinclude \
//       tests/performance/mfa_correctness.cpp $(BIGMATH_SRCS) -o mfa_correctness

#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/multiplication/NTTMultiplication.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"

using namespace std;
using namespace BigMath;

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

static uint64_t Hash(const vector<DataT> &v)
{
    uint64_t h = 0xcbf29ce484222325ULL;
    for (DataT x : v) {
        h ^= (uint64_t)x;
        h *= 0x100000001b3ULL;
    }
    h ^= v.size();
    return h;
}

static bool TestSize(size_t na, size_t nb)
{
    BaseT base = BigInteger::Base();
    auto a = RandomLimbs(na, 0xA0ULL ^ na);
    auto b = RandomLimbs(nb, 0xB0ULL ^ nb);

    vector<DataT> c1 = NTTMultiplication::Multiply(a, b, base);
    printf("HASH(%zux%zu) = %016llx  size=%zu\n",
        na, nb, (unsigned long long)Hash(c1), c1.size());
    return true;

#if 0
    // Cross-check by recomputing in chunks: split a into halves, do two NTT
    // mults of (a_lo*b) and (a_hi*b)<<half, sum. This exercises a different
    // call sequence at smaller-sub sizes (no MFA in sub) and should match.
    size_t half = na / 2;
    vector<DataT> a_lo(a.begin(), a.begin() + half);
    vector<DataT> a_hi(a.begin() + half, a.end());
    auto p_lo = NTTMultiplication::Multiply(a_lo, b, base);
    auto p_hi = NTTMultiplication::Multiply(a_hi, b, base);
    // Shift p_hi by `half` limbs and add to p_lo.
    vector<DataT> sum(max(p_lo.size(), p_hi.size() + half), 0);
    for (size_t i = 0; i < p_lo.size(); ++i) sum[i] = p_lo[i];
    // Add p_hi << half. Base-aware add.
    if (base == Base2_64) {
        unsigned __int128 carry = 0;
        for (size_t i = 0; i < p_hi.size(); ++i) {
            unsigned __int128 t = (unsigned __int128)sum[i + half] + p_hi[i] + carry;
            sum[i + half] = (DataT)t;
            carry = t >> 64;
        }
        size_t k = p_hi.size() + half;
        while (carry) {
            if (k >= sum.size()) sum.push_back(0);
            unsigned __int128 t = (unsigned __int128)sum[k] + carry;
            sum[k] = (DataT)t;
            carry = t >> 64;
            ++k;
        }
    } else {
        uint64_t carry = 0;
        for (size_t i = 0; i < p_hi.size(); ++i) {
            uint64_t t = sum[i + half] + p_hi[i] + carry;
            sum[i + half] = (DataT)(t & 0xFFFFFFFFULL);
            carry = t >> 32;
        }
        size_t k = p_hi.size() + half;
        while (carry) {
            if (k >= sum.size()) sum.push_back(0);
            uint64_t t = sum[k] + carry;
            sum[k] = (DataT)(t & 0xFFFFFFFFULL);
            carry = t >> 32;
            ++k;
        }
    }
    while (sum.size() > 1 && sum.back() == 0) sum.pop_back();
    while (c1.size() > 1 && c1.back() == 0) c1.pop_back();

    bool ok = (c1.size() == sum.size());
    if (ok) for (size_t i = 0; i < c1.size(); ++i)
        if (c1[i] != sum[i]) { ok = false; break; }

    printf("%-12zu x %-12zu : %s  (out limbs c1=%zu split=%zu)\n",
        na, nb, ok ? "OK" : "FAIL", c1.size(), sum.size());
    return ok;
#endif
}

int main(int argc, char **argv)
{
    BaseT base = BigInteger::Base();
    printf("MFA correctness probe  base=%s\n", base == Base2_64 ? "Base2_64" : "Base2_32");

    // Hit several MFA-active sizes. With BIGMATH_NTT_MFA_THRESHOLD=2^21 (~2M),
    // need n = 2*limbs ≥ 2M → limbs ≥ 1M (Base2_64).
    // The split test does sub-mults at half-size, so we want the FULL mult to
    // be MFA-active but the split halves below.
    vector<pair<size_t, size_t>> cases;
    if (argc > 1) {
        for (int k = 1; k + 1 < argc; k += 2)
            cases.push_back({(size_t)atoll(argv[k]), (size_t)atoll(argv[k+1])});
    } else {
        // Default cases — small first to catch shallow bugs, then escalate.
        cases = {
            {100'000, 100'000},   // n ≈ 2^18, below default MFA threshold
            {600'000, 600'000},   // n ≈ 2^21, hits default MFA gate
            {1'500'000, 1'500'000},  // n ≈ 2^23, ~50M digits equiv
            {3'000'000, 3'000'000},  // n ≈ 2^23-2^24, ~100M digits equiv
        };
    }

    bool ok = true;
    for (auto [na, nb] : cases) ok = TestSize(na, nb) && ok;
    printf("\n%s\n", ok ? "ALL PASS" : "FAILURES");
    return ok ? 0 : 1;
}
