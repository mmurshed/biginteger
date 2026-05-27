// MFA forward+inverse identity probe.
// Compares directly against non-MFA forward/inverse from the same file.

#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"

using namespace std;
using namespace BigMath;
using namespace BigMath::NttCrt;

static bool TestRoundTrip(Int n)
{
    mt19937_64 gen(0xDEADBEEFULL ^ (uint64_t)n);
    uniform_int_distribution<uint32_t> dis(0, P1 - 1);

    vector<UInt> orig(n), buf(n), scratch(n);
    for (Int i = 0; i < n; ++i) orig[i] = dis(gen);

    // MFA round trip
    buf = orig;
    ForwardMFA<F1, G1>(buf.data(), n, scratch.data());
    InverseMFA<F1, G1>(buf.data(), n, scratch.data());

    bool ok = true;
    Int firstFail = -1;
    for (Int i = 0; i < n; ++i)
    {
        if (buf[i] != orig[i])
        {
            ok = false;
            firstFail = i;
            break;
        }
    }
    printf("n=%-12d : %s", n, ok ? "OK" : "FAIL");
    if (!ok)
        printf("  firstFail=%d  got=%u expected=%u", firstFail, buf[firstFail], orig[firstFail]);
    printf("\n");
    return ok;
}

static bool TestNonMfaRoundTrip(Int n)
{
    mt19937_64 gen(0xCAFEBABEULL ^ (uint64_t)n);
    uniform_int_distribution<uint32_t> dis(0, P1 - 1);

    vector<UInt> orig(n), buf(n);
    for (Int i = 0; i < n; ++i) orig[i] = dis(gen);

    const auto &plan = GetPlan<F1, G1>(n);
    buf = orig;
    Forward<F1>(buf, plan);
    Inverse<F1>(buf, plan);

    bool ok = true;
    for (Int i = 0; i < n; ++i)
        if (buf[i] != orig[i]) { ok = false; break; }
    printf("nonMFA n=%-7d : %s\n", n, ok ? "OK" : "FAIL");
    return ok;
}

int main(int argc, char **argv)
{
    vector<Int> sizes;
    if (argc > 1)
        for (int k = 1; k < argc; ++k) sizes.push_back(atoi(argv[k]));
    else
        sizes = {1<<14, 1<<18, 1<<20, 1<<21, 1<<22, 1<<23};

    bool ok = true;
    for (Int n : sizes)
    {
        ok = TestNonMfaRoundTrip(n) && ok;
        ok = TestRoundTrip(n) && ok;
    }
    printf("\n%s\n", ok ? "ALL PASS" : "FAILURES");
    return ok ? 0 : 1;
}
