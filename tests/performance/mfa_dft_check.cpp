// Compare MFA forward against brute-force DFT (up to bit-reversal permutation).

#include <cstdio>
#include <random>
#include <vector>

#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"

using namespace std;
using namespace BigMath;
using namespace BigMath::NttCrt;

// Brute-force DFT mod P1 using planN1's root: roots[1] = w_n^1.
static vector<UInt> BruteDft(const vector<UInt> &x, const Plan<F1> &plan)
{
    Int n = (Int)x.size();
    UInt w = plan.forwardRoots[1];  // w_n^1
    vector<UInt> X((SizeT)n);
    for (Int k = 0; k < n; ++k)
    {
        UInt acc = 0;
        UInt cur = 1;  // w^{n*k} for varying n, starting at n=0
        UInt step = F1::Power(w, (UInt)k);  // w^k
        for (Int j = 0; j < n; ++j)
        {
            acc = F1::Add(acc, F1::Mul(x[(SizeT)j], cur));
            cur = F1::Mul(cur, step);
        }
        X[(SizeT)k] = acc;
    }
    return X;
}

static bool TestMul(Int n)
{
    mt19937_64 gen(0xABBA5EEDULL ^ (uint64_t)n);
    uniform_int_distribution<uint32_t> dis(0, P1 - 1);
    vector<UInt> a((SizeT)n / 2), b((SizeT)n / 2);
    for (auto &v : a) v = dis(gen);
    for (auto &v : b) v = dis(gen);

    // Brute cyclic convolution.
    vector<UInt> brute((SizeT)n, 0);
    for (Int i = 0; i < (Int)a.size(); ++i)
        for (Int j = 0; j < (Int)b.size(); ++j) {
            Int k = (i + j) % n;
            brute[(SizeT)k] = F1::Add(brute[(SizeT)k], F1::Mul(a[(SizeT)i], b[(SizeT)j]));
        }

    // MFA mul.
    vector<UInt> aBuf((SizeT)n, 0), bBuf((SizeT)n, 0), aSc((SizeT)n), bSc((SizeT)n);
    for (Int i = 0; i < (Int)a.size(); ++i) aBuf[(SizeT)i] = a[(SizeT)i];
    for (Int i = 0; i < (Int)b.size(); ++i) bBuf[(SizeT)i] = b[(SizeT)i];
    ForwardMFA<F1, G1>(aBuf.data(), n, aSc.data());
    ForwardMFA<F1, G1>(bBuf.data(), n, bSc.data());
    for (Int i = 0; i < n; ++i) aBuf[(SizeT)i] = F1::Mul(aBuf[(SizeT)i], bBuf[(SizeT)i]);
    InverseMFA<F1, G1>(aBuf.data(), n, aSc.data());

    Int wrong = 0;
    for (Int i = 0; i < n; ++i)
        if (aBuf[(SizeT)i] != brute[(SizeT)i]) ++wrong;
    printf("MUL n=%-7d : wrong=%-5d / %d\n", n, wrong, n);
    return wrong == 0;
}

int main(int argc, char **argv)
{
    if (argc > 1 && string(argv[1]) == "mul") {
        for (int k = 2; k < argc; ++k) TestMul(atoi(argv[k]));
        return 0;
    }
    Int n = argc > 1 ? atoi(argv[1]) : 128;
    mt19937_64 gen(0xCAFEDEEDULL);
    uniform_int_distribution<uint32_t> dis(0, P1 - 1);
    vector<UInt> orig((SizeT)n), buf((SizeT)n), scratch((SizeT)n);
    for (Int i = 0; i < n; ++i) orig[(SizeT)i] = dis(gen);

    const auto &plan = GetPlan<F1, G1>(n);
    auto dft = BruteDft(orig, plan);

    // MFA forward
    buf = orig;
    ForwardMFA<F1, G1>(buf.data(), n, scratch.data());

    // Verify: MFA output should equal DFT at some permuted index.
    // For each storage position p, find k such that dft[k] == buf[p].
    // Build the permutation.
    printf("n=%d\n", n);
    Int mismatches = 0;
    Int show = 16;
    for (Int p = 0; p < n; ++p)
    {
        // Find a matching dft index.
        Int found = -1;
        for (Int k = 0; k < n; ++k)
            if (dft[(SizeT)k] == buf[(SizeT)p]) { found = k; break; }
        if (found < 0) {
            if (mismatches < show)
                printf("  p=%-4d storage=%u  no DFT index matches\n", p, buf[(SizeT)p]);
            ++mismatches;
        } else if (p < show) {
            printf("  p=%-4d → k=%-4d  (storage=%u)\n", p, found, buf[(SizeT)p]);
        }
    }
    printf("\nMismatches: %d / %d\n", mismatches, n);

    // Compare against expected bit-reversed mapping.
    Int logn = __builtin_ctz((unsigned)n);
    Int wrongMap = 0;
    for (Int p = 0; p < n; ++p)
    {
        Int br = 0;
        for (Int b = 0; b < logn; ++b)
            br |= ((p >> b) & 1) << (logn - 1 - b);
        if (buf[(SizeT)p] != dft[(SizeT)br]) {
            if (wrongMap < show)
                printf("  WRONG: p=%-4d  storage=%u  expected dft[br(p)=%d]=%u\n",
                       p, buf[(SizeT)p], br, dft[(SizeT)br]);
            ++wrongMap;
        }
    }
    printf("Bit-rev mapping mismatches: %d / %d\n", wrongMap, n);
    return mismatches == 0 ? 0 : 1;
}
