// Stage 2: GMP reference run (cached per machine).
//
//   run_gmp <dataset_dir> <out_csv> ["machine label"]
//
// Times GMP on every manifest entry and records, for the exact-integer ops, a
// base-independent value hash so BigMath can be checked against it later. The
// result CSV is keyed to this machine by the driver, so GMP only runs once per
// machine per dataset.

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <gmp.h>

#include "bench_common.h"

using namespace bench;

// FNV-1a over the integer's LE 32-bit word stream, via mpz_export. Matches the
// BigMath limb hash bit-for-bit.
static std::string HashMpz(const mpz_t z) {
    Fnv64 h;
    size_t count = 0;
    // order=-1 (LS word first), size=4, endian=-1 (LE within word), nails=0.
    void *buf = mpz_export(nullptr, &count, -1, 4, -1, 0, z);
    if (buf) {
        h.Bytes((const uint8_t *)buf, count * 4);
        free(buf);
    }
    char out[20];
    snprintf(out, sizeof(out), "%016llx", (unsigned long long)h.h);
    return out;
}

static unsigned long PrecBits(int digits) {
    return (unsigned long)(digits * 3.322) + 64;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <dataset_dir> <out_csv> [machine]\n", argv[0]);
        return 2;
    }
    std::string dir = argv[1];
    std::string out_csv = argv[2];
    std::string machine = argc > 3 ? argv[3] : "unknown";

    std::vector<Entry> entries = ReadManifest(dir);
    std::vector<GmpResult> results(entries.size());

    fprintf(stderr, "[gmp] %s, %zu entries\n", gmp_version, entries.size());

    for (size_t i = 0; i < entries.size(); ++i) {
        const Entry &e = entries[i];
        GmpResult &r = results[i];
        r.hash = "-";
        fprintf(stderr, "  e%-4d %-9s %-22s ...\n", e.id, OpName(e.op), e.Label().c_str());

        switch (e.op) {
            case Op::Mul: {
                std::string sa = LoadOperand(dir, e.id, 0), sb = LoadOperand(dir, e.id, 1);
                mpz_t a, b, c; mpz_init_set_str(a, sa.c_str(), 10);
                mpz_init_set_str(b, sb.c_str(), 10); mpz_init(c);
                r.gmp_ms = BestMs([&] { mpz_mul(c, a, b); }, e.reps);
                r.hash = HashMpz(c);
                mpz_clear(a); mpz_clear(b); mpz_clear(c);
                break;
            }
            case Op::Div: {
                std::string sa = LoadOperand(dir, e.id, 0), sb = LoadOperand(dir, e.id, 1);
                mpz_t a, b, q, rem; mpz_init_set_str(a, sa.c_str(), 10);
                mpz_init_set_str(b, sb.c_str(), 10); mpz_init(q); mpz_init(rem);
                r.gmp_ms = BestMs([&] { mpz_tdiv_qr(q, rem, a, b); }, e.reps);
                r.hash = HashMpz(q);
                mpz_clear(a); mpz_clear(b); mpz_clear(q); mpz_clear(rem);
                break;
            }
            case Op::Parse: {
                std::string sa = LoadOperand(dir, e.id, 0);
                r.gmp_ms = BestMs([&] {
                    mpz_t g; mpz_init(g); mpz_set_str(g, sa.c_str(), 10); mpz_clear(g);
                }, e.reps);
                break;
            }
            case Op::ToStr: {
                std::string sa = LoadOperand(dir, e.id, 0);
                mpz_t g; mpz_init_set_str(g, sa.c_str(), 10);
                r.gmp_ms = BestMs([&] {
                    char *o = mpz_get_str(nullptr, 10, g); if (!o) abort(); free(o);
                }, e.reps);
                mpz_clear(g);
                break;
            }
            case Op::DecAdd: case Op::DecMul: case Op::DecDiv: {
                std::string sa = LoadOperand(dir, e.id, 0), sb = LoadOperand(dir, e.id, 1);
                unsigned long pb = PrecBits(e.a_digits + e.b_digits + e.scale);
                mpf_t a, b, c; mpf_init2(a, pb); mpf_init2(b, pb); mpf_init2(c, pb);
                mpf_set_str(a, sa.c_str(), 10); mpf_set_str(b, sb.c_str(), 10);
                if (e.op == Op::DecAdd)
                    r.gmp_ms = BestMs([&] { mpf_add(c, a, b); }, e.reps);
                else if (e.op == Op::DecMul)
                    r.gmp_ms = BestMs([&] { mpf_mul(c, a, b); }, e.reps);
                else
                    r.gmp_ms = BestMs([&] { mpf_div(c, a, b); }, e.reps);
                mpf_clear(a); mpf_clear(b); mpf_clear(c);
                break;
            }
            case Op::DecParse: {
                std::string sa = LoadOperand(dir, e.id, 0);
                unsigned long pb = PrecBits(e.a_digits);
                r.gmp_ms = BestMs([&] {
                    mpf_t g; mpf_init2(g, pb); mpf_set_str(g, sa.c_str(), 10); mpf_clear(g);
                }, e.reps);
                break;
            }
            case Op::DecToStr: {
                std::string sa = LoadOperand(dir, e.id, 0);
                unsigned long pb = PrecBits(e.a_digits);
                mpf_t g; mpf_init2(g, pb); mpf_set_str(g, sa.c_str(), 10);
                r.gmp_ms = BestMs([&] {
                    mp_exp_t exp; char *o = mpf_get_str(nullptr, &exp, 10, 0, g);
                    if (!o) abort(); free(o);
                }, e.reps);
                mpf_clear(g);
                break;
            }
        }
    }

    std::string comment = "GMP " + std::string(gmp_version) + " | machine=" + machine;
    WriteGmpResults(out_csv, entries, results, comment);
    fprintf(stderr, "[gmp] wrote %s\n", out_csv.c_str());
    return 0;
}
