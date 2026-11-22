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

#include "case_measurement.h"

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

static unsigned long PrecBits(long long digits) {
    return (unsigned long)(digits * 3.322) + 64;
}

// The isolated runs load the exact same integer values from GMP raw files.
static void LoadInteger(mpz_t value, const std::string &dir, int id, int slot, bool raw) {
    mpz_init(value);
    if (raw) {
        const std::string path = RawOperandPath(dir, id, slot);
        FILE *input = std::fopen(path.c_str(), "rb");
        if (!input) throw std::runtime_error("cannot open " + path);
        const size_t bytes = mpz_inp_raw(value, input);
        std::fclose(input);
        if (!bytes) throw std::runtime_error("invalid GMP raw operand " + path);
    } else if (mpz_set_str(value, LoadOperand(dir, id, slot).c_str(), 10) != 0) {
        throw std::runtime_error("invalid decimal integer");
    }
}

static int ExportInteger(const std::string &dir, int id, int slot, const char *path) {
    mpz_t n;
    LoadInteger(n, dir, id, slot, false);
    FILE *output = std::fopen(path, "wb");
    if (!output) throw std::runtime_error("cannot open raw operand for writing");
    const size_t bytes = mpz_out_raw(output, n);
    const int close_status = std::fclose(output);
    mpz_clear(n);
    return bytes && !close_status ? 0 : 1;
}

static void MeasureCase(const std::string &dir, const Entry &e, CaseMeasurement &result) {
    switch (e.op) {
        case Op::Mul: {
            mpz_t a, b, c;
            LoadInteger(a, dir, e.id, 0, result.memory);
            LoadInteger(b, dir, e.id, 1, result.memory);
            mpz_init(c);
            result.Run(e.reps, [&] { mpz_mul(c, a, b); });
            result.hash = HashMpz(c);
            mpz_clear(a); mpz_clear(b); mpz_clear(c);
            break;
        }
        case Op::Div:
        case Op::DivQR: {
            mpz_t a, b, q, rem;
            LoadInteger(a, dir, e.id, 0, result.memory);
            LoadInteger(b, dir, e.id, 1, result.memory);
            mpz_init(q); mpz_init(rem);
            if (e.op == Op::Div)
                result.Run(e.reps, [&] { mpz_tdiv_q(q, a, b); });
            else
                result.Run(e.reps, [&] { mpz_tdiv_qr(q, rem, a, b); });
            result.hash = HashMpz(q);
            if (e.op == Op::DivQR) result.hash += "/" + HashMpz(rem);
            mpz_clear(a); mpz_clear(b); mpz_clear(q); mpz_clear(rem);
            break;
        }
        case Op::Parse: {
            std::string sa = LoadOperand(dir, e.id, 0);
            result.Run(e.reps, [&] {
                mpz_t g; mpz_init(g); mpz_set_str(g, sa.c_str(), 10); mpz_clear(g);
            });
            break;
        }
        case Op::ToStr: {
            mpz_t g;
            LoadInteger(g, dir, e.id, 0, result.memory);
            result.Run(e.reps, [&] {
                char *o = mpz_get_str(nullptr, 10, g); if (!o) abort(); free(o);
            });
            mpz_clear(g);
            break;
        }
        case Op::DecAdd: case Op::DecMul: case Op::DecDiv: {
            std::string sa = LoadOperand(dir, e.id, 0), sb = LoadOperand(dir, e.id, 1);
            unsigned long pb = PrecBits(e.a_digits + e.b_digits + e.scale);
            mpf_t a, b, c; mpf_init2(a, pb); mpf_init2(b, pb); mpf_init2(c, pb);
            mpf_set_str(a, sa.c_str(), 10); mpf_set_str(b, sb.c_str(), 10);
            if (e.op == Op::DecAdd)
                result.Run(e.reps, [&] { mpf_add(c, a, b); });
            else if (e.op == Op::DecMul)
                result.Run(e.reps, [&] { mpf_mul(c, a, b); });
            else
                result.Run(e.reps, [&] { mpf_div(c, a, b); });
            mpf_clear(a); mpf_clear(b); mpf_clear(c);
            break;
        }
        case Op::DecParse: {
            std::string sa = LoadOperand(dir, e.id, 0);
            unsigned long pb = PrecBits(e.a_digits);
            result.Run(e.reps, [&] {
                mpf_t g; mpf_init2(g, pb); mpf_set_str(g, sa.c_str(), 10); mpf_clear(g);
            });
            break;
        }
        case Op::DecToStr: {
            std::string sa = LoadOperand(dir, e.id, 0);
            unsigned long pb = PrecBits(e.a_digits);
            mpf_t g; mpf_init2(g, pb); mpf_set_str(g, sa.c_str(), 10);
            result.Run(e.reps, [&] {
                mp_exp_t exp; char *o = mpf_get_str(nullptr, &exp, 10, 0, g);
                if (!o) abort(); free(o);
            });
            mpf_clear(g);
            break;
        }
    }
}

int main(int argc, char **argv) {
    if (argc == 4 && std::string(argv[1]) == "--memory-case")
        return RunMemoryCase(argv[2], argv[3], MeasureCase);
    if (argc == 6 && std::string(argv[1]) == "--export-integer")
        return ExportInteger(argv[2], std::stoi(argv[3]), std::stoi(argv[4]), argv[5]);
    if (argc < 3) {
        fprintf(stderr, "usage: %s <dataset_dir> <out_csv> [machine]\n", argv[0]);
        return 2;
    }
    std::string dir = argv[1];
    std::string out_csv = argv[2];
    std::string machine = argc > 3 ? argv[3] : "unknown";

    std::vector<Entry> entries = TimingCases(dir);
    std::vector<GmpResult> results(entries.size());

    fprintf(stderr, "[gmp] %s, %zu entries\n", gmp_version, entries.size());

    for (size_t i = 0; i < entries.size(); ++i) {
        const Entry &e = entries[i];
        fprintf(stderr, "  e%-4d %-9s %-22s ...\n", e.id, OpName(e.op), e.Label().c_str());

        CaseMeasurement measurement;
        MeasureCase(dir, e, measurement);
        results[i] = {measurement.first_ms, measurement.warm_ms, measurement.hash};
    }

    std::string comment = "GMP " + std::string(gmp_version) + " | machine=" + machine;
    WriteGmpResults(out_csv, entries, results, comment);
    fprintf(stderr, "[gmp] wrote %s\n", out_csv.c_str());
    return 0;
}
