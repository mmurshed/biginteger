// Stage 3: BigMath benchmark run.
//
//   run_bigmath <dataset_dir> <gmp_csv> <out_csv> <run_id> <machine> <profile> <history_md>
//
// Times BigMath on the dataset, verifies the exact-integer ops against the
// cached GMP hashes, and reports BigMath/GMP ratios. Writes a per-run CSV
// (incrementing run id) and appends the entries flagged `historical` in the
// manifest to a long-lived comparison table.

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "bench_common.h"

#include "biginteger/BigInteger.h"
#include "biginteger/ops/Operations.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "bigdecimal/BigDecimal.h"

using namespace bench;
using namespace BigMath;

// FNV-1a over the BigInteger value's little-endian 32-bit word stream, matching
// HashMpz (mpz_export, size=4, order=-1, endian=-1). BigMath limbs are base-2^32
// or base-2^64 depending on build (Base2_64 == 0 sentinel = full 64-bit limbs):
// emit one 32-bit word per limb in base-2^32, two (low then high) in base-2^64,
// then drop trailing zero words so the stream matches GMP's trimmed export.
static std::string HashBigInteger(const BigInteger &v) {
    const std::vector<DataT> &limbs = v.GetInteger();
    const bool limb64 = (BigInteger::Base() == Base2_64);
    std::vector<uint32_t> words;
    words.reserve(limbs.size() * (limb64 ? 2 : 1));
    for (DataT limb : limbs) {
        words.push_back((uint32_t)limb);
        if (limb64) words.push_back((uint32_t)(limb >> 32));
    }
    while (!words.empty() && words.back() == 0) words.pop_back();
    Fnv64 h;
    for (uint32_t w : words) h.Word32(w);
    char out[20];
    snprintf(out, sizeof(out), "%016llx", (unsigned long long)h.h);
    return out;
}

struct Row {
    Entry e;
    double bm_ms = 0;
    double gmp_ms = 0;
    std::string correct = "n/a"; // ok / FAIL / n/a
};

int main(int argc, char **argv) {
    if (argc < 8) {
        fprintf(stderr, "usage: %s <dataset_dir> <gmp_csv> <out_csv> <run_id> "
                        "<machine> <profile> <history_md>\n", argv[0]);
        return 2;
    }
    std::string dir = argv[1], gmp_csv = argv[2], out_csv = argv[3];
    int run_id = std::stoi(argv[4]);
    std::string machine = argv[5], profile = argv[6], history_md = argv[7];

    std::vector<Entry> entries = ReadManifest(dir);
    std::vector<GmpResult> gmp = ReadGmpResults(gmp_csv);
    if (gmp.size() != entries.size()) {
        fprintf(stderr, "[bigmath] gmp results (%zu) != manifest entries (%zu); "
                        "regenerate the GMP run for this dataset\n",
                gmp.size(), entries.size());
        return 1;
    }

    printf("%-9s %-22s %14s %14s %10s %8s\n",
           "op", "size", "BigMath ms", "GMP ms", "BM/GMP", "check");
    printf("%-9s %-22s %14s %14s %10s %8s\n",
           "----", "----", "----------", "------", "------", "-----");
    fflush(stdout);

    std::vector<Row> rows(entries.size());
    int failures = 0;

    for (size_t i = 0; i < entries.size(); ++i) {
        const Entry &e = entries[i];
        Row &row = rows[i];
        row.e = e;
        row.gmp_ms = gmp[i].gmp_ms;
        fprintf(stderr, "  e%-4d %-9s %-22s ...\n", e.id, OpName(e.op), e.Label().c_str());

        switch (e.op) {
            case Op::Mul: {
                BigInteger a = Parse(LoadOperand(dir, e.id, 0).c_str());
                BigInteger b = Parse(LoadOperand(dir, e.id, 1).c_str());
                BigInteger c;
                row.bm_ms = BestMs([&] { c = a * b; if (c.GetInteger().empty()) abort(); }, e.reps);
                row.correct = (HashBigInteger(c) == gmp[i].hash) ? "ok" : "FAIL";
                break;
            }
            case Op::Div: {
                BigInteger a = Parse(LoadOperand(dir, e.id, 0).c_str());
                BigInteger b = Parse(LoadOperand(dir, e.id, 1).c_str());
                BigInteger q;
                row.bm_ms = BestMs([&] { q = a / b; if (q.GetInteger().empty()) abort(); }, e.reps);
                row.correct = (HashBigInteger(q) == gmp[i].hash) ? "ok" : "FAIL";
                break;
            }
            case Op::Parse: {
                std::string sa = LoadOperand(dir, e.id, 0);
                row.bm_ms = BestMs([&] {
                    BigInteger b = Parse(sa.c_str()); if (b.GetInteger().empty()) abort();
                }, e.reps);
                break;
            }
            case Op::ToStr: {
                BigInteger b = Parse(LoadOperand(dir, e.id, 0).c_str());
                row.bm_ms = BestMs([&] {
                    std::string o = ToString(b); if (o.empty()) abort();
                }, e.reps);
                break;
            }
            case Op::DecAdd: case Op::DecMul: case Op::DecDiv: {
                BigDecimal a = BigDecimal::FromString(LoadOperand(dir, e.id, 0));
                BigDecimal b = BigDecimal::FromString(LoadOperand(dir, e.id, 1));
                BigDecimal c;
                if (e.op == Op::DecAdd)
                    row.bm_ms = BestMs([&] { c = a + b; }, e.reps);
                else if (e.op == Op::DecMul)
                    row.bm_ms = BestMs([&] { c = a * b; }, e.reps);
                else
                    row.bm_ms = BestMs([&] { c = a.Divide(b, e.scale, RoundingMode::HALF_EVEN); }, e.reps);
                if (c.Scale() < -1000000) abort();
                break;
            }
            case Op::DecParse: {
                std::string sa = LoadOperand(dir, e.id, 0);
                row.bm_ms = BestMs([&] {
                    BigDecimal b = BigDecimal::FromString(sa); if (b.Scale() < -1000000) abort();
                }, e.reps);
                break;
            }
            case Op::DecToStr: {
                BigDecimal b = BigDecimal::FromString(LoadOperand(dir, e.id, 0));
                row.bm_ms = BestMs([&] {
                    std::string o = b.ToPlainString(); if (o.empty()) abort();
                }, e.reps);
                break;
            }
        }
        if (row.correct == "FAIL") ++failures;

        double ratio = row.gmp_ms > 0 ? row.bm_ms / row.gmp_ms : 0.0;
        printf("%-9s %-22s %14.3f %14.3f %9.2fx %8s\n",
               OpName(e.op), e.Label().c_str(), row.bm_ms, row.gmp_ms, ratio,
               row.correct.c_str());
        fflush(stdout);
    }

    // ── per-run CSV ──────────────────────────────────────────────────────────
    {
        std::ofstream out(out_csv);
        out << "# run_id=" << run_id << " machine=" << machine
            << " profile=" << profile << "\n";
        out << "id,op,label,reps,bigmath_ms,gmp_ms,ratio,correct\n";
        for (const Row &r : rows) {
            double ratio = r.gmp_ms > 0 ? r.bm_ms / r.gmp_ms : 0.0;
            out << r.e.id << ',' << OpName(r.e.op) << ',' << r.e.Label() << ','
                << r.e.reps << ',' << r.bm_ms << ',' << r.gmp_ms << ','
                << ratio << ',' << r.correct << '\n';
        }
        fprintf(stderr, "[bigmath] wrote %s\n", out_csv.c_str());
    }

    // ── append flagged rows to the long-lived history table ──────────────────
    {
        bool fresh = !std::ifstream(history_md).good();
        std::ofstream hist(history_md, std::ios::app);
        if (fresh) {
            hist << "# BigMath vs GMP — historical comparison\n\n"
                 << "Selected large-size results, one block per run. Lower ratio is better "
                 << "(BigMath ms / GMP ms).\n";
        }
        hist << "\n## run " << run_id << " — machine=" << machine
             << " profile=" << profile << "\n\n"
             << "| op | size | BigMath ms | GMP ms | BM/GMP | check |\n"
             << "|----|------|-----------:|-------:|-------:|:-----:|\n";
        for (const Row &r : rows) {
            if (!r.e.historical) continue;
            double ratio = r.gmp_ms > 0 ? r.bm_ms / r.gmp_ms : 0.0;
            char line[256];
            snprintf(line, sizeof(line), "| %s | %s | %.3f | %.3f | %.2fx | %s |\n",
                     OpName(r.e.op), r.e.Label().c_str(), r.bm_ms, r.gmp_ms, ratio,
                     r.correct.c_str());
            hist << line;
        }
        fprintf(stderr, "[bigmath] appended history -> %s\n", history_md.c_str());
    }

    if (failures) {
        fprintf(stderr, "[bigmath] %d correctness FAILURE(s) vs GMP\n", failures);
        return 1;
    }
    return 0;
}
