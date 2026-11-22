// Stage 3: BigMath benchmark run.
//
//   run_bigmath <dataset_dir> <gmp_csv> <out_csv> <run_id> <machine> <profile> <history_md>
//
// Times BigMath on the dataset, verifies the exact-integer ops against the
// cached GMP hashes, and reports BigMath/GMP ratios. Writes a per-run CSV
// (incrementing run id) and appends the entries flagged `historical` in the
// manifest to a long-lived comparison table.

#include <charconv>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include "case_measurement.h"

#include "biginteger/BigInteger.h"
#include "biginteger/ops/Operations.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/algorithms/multiplication/FermatOuterMultiplication.h"
#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"
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
    double first_ms = 0, gmp_first_ms = 0;
    double bm_ms = 0;
    double gmp_ms = 0;
    std::string correct = "n/a"; // ok / FAIL / n/a
};

// GMP raw magnitudes have a four-byte big-endian length followed by that
// many big-endian bytes. Decode into the final limb vector in bounded chunks.
static BigInteger LoadInteger(const std::string &dir, int id, int slot, bool raw) {
    if (!raw) return Parse(LoadOperand(dir, id, slot).c_str());
    const std::string path = RawOperandPath(dir, id, slot);
    std::ifstream input(path, std::ios::binary);
    unsigned char header[4];
    if (!input.read(reinterpret_cast<char *>(header), 4) || (header[0] & 0x80))
        throw std::runtime_error("expected a nonnegative GMP raw operand: " + path);
    size_t bytes = 0;
    for (unsigned char b : header) bytes = (bytes << 8) | b;
    const size_t limb_bytes = BigInteger::Base() == Base2_64 ? 8 : 4;
    std::vector<DataT> limbs((bytes + limb_bytes - 1) / limb_bytes);
    unsigned char chunk[32768];
    size_t left = bytes;
    while (left) {
        const size_t count = std::min(left, sizeof(chunk));
        if (!input.read(reinterpret_cast<char *>(chunk), count))
            throw std::runtime_error("truncated raw operand: " + path);
        for (size_t i = 0; i < count; ++i) {
            const size_t pos = --left;
            limbs[pos / limb_bytes] |= DataT(chunk[i]) << (8 * (pos % limb_bytes));
        }
    }
    return BigInteger(std::move(limbs), false);
}

template <class Product>
static void MeasureProduct(const std::string &dir, const Entry &e, CaseMeasurement &result,
                           const Product &product) {
    BigInteger a = LoadInteger(dir, e.id, 0, result.memory);
    BigInteger b = LoadInteger(dir, e.id, 1, result.memory);
    BigInteger c;
    result.Run(e.reps, [&] { c = product(a, b); if (c.Limbs().empty()) abort(); });
    result.hash = HashBigInteger(c);
}

static BigInteger MultiplyForMemoryTest(const BigInteger &a, const BigInteger &b,
                                        std::string_view mode, unsigned setting) {
    if (BigInteger::Base() != Base2_64)
        throw std::invalid_argument("Fermat comparisons require 64-bit limbs");
    const auto left = a.AsVector(), right = b.AsVector();
    const auto &x = left.get(), &y = right.get();
    if (mode == "--fermat") {
        const auto p = setting ? FermatOuter::MakeParameters(x.size(), y.size(), setting)
                               : FermatOuter::SelectParameters(x.size(), y.size());
        if (!p) throw std::invalid_argument("Fermat ring cannot hold the requested product");
        return BigInteger(FermatOuter::MultiplyFermat(x, y, *p), false);
    }

#if BIGMATH_NTT_WIDE && BIGMATH_NTT_TRUNC && BIGMATH_NTT_MFA && BIGMATH_NTT_MFA_FUSE
    if (setting < 2 || x.size() != y.size() || setting > x.size())
        throw std::invalid_argument("wide-block comparisons require equal lengths and at least two blocks");
    const size_t block = FermatOuter::DivideUp(x.size(), setting);
    if (FermatOuter::DivideUp(x.size(), block) != setting)
        throw std::invalid_argument("operand cannot be spread across that many equal-width blocks");
    // Choose the block geometry without changing the library's global
    // thresholds. Each eligible pair uses the ordinary wide-prime path.
    auto partial = [](const auto &lhs, const auto &rhs) {
        if (2 * (lhs.size() + rhs.size()) - 1 > (size_t)NttCrt::Wide::WideResidueLength)
            return NttCrt::Wide::MultiplyWide(lhs, rhs, Base2_64);
        return NttCrt::Multiply(lhs, rhs, Base2_64);
    };
    return BigInteger(NttCrt::MultiplyBlocked(x, y, Base2_64, 2 * block, block, partial), false);
#else
    throw std::invalid_argument("wide-block comparisons require the fused wide CRT build");
#endif
}

static void MeasureCase(const std::string &dir, const Entry &e, CaseMeasurement &result) {
    switch (e.op) {
        case Op::Mul: {
            MeasureProduct(dir, e, result, [](const auto &a, const auto &b) { return a * b; });
            break;
        }
        case Op::Div:
        case Op::DivQR: {
            BigInteger a = LoadInteger(dir, e.id, 0, result.memory);
            BigInteger b = LoadInteger(dir, e.id, 1, result.memory);
            std::string answer;
            if (e.op == Op::Div) {
                BigInteger q;
                result.Run(e.reps, [&] { q = a / b; });
                answer = HashBigInteger(q);
            } else {
                std::pair<BigInteger, BigInteger> qr;
                result.Run(e.reps, [&] { qr = DivideAndRemainder(a, b); });
                answer = HashBigInteger(qr.first) + "/" + HashBigInteger(qr.second);
            }
            result.hash = answer;
            break;
        }
        case Op::Parse: {
            std::string sa = LoadOperand(dir, e.id, 0);
            result.Run(e.reps, [&] {
                BigInteger b = Parse(sa.c_str()); if (b.Limbs().empty()) abort();
            });
            break;
        }
        case Op::ToStr: {
            BigInteger b = LoadInteger(dir, e.id, 0, result.memory);
            result.Run(e.reps, [&] {
                std::string o = ToString(b); if (o.empty()) abort();
            });
            break;
        }
        case Op::DecAdd: case Op::DecMul: case Op::DecDiv: {
            BigDecimal a = BigDecimal::FromString(LoadOperand(dir, e.id, 0));
            BigDecimal b = BigDecimal::FromString(LoadOperand(dir, e.id, 1));
            BigDecimal c;
            if (e.op == Op::DecAdd)
                result.Run(e.reps, [&] { c = a + b; });
            else if (e.op == Op::DecMul)
                result.Run(e.reps, [&] { c = a * b; });
            else
                result.Run(e.reps, [&] { c = a.Divide(b, e.scale, RoundingMode::HALF_EVEN); });
            if (c.Scale() < -1000000) abort();
            break;
        }
        case Op::DecParse: {
            std::string sa = LoadOperand(dir, e.id, 0);
            result.Run(e.reps, [&] {
                BigDecimal b = BigDecimal::FromString(sa); if (b.Scale() < -1000000) abort();
            });
            break;
        }
        case Op::DecToStr: {
            BigDecimal b = BigDecimal::FromString(LoadOperand(dir, e.id, 0));
            result.Run(e.reps, [&] {
                std::string o = b.ToPlainString(); if (o.empty()) abort();
            });
            break;
        }
    }
}

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--thread-mode") {
        std::puts(BIGMATH_USE_THREADS ? "threaded" : "single");
        return 0;
    }
    if (argc == 4 && std::string_view(argv[1]) == "--memory-case")
        return RunMemoryCase(argv[2], argv[3], MeasureCase);
    if (argc == 6 && std::string_view(argv[1]) == "--memory-case") {
        const std::string_view mode(argv[4]), input(argv[5]);
        unsigned setting = 0;
        const auto parsed = std::from_chars(input.data(), input.data() + input.size(), setting);
        if ((mode != "--fermat" && mode != "--wide-blocks") || parsed.ec != std::errc{} ||
            parsed.ptr != input.data() + input.size()) {
            std::fprintf(stderr, "expected --fermat EXPONENT or --wide-blocks COUNT\n");
            return 2;
        }
        return RunMemoryCase(argv[2], argv[3], [&](const auto &dir, const Entry &e, CaseMeasurement &result) {
            if (e.op != Op::Mul) throw std::invalid_argument("Fermat comparisons require multiplication");
            MeasureProduct(dir, e, result, [&](const auto &a, const auto &b) {
                return MultiplyForMemoryTest(a, b, mode, setting);
            });
        });
    }
    if (argc < 8) {
        fprintf(stderr, "usage: %s <dataset_dir> <gmp_csv> <out_csv> <run_id> "
                        "<machine> <profile> <history_md>\n", argv[0]);
        return 2;
    }
    std::string dir = argv[1], gmp_csv = argv[2], out_csv = argv[3];
    int run_id = std::stoi(argv[4]);
    std::string machine = argv[5], profile = argv[6], history_md = argv[7];

    std::vector<Entry> entries = TimingCases(dir);
    std::vector<GmpResult> gmp = ReadGmpResults(gmp_csv, entries);
    if (gmp.size() != entries.size()) {
        fprintf(stderr, "[bigmath] gmp results (%zu) != manifest entries (%zu); "
                        "regenerate the GMP run for this dataset\n",
                gmp.size(), entries.size());
        return 1;
    }

    printf("%-9s %-22s %14s %14s %10s %8s\n",
           "op", "size", "BM warm ms", "GMP warm ms", "BM/GMP", "check");
    printf("%-9s %-22s %14s %14s %10s %8s\n",
           "----", "----", "----------", "------", "------", "-----");
    fflush(stdout);

    std::vector<Row> rows(entries.size());
    int failures = 0;

    for (size_t i = 0; i < entries.size(); ++i) {
        const Entry &e = entries[i];
        Row &row = rows[i];
        row.e = e;
        row.gmp_first_ms = gmp[i].first_ms;
        row.gmp_ms = gmp[i].gmp_ms;
        fprintf(stderr, "  e%-4d %-9s %-22s ...\n", e.id, OpName(e.op), e.Label().c_str());

        CaseMeasurement measurement;
        MeasureCase(dir, e, measurement);
        row.first_ms = measurement.first_ms;
        row.bm_ms = measurement.warm_ms;
        if (HasAnswerHash(e.op))
            row.correct = measurement.hash == gmp[i].hash ? "ok" : "FAIL";
        if (row.correct == "FAIL") ++failures;

        double ratio = row.gmp_ms > 0 ? row.bm_ms / row.gmp_ms : std::numeric_limits<double>::quiet_NaN();
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
        out << "id,op,label,warm_reps,bigmath_first_ms,gmp_first_ms,bigmath_warm_ms,gmp_warm_ms,warm_ratio,correct\n";
        for (const Row &r : rows) {
            double ratio = r.gmp_ms > 0 ? r.bm_ms / r.gmp_ms : std::numeric_limits<double>::quiet_NaN();
            out << r.e.id << ',' << OpName(r.e.op) << ',' << r.e.Label() << ','
                << r.e.reps << ',' << r.first_ms << ',' << r.gmp_first_ms << ',' << r.bm_ms << ',' << r.gmp_ms << ','
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
                 << "(warm BigMath ms / warm GMP ms).\n";
        }
        hist << "\n## run " << run_id << " — machine=" << machine
             << " profile=" << profile << "\n\n"
             << "| op | size | BigMath warm ms | GMP warm ms | BM/GMP | check |\n"
             << "|----|------|-----------:|-------:|-------:|:-----:|\n";
        for (const Row &r : rows) {
            if (!r.e.historical) continue;
            double ratio = r.gmp_ms > 0 ? r.bm_ms / r.gmp_ms : std::numeric_limits<double>::quiet_NaN();
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
