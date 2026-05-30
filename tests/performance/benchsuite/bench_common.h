// Shared definitions for the split BigMath-vs-GMP benchmark suite.
//
// The suite is split into three independent stages so the slow parts can be
// cached and reused:
//   1. gen_dataset  — materialize random operands once, save to a directory
//                     (the driver then zips it).
//   2. run_gmp      — time GMP on the dataset once per machine, record answer
//                     hashes for correctness cross-check.
//   3. run_bigmath  — time BigMath on the same dataset every run, verify
//                     against the cached GMP hashes, append to history.
//
// All three agree on the workload through the on-disk manifest (manifest.csv),
// which gen_dataset writes and the other two read. The manifest — not this
// header — is the source of truth at run time, so a dataset stays valid even
// if the workload tables below later change.

#ifndef BIGMATH_BENCH_COMMON_H
#define BIGMATH_BENCH_COMMON_H

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace bench {

// ─── operation taxonomy ────────────────────────────────────────────────────
enum class Op {
    Mul,        // BigInteger a * b           (2 operands)
    Div,        // BigInteger a / b           (2 operands)
    Parse,      // decimal string -> BigInteger
    ToStr,      // BigInteger -> decimal string
    DecAdd,     // BigDecimal a + b           (2 operands)
    DecMul,     // BigDecimal a * b           (2 operands)
    DecDiv,     // BigDecimal a / b @ scale   (2 operands)
    DecParse,   // decimal string -> BigDecimal
    DecToStr,   // BigDecimal -> decimal string
};

inline const char *OpName(Op op) {
    switch (op) {
        case Op::Mul:      return "mul";
        case Op::Div:      return "div";
        case Op::Parse:    return "parse";
        case Op::ToStr:    return "tostr";
        case Op::DecAdd:   return "dec_add";
        case Op::DecMul:   return "dec_mul";
        case Op::DecDiv:   return "dec_div";
        case Op::DecParse: return "dec_parse";
        case Op::DecToStr: return "dec_tostr";
    }
    return "?";
}

inline Op OpFromName(const std::string &s) {
    if (s == "mul")       return Op::Mul;
    if (s == "div")       return Op::Div;
    if (s == "parse")     return Op::Parse;
    if (s == "tostr")     return Op::ToStr;
    if (s == "dec_add")   return Op::DecAdd;
    if (s == "dec_mul")   return Op::DecMul;
    if (s == "dec_div")   return Op::DecDiv;
    if (s == "dec_parse") return Op::DecParse;
    if (s == "dec_tostr") return Op::DecToStr;
    throw std::runtime_error("unknown op: " + s);
}

// Number of decimal-string operand files a given op consumes.
inline int OperandCount(Op op) {
    switch (op) {
        case Op::Mul: case Op::Div:
        case Op::DecAdd: case Op::DecMul: case Op::DecDiv:
            return 2;
        default:
            return 1;
    }
}

// Whether the op compares two big integers exactly (so a value hash is worth
// recording). The decimal/float ops compare GMP's mpf against BigDecimal, which
// round differently, so they are timed but not hash-checked.
inline bool HasAnswerHash(Op op) { return op == Op::Mul || op == Op::Div; }

// ─── one benchmark entry ────────────────────────────────────────────────────
struct Entry {
    int  id        = 0;
    Op   op        = Op::Mul;
    int  a_digits  = 0;   // primary operand digit count (int part for Dec ops)
    int  b_digits  = 0;   // secondary operand digits (frac part for 1-operand Dec ops)
    int  scale     = 0;   // DecDiv target scale; 0 elsewhere
    int  reps      = 1;   // timed repetitions; report the minimum
    int  historical = 0;  // 1 => row is mirrored into the long-lived history table

    // Cosmetic label rebuilt from the numeric fields (manifest stays compact).
    std::string Label() const {
        char buf[96];
        switch (op) {
            case Op::DecAdd: case Op::DecMul:
                snprintf(buf, sizeof(buf), "%d.%d digits", a_digits, b_digits);
                break;
            case Op::DecDiv:
                snprintf(buf, sizeof(buf), "%d.%d / %d dp", a_digits, b_digits, scale);
                break;
            default:
                if (b_digits > 0 && OperandCount(op) == 2)
                    snprintf(buf, sizeof(buf), "%dx%d digits", a_digits, b_digits);
                else
                    snprintf(buf, sizeof(buf), "%d digits", a_digits);
        }
        return buf;
    }
};

// ─── repetition policy (shared so GMP and BigMath agree) ────────────────────
inline int IterationsForDigits(int digits) {
    if (digits <= 2000)   return 7;
    if (digits <= 20000)  return 5;
    if (digits <= 100000) return 3;
    return 1;
}

// ─── deterministic operand generation ───────────────────────────────────────
// Random digit string of the given length, no leading zero.
inline std::string GenerateRandomDigits(int length, uint64_t seed) {
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<int> dis(0, 9);
    std::uniform_int_distribution<int> first_dis(1, 9);
    std::string s;
    s.reserve(length);
    s += char('0' + first_dis(gen));
    for (int i = 1; i < length; ++i)
        s += char('0' + dis(gen));
    return s;
}

// Stable per-(entry,operand) seed. Fixed mixing constant => reproducible
// datasets across machines and runs.
inline uint64_t OperandSeed(int id, int operand_index) {
    uint64_t h = 0x9E3779B97F4A7C15ull;
    h ^= (uint64_t)id * 0xD1B54A32D192ED03ull;
    h ^= (uint64_t)(operand_index + 1) * 0xCBF29CE484222325ull;
    h ^= h >> 29; h *= 0xBF58476D1CE4E5B9ull; h ^= h >> 32;
    return h;
}

// Build the decimal-string operand for an entry slot, exactly as the dataset
// generator does. Used only by gen_dataset; the run stages read files instead.
inline std::string BuildOperandString(const Entry &e, int operand_index) {
    uint64_t seed = OperandSeed(e.id, operand_index);
    switch (e.op) {
        case Op::Mul: case Op::Div:
        case Op::Parse: case Op::ToStr:
            return GenerateRandomDigits(operand_index == 0 ? e.a_digits : e.b_digits, seed);

        case Op::DecAdd: case Op::DecMul: case Op::DecDiv: {
            // a_digits = integer-part length, b_digits = fractional-part length.
            std::string ip = GenerateRandomDigits(e.a_digits, seed);
            std::string fp = GenerateRandomDigits(e.b_digits, seed ^ 0xF00Dull);
            return ip + "." + fp;
        }
        case Op::DecParse: case Op::DecToStr: {
            std::string s = GenerateRandomDigits(e.a_digits, seed);
            s.insert(s.size() / 2, ".");
            return s;
        }
    }
    return "";
}

// ─── value hash (base-independent, cheap even at 200M digits) ────────────────
// FNV-1a over the integer's little-endian 32-bit word stream. BigMath limbs are
// base-2^32 values held in 64-bit slots, so the low 32 bits of each limb ARE
// the LE words; GMP produces the identical stream via mpz_export(order=-1,
// size=4, endian=-1). No decimal conversion required, so hashing never touches
// the slow ToString path.
struct Fnv64 {
    uint64_t h = 0xCBF29CE484222325ull;
    void Byte(uint8_t b) { h ^= b; h *= 0x100000001B3ull; }
    void Word32(uint32_t w) {
        Byte((uint8_t)(w)); Byte((uint8_t)(w >> 8));
        Byte((uint8_t)(w >> 16)); Byte((uint8_t)(w >> 24));
    }
    void Bytes(const uint8_t *p, size_t n) { for (size_t i = 0; i < n; ++i) Byte(p[i]); }
};

// ─── manifest I/O ───────────────────────────────────────────────────────────
// CSV columns: id,op,a_digits,b_digits,scale,reps,historical
inline std::string ManifestPath(const std::string &dir) { return dir + "/manifest.csv"; }
inline std::string OperandPath(const std::string &dir, int id, int operand_index) {
    return dir + "/e" + std::to_string(id) + "_" + std::to_string(operand_index) + ".txt";
}

inline void WriteManifest(const std::string &dir, const std::vector<Entry> &entries) {
    std::ofstream out(ManifestPath(dir));
    if (!out) throw std::runtime_error("cannot write manifest in " + dir);
    out << "id,op,a_digits,b_digits,scale,reps,historical\n";
    for (const Entry &e : entries) {
        out << e.id << ',' << OpName(e.op) << ',' << e.a_digits << ',' << e.b_digits
            << ',' << e.scale << ',' << e.reps << ',' << e.historical << '\n';
    }
}

inline std::vector<Entry> ReadManifest(const std::string &dir) {
    std::ifstream in(ManifestPath(dir));
    if (!in) throw std::runtime_error("cannot read manifest in " + dir +
                                      " (run gen_dataset first)");
    std::vector<Entry> entries;
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> f;
        size_t start = 0;
        for (size_t i = 0; i <= line.size(); ++i) {
            if (i == line.size() || line[i] == ',') {
                f.push_back(line.substr(start, i - start));
                start = i + 1;
            }
        }
        if (f.size() < 7) throw std::runtime_error("bad manifest line: " + line);
        Entry e;
        e.id         = std::stoi(f[0]);
        e.op         = OpFromName(f[1]);
        e.a_digits   = std::stoi(f[2]);
        e.b_digits   = std::stoi(f[3]);
        e.scale      = std::stoi(f[4]);
        e.reps       = std::stoi(f[5]);
        e.historical = std::stoi(f[6]);
        entries.push_back(e);
    }
    return entries;
}

inline std::string LoadOperand(const std::string &dir, int id, int operand_index) {
    std::ifstream in(OperandPath(dir, id, operand_index), std::ios::binary);
    if (!in) throw std::runtime_error("missing operand file: " +
                                      OperandPath(dir, id, operand_index));
    std::string s((std::istreambuf_iterator<char>(in)),
                  std::istreambuf_iterator<char>());
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return s;
}

// ─── per-machine result CSV I/O ─────────────────────────────────────────────
// GMP results: id,op,label,reps,gmp_ms,hash
struct GmpResult { double gmp_ms = 0; std::string hash; };

inline void WriteGmpResults(const std::string &path,
                            const std::vector<Entry> &entries,
                            const std::vector<GmpResult> &res,
                            const std::string &header_comment) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write gmp results: " + path);
    out << "# " << header_comment << "\n";
    out << "id,op,label,reps,gmp_ms,hash\n";
    for (size_t i = 0; i < entries.size(); ++i) {
        const Entry &e = entries[i];
        out << e.id << ',' << OpName(e.op) << ',' << e.Label() << ',' << e.reps
            << ',' << res[i].gmp_ms << ',' << res[i].hash << '\n';
    }
}

inline std::vector<GmpResult> ReadGmpResults(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read gmp results: " + path +
                                      " (run run_gmp first)");
    std::vector<GmpResult> res;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("id,", 0) == 0) continue; // header
        // split (label may contain spaces but never commas)
        std::vector<std::string> f;
        size_t start = 0;
        for (size_t i = 0; i <= line.size(); ++i) {
            if (i == line.size() || line[i] == ',') {
                f.push_back(line.substr(start, i - start));
                start = i + 1;
            }
        }
        if (f.size() < 6) continue;
        GmpResult r;
        r.gmp_ms = std::stod(f[4]);
        r.hash   = f[5];
        res.push_back(r);
    }
    return res;
}

// ─── timing ─────────────────────────────────────────────────────────────────
template <typename F>
inline double BestMs(F &&fn, int reps) {
    double best = 1e300;
    for (int i = 0; i < reps; ++i) {
        auto t0 = std::chrono::high_resolution_clock::now();
        fn();
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        if (ms < best) best = ms;
    }
    return best;
}

} // namespace bench

#endif // BIGMATH_BENCH_COMMON_H
