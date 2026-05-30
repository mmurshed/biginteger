// The benchmark workload tables. Only gen_dataset includes this; once a dataset
// is materialized the manifest is the source of truth for the run stages.
//
// Three profiles trade dataset size / wall time for coverage:
//   quick   — seconds, tiny dataset; for testing the pipeline itself
//   default — minutes, modest dataset; the everyday regression set
//   full    — the historical 200M-digit sweep from bench_vs_gmp.cpp (gigabytes)

#ifndef BIGMATH_BENCH_WORKLOAD_H
#define BIGMATH_BENCH_WORKLOAD_H

#include <string>
#include <vector>

#include "bench_common.h"

namespace bench {

inline int RepsFor(Op op, int a_digits, int b_digits) {
    int d = (op == Op::DecAdd || op == Op::DecMul || op == Op::DecDiv ||
             op == Op::DecParse || op == Op::DecToStr)
                ? a_digits + b_digits
                : std::max(a_digits, b_digits);
    int reps = IterationsForDigits(d);
    if ((op == Op::ToStr || op == Op::DecToStr) && d >= 100000) reps = 1;
    return reps;
}

struct Tables {
    std::vector<std::pair<int, int>> mul_bal, mul_skew, div_bal, div_skew;
    std::vector<int> parse, tostr;
    std::vector<std::pair<int, int>> dec_simple; // (int,frac) for add/mul/div
    std::vector<std::tuple<int, int, int>> dec_div_scales; // (int,frac,scale)
    std::vector<int> dec_parse, dec_tostr;
    int hist_mul_min = 0;  // mark balanced mul >= this as historical
    int hist_div_min = 0;  // mark skewed   div >= this as historical
};

inline Tables TablesFor(const std::string &profile) {
    Tables t;
    if (profile == "quick") {
        t.mul_bal  = {{1000,1000},{10000,10000}};
        t.mul_skew = {{10000,1000}};
        t.div_bal  = {{1000,1000},{10000,10000}};
        t.div_skew = {{40000,10000}};
        t.parse    = {1000,10000};
        t.tostr    = {1000,10000};
        t.dec_simple    = {{100,10},{1000,100}};
        t.dec_div_scales = {{2000,200,0},{2000,200,100}};
        t.dec_parse = {100,1000};
        t.dec_tostr = {100,1000};
        t.hist_mul_min = 10000;
        t.hist_div_min = 40000;
    } else if (profile == "full") {
        for (int d : {1000,5000,10000,50000,100000,500000,1000000,2000000,
                      5000000,10000000,20000000,50000000,100000000,200000000})
            t.mul_bal.push_back({d, d});
        t.mul_skew = {{100000,10000},{500000,50000},{1000000,100000},
                      {2000000,200000},{5000000,500000},{10000000,1000000},
                      {20000000,2000000},{50000000,5000000},{100000000,10000000},
                      {200000000,20000000}};
        for (int d : {1000,5000,10000,50000,100000,500000,1000000,5000000})
            t.div_bal.push_back({d, d});
        t.div_skew = {{40000,10000},{100000,10000},{200000,50000},{500000,100000},
                      {1000000,200000},{2000000,500000},{5000000,1000000},
                      {10000000,2000000},{20000000,4000000},{50000000,10000000},
                      {100000000,20000000},{200000000,40000000}};
        for (int d : {1000,10000,50000,100000,500000,1000000,2000000,5000000,
                      10000000,20000000,50000000})
            t.parse.push_back(d);
        for (int d : {1000,10000,50000,100000,200000,500000,1000000,2000000,
                      5000000,10000000,20000000})
            t.tostr.push_back(d);
        for (int d : {100,1000,5000,20000}) t.dec_simple.push_back({d, d/10});
        for (int s : {0,10,100,1000,5000}) t.dec_div_scales.push_back({2000,200,s});
        t.dec_parse = {100,1000,10000,50000};
        t.dec_tostr = {100,1000,10000,50000};
        t.hist_mul_min = 10000000;
        t.hist_div_min = 10000000;
    } else { // default
        for (int d : {1000,10000,100000,1000000}) t.mul_bal.push_back({d, d});
        t.mul_skew = {{100000,10000},{1000000,100000}};
        for (int d : {1000,10000,100000,1000000}) t.div_bal.push_back({d, d});
        t.div_skew = {{100000,10000},{1000000,200000}};
        t.parse = {1000,10000,100000,1000000};
        t.tostr = {1000,10000,100000,1000000};
        for (int d : {100,1000,5000,20000}) t.dec_simple.push_back({d, d/10});
        for (int s : {0,100,1000}) t.dec_div_scales.push_back({2000,200,s});
        t.dec_parse = {100,1000,10000};
        t.dec_tostr = {100,1000,10000};
        t.hist_mul_min = 1000000;
        t.hist_div_min = 1000000;
    }
    return t;
}

inline std::vector<Entry> BuildWorkload(const std::string &profile) {
    Tables t = TablesFor(profile);
    std::vector<Entry> out;
    int id = 0;
    auto add = [&](Op op, int a, int b, int scale, int hist) {
        Entry e;
        e.id = id++;
        e.op = op;
        e.a_digits = a;
        e.b_digits = b;
        e.scale = scale;
        e.reps = RepsFor(op, a, b);
        e.historical = hist;
        out.push_back(e);
    };

    for (auto [a, b] : t.mul_bal)  add(Op::Mul, a, b, 0, a >= t.hist_mul_min ? 1 : 0);
    for (auto [a, b] : t.mul_skew) add(Op::Mul, a, b, 0, 0);
    for (auto [a, b] : t.div_bal)  add(Op::Div, a, b, 0, 0);
    for (auto [a, b] : t.div_skew) add(Op::Div, a, b, 0, a >= t.hist_div_min ? 1 : 0);
    for (int d : t.parse) add(Op::Parse, d, 0, 0, 0);
    for (int d : t.tostr) add(Op::ToStr, d, 0, 0, 0);
    for (auto [a, b] : t.dec_simple) { add(Op::DecAdd, a, b, 0, 0); add(Op::DecMul, a, b, 0, 0); }
    for (auto [a, b] : t.dec_simple) add(Op::DecDiv, a, b, b, 0);
    for (auto [a, b, s] : t.dec_div_scales) add(Op::DecDiv, a, b, s, 0);
    for (int d : t.dec_parse) add(Op::DecParse, d, 0, 0, 0);
    for (int d : t.dec_tostr) add(Op::DecToStr, d, 0, 0, 0);

    return out;
}

} // namespace bench

#endif // BIGMATH_BENCH_WORKLOAD_H
