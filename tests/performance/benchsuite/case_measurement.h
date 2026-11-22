// First and warm calls in a single benchmark process. With --memory-case the
// driver execs a fresh process for each row, so the OS peak belongs to one case.
#pragma once
#include "bench_common.h"

#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#endif

namespace bench {
inline uint64_t PeakResidentBytes() {
#if defined(__APPLE__) || defined(__linux__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0)
        throw std::runtime_error("getrusage failed");
#ifdef __APPLE__
    return static_cast<uint64_t>(usage.ru_maxrss);
#else
    return static_cast<uint64_t>(usage.ru_maxrss) * 1024;
#endif
#else
    throw std::runtime_error("peak RSS measurement requires Linux or macOS");
#endif
}

struct CaseMeasurement {
    bool memory = false;
    double first_ms = 0, warm_ms = 0;
    uint64_t setup_peak = 0, first_peak = 0, warm_peak = 0;
    std::string hash = "-";

    template <class F>
    void Run(int reps, F operation) {
        if (reps < 1) throw std::invalid_argument("expected at least one warm call");
        if (memory) setup_peak = PeakResidentBytes();
        auto timed = [&] {
            const auto start = std::chrono::steady_clock::now();
            operation();
            return std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start).count();
        };
        first_ms = timed();
        if (memory) first_peak = PeakResidentBytes();
        warm_ms = timed();
        for (int i = 1; i < reps; ++i) warm_ms = std::min(warm_ms, timed());
        if (memory) warm_peak = PeakResidentBytes();
    }
};

inline std::string RawOperandPath(const std::string &dir, int id, int slot) {
    return dir + "/e" + std::to_string(id) + "_" + std::to_string(slot) + ".raw";
}

template <class Measure>
int RunMemoryCase(const std::string &dir, const std::string &index, Measure measure) {
    const Entry entry = TimingCases(dir).at(std::stoul(index));
    CaseMeasurement result;
    result.memory = true;
    try {
        measure(dir, entry, result);
    } catch (const std::invalid_argument &e) {
        // Only this explicit library limit means "unsupported". All other
        // failures must fail the benchmark rather than turn into missing data.
        if (entry.op != Op::Mul ||
            std::string(e.what()).find("exceed the CRT NTT length ceiling") == std::string::npos)
            throw;
        std::printf("crt_length_limit,,,,,,\n");
        return 0;
    }
    std::printf("ok,%.9g,%.9g,%llu,%llu,%llu,%s\n", result.first_ms, result.warm_ms,
                (unsigned long long)result.setup_peak, (unsigned long long)result.first_peak,
                (unsigned long long)result.warm_peak, result.hash.c_str());
    return 0;
}
} // namespace bench
