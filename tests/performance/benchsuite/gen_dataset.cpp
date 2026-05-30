// Stage 1: dataset generation.
//
//   gen_dataset <dataset_dir> <profile>
//
// Writes manifest.csv plus one decimal-string file per operand into
// <dataset_dir>. The driver zips the directory afterward. Operands are fully
// determined by (entry id, operand index) so the dataset is reproducible on any
// machine — the GMP and BigMath stages both operate on these exact files.

#include <cstdio>
#include <string>
#include <vector>

#include "bench_common.h"
#include "workload.h"

using namespace bench;

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <dataset_dir> <profile>\n", argv[0]);
        return 2;
    }
    std::string dir = argv[1];
    std::string profile = argv[2];

    std::vector<Entry> entries = BuildWorkload(profile);
    fprintf(stderr, "[gen] profile=%s entries=%zu dir=%s\n",
            profile.c_str(), entries.size(), dir.c_str());

    WriteManifest(dir, entries);

    size_t total_bytes = 0;
    for (const Entry &e : entries) {
        int n = OperandCount(e.op);
        for (int k = 0; k < n; ++k) {
            std::string s = BuildOperandString(e, k);
            std::string path = OperandPath(dir, e.id, k);
            FILE *f = fopen(path.c_str(), "wb");
            if (!f) { fprintf(stderr, "[gen] cannot write %s\n", path.c_str()); return 1; }
            fwrite(s.data(), 1, s.size(), f);
            fclose(f);
            total_bytes += s.size();
        }
        fprintf(stderr, "  e%-4d %-9s %-22s\n", e.id, OpName(e.op), e.Label().c_str());
    }

    fprintf(stderr, "[gen] done: %zu entries, %.1f MB of operands\n",
            entries.size(), total_bytes / 1e6);
    return 0;
}
