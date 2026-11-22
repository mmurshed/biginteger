"""Compare matched dispatch timings from one or more validation artifacts."""

import argparse
import csv
import math
from pathlib import Path
from statistics import geometric_mean, median


def load_samples(directory, variant):
    samples = {}
    for index in range(1, 4):
        with (directory / f"{variant}-{index}.csv").open(newline="") as source:
            rows = {}
            for row in csv.DictReader(source):
                key = tuple(row[field] for field in ("operation", "state", "a_limbs", "b_limbs"))
                value = float(row["median_us"])
                if key in rows or not math.isfinite(value) or value <= 0:
                    raise ValueError(f"Invalid measurement in {directory}: {key} = {value}")
                rows[key] = value
            if not rows or (samples and samples.keys() != rows.keys()):
                raise ValueError(f"Unmatched measurements in {directory}, {variant}-{index}")
            for key, value in rows.items():
                samples.setdefault(key, []).append(value)
    return samples


def compare(directory):
    before = load_samples(directory, "baseline")
    after = load_samples(directory, "candidate")
    if before.keys() != after.keys():
        raise ValueError(f"Baseline and candidate cover different shapes in {directory}")
    groups = {}
    with (directory / "comparison.csv").open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(("operation", "state", "a_limbs", "b_limbs", "baseline_us",
                         "candidate_us", "paired_speedup"))
        for key, baseline in before.items():
            candidate = after[key]
            speedup = median(a / b for a, b in zip(baseline, candidate))
            writer.writerow((*key, median(baseline), median(candidate), f"{speedup:.4f}"))
            groups.setdefault(key[:2], []).append(speedup)
    print(f"### {directory.name}\n")
    print("| Case | State | Geometric mean speedup | Slowest case |")
    print("| --- | --- | ---: | ---: |")
    for (operation, state), ratios in groups.items():
        print(f"| {operation} | {state} | {geometric_mean(ratios):.3f}× | {min(ratios):.3f}× |")
    print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", nargs="+", type=Path)
    for artifact in parser.parse_args().artifacts:
        compare(artifact)
