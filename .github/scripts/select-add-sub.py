"""Select conservative add/sub cutoffs from each platform's three surveys."""

import argparse
from pathlib import Path
import re


MACROS = (
    "BIGMATH_ADD_PARALLEL_WARM_LIMBS",
    "BIGMATH_ADD_PARALLEL_PARKED_LIMBS",
    "BIGMATH_SUB_PARALLEL_WARM_LIMBS",
    "BIGMATH_SUB_PARALLEL_PARKED_LIMBS",
)


def read_header(path):
    text = path.read_text()
    values = {}
    for macro in MACROS:
        matches = re.findall(rf"^#define {macro} ([0-9]+)$", text, re.MULTILINE)
        if len(matches) != 1:
            raise ValueError(f"Expected one definition of {macro} in {path}")
        values[macro] = int(matches[0])
    return values


def select(directory):
    samples = [read_header(directory / f"survey-{index}.h") for index in range(1, 4)]
    cutoffs = {}
    for macro in MACROS:
        choices = [sample[macro] for sample in samples]
        # Each individual survey establishes a sustained speedup from its
        # cutoff. Choose the largest cutoff that all three surveys support.
        cutoffs[macro] = max(choices) if all(choices) else 0
    return cutoffs


def write_header(path, cutoffs):
    lines = ["#pragma once", "// Selected from three matching platform surveys.",
             "#if BIGMATH_LIMB_64 == 1"]
    for macro, value in cutoffs.items():
        lines.extend((f"#ifndef {macro}", f"#define {macro} {value}", "#endif"))
    lines.append("#endif")
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("artifacts", nargs="+", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    print("### Addition and subtraction cutoffs\n")
    print("Limbs per operand. Zero keeps that worker state serial.\n")
    print("| Platform | Add warm | Add parked | Subtract warm | Subtract parked |")
    print("| --- | ---: | ---: | ---: | ---: |")
    for artifact in args.artifacts:
        platform = artifact.name.removeprefix("survey-")
        selected = select(artifact)
        write_header(args.output / f"{platform}.h", selected)
        print(f"| {platform} | " + " | ".join(str(n) for n in selected.values()) + " |")
