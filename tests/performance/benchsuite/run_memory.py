#!/usr/bin/env python3
"""Measure each operation/library in a fresh process, on matched operands."""
import argparse
import csv
import math
from pathlib import Path
import random
import subprocess
import tempfile

METRICS = ("status", "first_ms", "warm_ms", "setup_peak_bytes", "first_peak_bytes",
           "warm_peak_bytes", "hash")
FIELDS = ("id", "op", "a_digits", "b_digits", "scale", "warm_reps", "source",
          "library", *METRICS, "check")
# Bracket transform doublings within the existing sweep, then test the ceiling.
EXTRA_DIGITS = (10_200_000, 20_300_000, 40_000_000, 41_000_000, 80_000_000,
                82_000_000, 160_000_000, 162_000_000, 320_000_000, 324_000_000,
                500_000_000, 750_000_000, 1_000_000_000,
                480_000_000, 484_000_000, 486_000_000, 520_000_000,
                2_500_000_000, 5_000_000_000)
INTEGER_INPUT = {"mul", "div", "div_q", "div_qr", "tostr"}
CHECKED = {"mul", "div_q", "div_qr"}


def generate_raw(path, digits, slot):
    # A value with floor(d*log2(10)) bits lies strictly between 10^(d-1)
    # and 10^d. Write the GMP raw header and a seeded big-endian magnitude.
    bits = math.floor(digits * math.log2(10))
    size = (bits + 7) // 8
    if size >= 1 << 31:
        raise ValueError("operand exceeds the signed GMP raw length field")
    high_bits = (bits - 1) % 8 + 1
    rng = random.Random(2 * digits + slot)
    with path.open("wb") as out:
        out.write(size.to_bytes(4, "big"))
        out.write(bytes([(1 << (high_bits - 1)) | rng.getrandbits(high_bits - 1)]))
        left = size - 1
        while left:
            count = min(left, 1 << 20)
            out.write(rng.randbytes(count))
            left -= count


def measure(dataset, scratch, output, builds, extended, balanced_digits, product_options, warm_reps):
    with (dataset / "manifest.csv").open() as f:
        reader = csv.DictReader(f)
        manifest_fields = reader.fieldnames
        entries = list(reader)
    extra_ids = set()
    if extended or balanced_digits:
        next_id = max((int(e["id"]) for e in entries), default=-1) + 1
        if balanced_digits:
            entries = []
        for idx, digits in enumerate(balanced_digits or EXTRA_DIGITS, next_id):
            extra_ids.add(str(idx))
            entries.append(dict(id=str(idx), op="mul", a_digits=str(digits),
                                b_digits=str(digits), scale="0", reps=str(warm_reps), historical="0"))
    with (scratch / "manifest.csv").open("w") as f:
        writer = csv.DictWriter(f, manifest_fields)
        writer.writeheader()
        writer.writerows(entries)

    writer = csv.DictWriter(output, FIELDS, lineterminator="\n")
    writer.writeheader()
    index = 0
    for entry in entries:
        eid, op = entry["id"], entry["op"]
        is_extra = eid in extra_ids
        count = 2 if op in {"mul", "div", "div_q", "div_qr", "dec_add", "dec_mul", "dec_div"} else 1
        fixtures = []
        for slot in range(count):
            stem = f"e{eid}_{slot}"
            if op in INTEGER_INPUT:
                path = scratch / (stem + ".raw")
                if is_extra:
                    generate_raw(path, int(entry["a_digits"]), slot)
                else:
                    subprocess.run([builds["gmp"], "--export-integer", dataset, eid,
                                    str(slot), path], check=True)
            else:
                path = scratch / (stem + ".txt")
                path.symlink_to(dataset / path.name)
            fixtures.append(path)

        modes = ("div_q", "div_qr") if op in {"div", "div_q"} else (op,)
        for mode in modes:
            print(f"  {eid} {mode} {entry['a_digits']} x {entry['b_digits']}", flush=True)
            ref = None
            for library, executable in builds.items():
                command = [executable, "--memory-case", scratch, str(index),
                           *product_options.get(library, ())]
                result = subprocess.check_output(command, text=True).strip().split(",")
                if len(result) != len(METRICS) or result[0] not in {"ok", "crt_length_limit"}:
                    raise RuntimeError(f"invalid {library} result: {result}")
                metrics = dict(zip(METRICS, result))
                if metrics["status"] == "ok":
                    setup = int(metrics["setup_peak_bytes"])
                    first = int(metrics["first_peak_bytes"])
                    warm = int(metrics["warm_peak_bytes"])
                    if not 0 < setup <= first <= warm:
                        raise RuntimeError(f"invalid RSS high-water marks: {metrics}")
                if library == "gmp":
                    if metrics["status"] != "ok":
                        raise RuntimeError("missing GMP reference")
                    ref = metrics["hash"]
                check = "n/a"
                if mode in CHECKED and metrics["status"] == "ok":
                    if ref == "-" or metrics["hash"] != ref:
                        raise RuntimeError(f"hash mismatch: case {eid} {mode}, {library}")
                    check = "reference" if library == "gmp" else "ok"
                writer.writerow(dict(id=eid, op=mode, a_digits=entry["a_digits"],
                                     b_digits=entry["b_digits"], scale=entry["scale"],
                                     warm_reps=entry["reps"], source="seeded_binary" if is_extra else "manifest",
                                     library=library, check=check, **metrics))
                output.flush()
            index += 1
        for path in fixtures:
            path.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    parser.add_argument("build", type=Path, help="default Release build with both runners")
    parser.add_argument("single_build", type=Path, help="Release build with BIGMATH_USE_THREADS=0")
    parser.add_argument("output", type=Path)
    cases = parser.add_mutually_exclusive_group()
    cases.add_argument("--extended", action="store_true", help="include balanced products through 5B digits")
    cases.add_argument("--balanced-digits", type=int, nargs="+", metavar="DIGITS",
                       help="measure only seeded balanced products at these decimal lengths")
    parser.add_argument("--balanced-warm-reps", type=int, metavar="N",
                        help="record the fastest of N warm calls for each --balanced-digits case (default 1)")
    parser.add_argument("--fermat-blocks", type=int,
                        help="compare the Fermat FFT with this many wide CRT blocks per operand")
    parser.add_argument("--fermat-exponents", type=int, nargs="+",
                        help="sweep FFT exponents for the Fermat path (0 selects automatically)")
    parser.add_argument("--fermat-threads", choices=("both", "default", "single"), default="both",
                        help="build modes to measure in the Fermat comparison")
    args = parser.parse_args()
    if args.fermat_blocks is not None and (args.fermat_blocks < 2 or not args.balanced_digits):
        parser.error("--fermat-blocks requires at least 2 blocks and --balanced-digits")
    if args.fermat_blocks is None and (args.fermat_exponents is not None or args.fermat_threads != "both"):
        parser.error("Fermat options require --fermat-blocks")
    args.fermat_exponents = args.fermat_exponents or [0]
    if any(k != 0 and not 6 <= k <= 14 for k in args.fermat_exponents):
        parser.error("Fermat exponents must be 0 (automatic) or 6 through 14")
    if len(set(args.fermat_exponents)) != len(args.fermat_exponents):
        parser.error("Fermat exponents must be distinct")
    if args.balanced_digits and any(d < 1 for d in args.balanced_digits):
        parser.error("--balanced-digits requires positive decimal lengths")
    if args.balanced_warm_reps is not None and (args.balanced_warm_reps < 1 or not args.balanced_digits):
        parser.error("--balanced-warm-reps requires --balanced-digits and a positive count")
    builds = {"gmp": args.build.resolve() / "bench_run_gmp",
              "bigmath_default": args.build.resolve() / "bench_run_bigmath",
              "bigmath_single": args.single_build.resolve() / "bench_run_bigmath"}
    product_options = {}
    for label, expected in [("default", "threaded"), ("single", "single")]:
        library = f"bigmath_{label}"
        executable = builds[library]
        if args.fermat_blocks is not None:
            del builds[library]
            if args.fermat_threads not in ("both", label):
                continue
            block_label = f"wide_blocks{args.fermat_blocks}_{label}"
            builds[block_label] = executable
            product_options[block_label] = ("--wide-blocks", str(args.fermat_blocks))
            for exponent in args.fermat_exponents:
                choice = f"k{exponent}" if exponent else "auto"
                fft_label = f"fermat_{label}_{choice}"
                builds[fft_label] = executable
                product_options[fft_label] = ("--fermat", str(exponent))
        mode = subprocess.check_output([executable, "--thread-mode"], text=True).strip()
        if mode != expected:
            parser.error(f"{library} requires a {expected} build, got {mode}")
    # Keep existing results intact if a child fails halfway through the sweep.
    with tempfile.TemporaryDirectory(dir=args.output.parent) as tmp:
        scratch = Path(tmp).resolve()
        partial = scratch / "results.csv"
        with partial.open("w", newline="") as out:
            measure(args.dataset.resolve(), scratch, out, builds, args.extended, args.balanced_digits,
                    product_options, args.balanced_warm_reps or 1)
        partial.replace(args.output)


if __name__ == "__main__":
    main()
