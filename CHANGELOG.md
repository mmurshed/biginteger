# Changelog

## Unreleased — Toom-3 tier in the Karatsuba recursion

Below the NTT entry, products and squares used plain Karatsuba. Inside the
pointer-based Karatsuba recursion (64-bit limbs), a near-balanced level whose
shorter operand has ≥ 168 limbs (`BIGMATH_MUL_TOOM33_THRESHOLD`) now splits
three ways (Toom-3: points 0, ±1, 2, ∞, Bodrato interpolation on unsigned
limbs); squares from 160 limbs (`BIGMATH_SQR_TOOM3_THRESHOLD`). This also
reaches the balanced pieces of sliced unbalanced products and the products
inside Burnikel-Ziegler division. Through the public API on an M1 Max:

| case | change |
|---|---:|
| product 192×192 / 288×288 | −5.9% / −5.6% |
| square 192 / 256 | −4.0% / −2.9% |
| division 800 ÷ 400, 1600 ÷ 400 limbs | −2.6%, −2.5% |
| product 256×256 | +1 to +2.5% |

A Toom-4 split (points 0, ±1, ±2, ½, ∞) is implemented and unit-tested but
off (`BIGMATH_MUL_TOOM44_THRESHOLD = 1000000`): on M1 it has no window below
the NTT entry. `tests/performance/toom_threshold_sweep.cpp` retunes both.

## Unreleased — small-number storage

`BigInteger` keeps values of up to 4 limbs (256 bits) inside the object
instead of in a heap-allocated `std::vector`, and the operators work on
operands of up to 16 limbs in a stack buffer. Through the public API on an
M1 Max, main → this change, ns per operation:

| limbs | `a*b` | `a+b` | `a-b` | `a/b` | `a%b` | copy | `a + 5` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 20.6 → 7.5 | 22.8 → 9.7 | 24.0 → 10.3 | 33.6 → 15.1 | 47.3 → 13.4 | 13.2 → 0.3 | 66.2 → 9.9 |
| 2 | 30.7 → 10.6 | 23.6 → 10.2 | 26.1 → 11.2 | 21.6 → 7.3 | 47.5 → 9.3 | 13.5 → 0.3 | 68.8 → 13.6 |
| 4 | 40.8 → 29.6 | 25.1 → 10.9 | 26.2 → 12.1 | 46.1 → 30.6 | 67.7 → 33.4 | 14.9 → 0.3 | 70.5 → 14.0 |
| 6 | 56.4 → 46.7 | 26.5 → 24.4 | 28.0 → 25.6 | 66.8 → 31.5 | 71.6 → 46.3 | 15.9 → 15.9 | 71.9 → 26.2 |
| 16 | 168.9 → 152.2 | 29.0 → 29.7 | 30.1 → 30.2 | 66.4 → 45.1 | 90.0 → 64.0 | 15.7 → 15.7 | 70.2 → 27.8 |

`BigIntegerBuilder::From(integer)` builds the limbs directly instead of
formatting and re-parsing the number, which is most of the `a + 5` gain.

### API change

- `BigInteger::GetInteger()` now returns `std::vector<DataT>` **by value**
  (a small value has no vector to refer to). Code that indexes it, iterates
  it or passes it to a function still compiles, but each call copies. Use
  the new `Limbs()` (`std::span<const DataT>`, no copy) to read limbs, or pass
  `x.AsVector()` to functions taking `std::vector<DataT> const &` (no copy for
  values on the heap). Do not bind `x.AsVector().get()` to a reference: for
  a small value it refers to a temporary (clang warns, via
  `[[clang::lifetimebound]]`).
- New: `BigInteger::FromLimbs(span, negative)`, `IsInline()`,
  `InlineLimbs`. Moved-from values are zero.

### Costs

- `sizeof(BigInteger)` grows from 32 to 64 bytes.
- Constructing a small `BigInteger` from a vector copies the limbs inline
  (+~1 ns). On one shape this shows: dividing two 1000-limb numbers whose
  quotient is 1–2 limbs goes from 66 to 72 ns in a steady loop (up to ~90 ns
  right after a run of large operations). The division itself is unchanged;
  the quotient's vector is now freed at construction rather than with the
  `BigInteger`, which interacts with allocator state.

## Unreleased — 2026-11-21 division pass (PRs #141–#149)

The same treatment for division: measure the whole (divisor length, ratio)
grid against GMP, then rebuild each layer. Quotient and remainder through
the public API, dividend 5× the divisor, before → after the pass: 40k digits
0.88 → 0.28 ms, 200k 1.67 → 1.18 ms, 1M 5.3 → 4.4 ms, 5M 39 → 25 ms, 20M
141 → 104 ms. Division now **beats GMP from 200k-digit dividends** (was
500k), at 0.19–0.67× on the canonical skewed rows, and trails by 1.4–2.0×
below that (was 2.0–2.9×). ToString, which is built on division, beats GMP
cold from 200k digits (was 1M). Tables in [BENCHMARK.md](BENCHMARK.md),
write-ups in [docs/DIVISION.md](docs/DIVISION.md).

### Performance

- **Schoolbook basecase** (#142): `FastDivision` on 64-bit limbs is a
  raw-buffer Knuth D loop on a flag-chain `SubMul1` kernel, with stack
  workspace up to 192 limbs. 2.2–3.5× of GMP's basecase → 1.1–1.6×.
- **Burnikel–Ziegler in place** (#143): rewritten in the shape of GMP's
  `mpn_dcpi1_div_qr` — it works inside the shifted dividend buffer, splits
  the divisor by pointer offset (no zero-padding for odd lengths) and bottoms
  out in the new basecase at 64 limbs. 1.5–2.4× faster; a quotient shorter
  than the divisor now costs what quotient-sized division costs.
- **Newton reciprocal** (#144): a guarded precision schedule (exact base of
  ≤ 32 limbs, each earlier precision `ceil((p+2)/2)`) with an unpadded Newton
  step replaces exact doubling, whose error had drifted to ~B^5.5 at 1024
  limbs and was being repaired by an extra full-size iteration. Contract
  `0 ≤ ⌊B^2n/D⌋ − R ≤ 3`, no refinement pass.
- **Dispatch re-sweep** (#145): one rule per regime instead of seven ratio
  bands. BZ from 96-limb divisors (was 512, with shape conditions), Newton
  from 1536 limbs at ratio ≥ 2 (768 at ratio ≥ 3), quotient-sized division
  when the quotient is shorter than the divisor but ≥ 3072 limbs. The old
  bands were more than 10% off the best algorithm on 69 of 232 grid cells,
  by up to 2.8×; 3 cells remain, all under 16%.
- **Single-limb divisor** (#146): `ClassicDivision::DivRem1` keeps the
  running remainder normalized across the loop. 6.4 → 4.9 ns per limb.
- **Quotient-only division** (#147): `Divide` / `operator/` with a quotient
  shorter than the divisor computes it from the operands' top limbs — a
  `(2δ+4)`-by-`(δ+4)` division instead of work proportional to `δ · nb`,
  with one full-width product only when the top remainder cannot decide.
  1024-limb divisor, 21-limb quotient: 12.4 → 0.58 µs.
- **Newton leading block** (#149): a leading block of ≤ 32 quotient limbs
  goes to schoolbook instead of a full-size Newton chunk. The normalization
  shift pushes one limb out of a `k·n`-limb dividend for one random input in
  three; those 8192 ÷ 4096-limb divisions are 16% faster.
- **Thread pool** (#141, #148): small NTT products (≈ 310–512 limbs per
  operand) are taken only while the pool is warm — isolated ones were up to
  2.2× slower than Karatsuba because parked workers take as long to wake as
  the whole multiply — and awake workers now pick a job up without a lock.
  Every worker used to take the pool mutex on every job; the convoy cost
  ~13 µs per warm dispatch (8 slices × 2 µs: 15.0 → 4.0 µs per job) and made
  NTT timings depend on unrelated code layout: the same multiply ran 15–25%
  slower in one build than in another.

### Tests

New unit suites for each change: `test_fast_division_basecase`,
`test_bz_in_place`, `test_newton_division` (reciprocal contract, zero
fallbacks, short leading block), `test_division_dispatch` (both sides of
every dispatch edge), `test_div_single_limb`, `test_division_quotient_only`
(`q·b + r` with `r = 0, 1, b − 1` and `q·b − 1`), and two pool tests (job
descriptions never mix; every worker answers every dispatch). The dispatcher
was also swept against GMP on 4612 shapes with no mismatch and no Newton
fallback. ASan + UBSan on the division suites, ThreadSanitizer on the pool,
NTT and Newton suites.

### Not done / follow-ups

- Divisors of a few hundred limbs at high ratio are still ~1.8–2.0× behind
  GMP (100k ÷ 10k digits). Both BZ and Newton land there; closing it needs a
  faster multiply below the NTT entry (GMP has Toom tiers and an assembly
  basecase), not another division algorithm.
- Quotient-only division with a short quotient is ~1.8× of GMP's
  `mpz_tdiv_q`; the top division could use an approximate-quotient basecase.
- Cold pool dispatch is unchanged (~45 µs until the last parked worker
  arrives). Per-worker condition variables were measured and gave nothing:
  the stagger is kernel wake latency.
- `bench_vs_gmp` times `a / b` against `mpz_tdiv_qr`. On the skewed rows both
  sides still do a full division, but the balanced rows are no longer
  comparable; the harness should time like against like.
- Ideas not pursued: a half-size ("mu") reciprocal for ratio ≈ 2; cyclic
  residues plus a low-limb CRT to halve Newton's transforms at power-of-two
  sizes; a 2/1 reciprocal without `__udivti3`.

## Unreleased — 2026-11-20 multiplication pass (PRs #129–#139)

A top-to-bottom pass over multiplication. Through the public API, old `main`
vs this tree on the same machine: multiplication is **1.5–4.5× faster at
every size** and now beats GMP from ~10k decimal digits (was ~500k);
division, parse and ToString — untouched, but built on `Multiply` — are
1.3–2.8× faster from 100k–200k digits. Tables in
[docs/MULTIPLICATION.md](docs/MULTIPLICATION.md) and
[BENCHMARK.md](BENCHMARK.md).

### Correctness

- **Karatsuba heap overflow on skewed operands** (#137). The split clamp made
  the recursion `la/lb` levels deep when one operand was ≥ 2× the other, and
  from about `la > 16·lb` the parked temporaries ran past the `16n`
  workspace — on shapes the dispatcher routed there (e.g. 1214 × 65 limbs).
  Such operands are now sliced into balanced pieces.
- **Thread pool unsafe under concurrent callers** (#133). A second user
  thread dispatching while a job was in flight overwrote it; two threads
  doing NTT-sized multiplies hung. One caller now owns the pool's slot and
  any other runs its work inline.
- Linux build of the unit tests restored (#130; ambiguous
  `BigIntegerBuilder::From` overload since #128).

### Performance

- **Limb kernels** (#136): `common/LimbKernels.h` — `add_n`/`sub_n`/`mul_1`/
  `addmul_1` in 4-limb straight-line blocks on the carry builtins, with the
  low and high product halves on separate carry chains. `addmul_1` 2.5–3.3 →
  1.64–1.76 cycles/limb, within ~1.2× of GMP's assembly, no assembly here.
  The documented "portable-C++ wall" below the NTT band was the loop idiom.
  Karatsuba leaf cutoff 32 → 24.
- **CRT NTT kernels**: NEON tail kernels so every transform size class is
  all-NEON (#131, two thirds of sizes ran 3–4 scalar layers); leaner Shoup
  butterfly — `vqdmulhq_s32` quotient and `vminq_u32` reduce, ~11 ops instead
  of ~19 (#132).
- **Threading**: pool rewrite — claim-based dispatch and a 100 µs spin before
  parking, `ParallelDo(6)` 22.5 → 2.9 µs (#133); Garner/carry finalize with
  NEON digits and chunk-parallel carry, 27 → 0.9 ms at n = 2^22 (#134);
  pack inside the forward units and a six-unit split inverse (#135).
- **Truncated transform length** (#138): products just above a power of two
  `H` use an `H + L` transform (mod `x^H + 1` and mod `x^L − 1`, recombined)
  instead of `2H` — 1.3–2.2× in that band, removing most of the power-of-two
  staircase below the MFA band.
- **Dispatch** (#129, #139): product vector moved instead of copied into the
  result `BigInteger`; NTT entry from an area rule (`la·lb ≥ 96000`,
  `min ≥ 64`, total ≥ 640 — balanced from ~310 limbs per operand instead of
  640); schoolbook only for `min ≤ 24`; square NTT threshold 640 → 352.

### Not done / follow-ups

- Division below ~200k-digit dividends is unchanged (`FastDivision`'s
  `submul` loops still use the old carry idiom; the new kernels apply
  directly), and the Newton/BZ dispatch frontier has not been re-swept
  against the faster multiply. *(Done in the 2026-11-21 division pass,
  above.)*
- 1–4-limb products are still 3–4× behind GMP: the result allocation
  dominates.
- The fused-MFA band (n ≥ 2^22) keeps its power-of-two staircase.
- No sanitizer job in CI; the Karatsuba overflow was only caught because a
  later change happened to make it crash.

## v13.0 — 2026-06-12

The 2026 performance campaign (PRs #82–#119), the full-codebase audit
(#120), and its follow-ups (#121–#125). BigMath now **beats GMP at
production scale (≥ 500k decimal digits) on multiplication, division,
parse, and ToString — in portable C++**. See README.md for the current
ratio tables and BENCHMARK.md for the canonical run.

### Performance

- **Multiplication**: NEON Shoup butterflies on aarch64 (radix-8/4/2 layers,
  ~4× per butterfly pass); 3-prime CRT NTT is now the effective default from
  1280 total limbs (Karatsuba below; Toom-3 dispatch window retired); fused
  MFA (Bailey 6-step) from 2^20 coefficients with row-chunked parallel
  stages and on-the-fly operand packing. Beats GMP from ~500k digits
  balanced and at all measured skews.
- **Division**: Burnikel–Ziegler odd-size bottom-padding fix (the "2^k+1
  family" pathology — up to 10× on ~half of real divisor sizes), BZ
  recursion basecase 512 → 128, re-swept Newton band frontier
  (8/1 @ 896 … 4/3 @ 131072), quotient-sized division for short-quotient
  shapes, cyclic (wrap-around) Newton products through the 2^26 CRT
  ceiling. Division beats GMP at all ≥ 20M-digit dividends.
- **Squaring**: routed to the CRT NTT self-multiply with duplicate-forward
  skip — 2.5–5.8× over the former Goldilocks NTTSquare (#122); threshold
  retuned 2048 → 640 limbs.
- **Decimal I/O**: divide-and-conquer ToString and parse fan their subtrees
  across the thread pool from 100k digits; both directions beat GMP from
  500k digits. Cold-parse Pow10 chain another −14…25% via the new square
  path.

### Correctness (audit #120, regression-tested in #121, fuzz-gated in #125)

- **FastDivision** returned silently wrong quotients when a remainder
  window's top two limbs equaled the divisor's top two (Möller-Granlund 3/2
  precondition violation; ~2^-128 on random operands, real on structured
  ones such as divisors near powers of two). Now special-cased exactly as
  GMP's `sbpi1_div_qr`.
- `a % b` took the divisor's sign into account; truncated-division remainder
  now follows the dividend only (`7 % -3` is `1`, was `-1`).
- CRT NTT length-ceiling guards ran after a narrowing cast and could be
  bypassed above 2^30 coefficients; KnuthDivision rewritten on 128-bit
  arithmetic (divided by zero under 64-bit limbs); carry/borrow fixes in
  vector `AddTo`, windowed `Add`, scalar `SubtractFrom`; `operator>>` sets
  `failbit` on non-numeric input; assorted dead/incorrect code removed.
- New CI gates: 285-test unit suite including audit regression pins and a
  seeded adversarial-limb fuzz (0/1/2^64−1-biased operands cross-checked
  against O(n²) references) — the input class behind every silent
  wrong-answer bug this cycle.

### Breaking changes

- `BigInteger::operator-()` is now value-returning and `const`; it
  previously negated **in place** and returned an lvalue reference.
  Replace `-x;` statements with `x = -x;`.
- `BigInteger(SizeT size, bool negative)` is now `explicit`. Implicit int
  conversions previously produced an N-limb **zero**, making expressions
  like `a + 5` and `a == 100` compile but silently misbehave; they no
  longer compile — use `BigIntegerBuilder::From(...)`.
- Public headers no longer inject `using namespace std;` into consumer
  translation units. Code that relied on the leak must qualify std names
  or add its own using-directive.
- Move construction/assignment are now real (previously deep-copied).

### Build

- CMake `include/` + `src/` split is the only build (`bigmath::bigmath`);
  `find_package(bigmath)` works against an installed prefix.
- Dispatch thresholds consolidated in
  `include/biginteger/build/DispatchThresholds.h`; `dispatch_tuner
  --emit-header` carries the complete division band set.

## v12.0 (untagged)

CMake build refactor (header-only → `include/` + `src/` split, static
library, CI matrix incl. Windows/MSYS2), 64-bit limb default, 3-prime CRT
NTT, multithreaded transforms, radix-4/8 fused butterflies, MFA, Newton
division with cached reciprocals, BigDecimal layer. Tagged retroactively
as part of the v13.0 history; see the merged PR record (#18–#81).

## v10.0 and earlier

Pre-campaign history: classic/Karatsuba/Toom/NTT multiplication,
Knuth/Burnikel–Ziegler division, calculator REPL. See git tags v6.7–v10.0.
