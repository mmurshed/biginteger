# Changelog

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
  against the faster multiply.
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
