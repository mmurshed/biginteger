# Changelog

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
