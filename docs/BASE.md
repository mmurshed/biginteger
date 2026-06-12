# Number representation: base and endianness

A technical reference for the foundational choice underlying every arithmetic algorithm in this BigInteger library: **little-endian limbs in base 2⁶⁴ stored as `uint64_t`** (the default since the 2026-05 limb refactor, PRs #18–#30). The historical layout — 32-bit limb values in `uint64_t` slots, base 2³² — remains available as a fallback via `-DBIGMATH_LIMB_64=0` and is documented here for context. This document explains the current choice, the history of how it got here, shows how it works in practice, lays out the consequences across the codebase, and explores the alternatives that were considered and rejected.

---

## Table of contents

1. [Scope](#scope)
2. [The choice in one paragraph](#the-choice-in-one-paragraph)
3. [Type definitions](#type-definitions)
4. [Defending the choice](#defending-the-choice)
   - [Why base 2^64 — and why it was base 2^32 for years](#why-base-264--and-why-it-was-base-232-for-years)
   - [Why little-endian](#why-little-endian)
   - [The legacy 32-in-64 layout (historical)](#the-legacy-32-in-64-layout-historical)
   - [Why a runtime `BaseT` parameter](#why-a-runtime-baset-parameter)
5. [How it works in practice](#how-it-works-in-practice)
   - [Single-limb scalar multiply](#single-limb-scalar-multiply)
   - [Carry propagation in addition](#carry-propagation-in-addition)
   - [Borrow propagation in subtraction](#borrow-propagation-in-subtraction)
   - [Schoolbook multiplication inner loop](#schoolbook-multiplication-inner-loop)
   - [NTT input split](#ntt-input-split)
   - [Decimal I/O via grouped-digit chunks](#decimal-io-via-grouped-digit-chunks)
6. [Benchmark consequences](#benchmark-consequences)
7. [Optimizations already implemented](#optimizations-already-implemented)
8. [Future opportunities](#future-opportunities)
9. [Explored but rejected](#explored-but-rejected)
10. [References](#references)

---

## Scope

This document explains why a single design decision — "limbs are true 64-bit values in `uint64_t` slots, stored low-limb-first, with a runtime `BaseT` parameter" — propagates through the entire arithmetic stack, what the alternative layouts cost, and how the library migrated from its original 32-bit-value layout to the current one.

Companion documents cover the algorithms built on this foundation:

- [MULTIPLICATION.md](MULTIPLICATION.md)
- [DIVISION.md](DIVISION.md)

This document is the prerequisite for both.

Assumed reader: a working C++ engineer familiar with `uint32_t` / `uint64_t` arithmetic, bit operations, and the basic shape of multi-precision algorithms.

---

## The choice in one paragraph

A `BigInteger` (`include/biginteger/BigInteger.h`) wraps a `std::vector<DataT>` where `DataT = uint64_t`, plus a sign bit. In the default build (`BIGMATH_LIMB_64=1`), each element of the vector holds a **full 64-bit limb value** (`0 ≤ limb < 2⁶⁴`); carry and product headroom comes from `unsigned __int128` (`ULong128`) accumulators in the inner loops, not from unused bits in the storage slot. Limbs are stored **little-endian**: `vec[0]` is the least-significant limb. The base is exposed as a runtime parameter `BaseT base`. Because the numeric value 2⁶⁴ does not fit in `BaseT` (`int64_t`), the production base is the **sentinel `Base2_64 = 0`** — 0 is never a legal numeric base, so code that only *compares* the base uses the sentinel directly, while code that needs to *use* the base as a shift or mask takes a separate `Base2_64` branch and operates on the limb width via the `LimbBits` / `LimbBase` / `LimbMask` helpers. Non-power-of-two bases (`Base10`, `Base100`, `Base100M`) share the same algorithm code paths for I/O and testing. Compiling with `-DBIGMATH_LIMB_64=0` restores the historical layout: 32-bit limb values in `uint64_t` slots, base `Base2_32`, with the upper 32 bits of each slot serving as carry headroom.

```
                    most significant ──→
        ┌──────────┬──────────┬──────────┬──────────┐
   a =  │  a[3]    │  a[2]    │  a[1]    │  a[0]    │       sign = false  ⇒  a ≥ 0
        └──────────┴──────────┴──────────┴──────────┘
            ↑                              ↑
         high limb                      low limb (vec[0])

   default (BIGMATH_LIMB_64=1): each limb is a uint64_t holding a value in [0, 2⁶⁴).
   carry headroom comes from __uint128_t accumulators in the inner loops.

   legacy (-DBIGMATH_LIMB_64=0): each uint64_t slot holds a value in [0, 2³²);
   the high 32 bits of each slot are headroom for carry accumulation.
```

---

## Type definitions

From `include/biginteger/common/Constants.h`:

| name | underlying | width | role |
|---|---|---|---|
| `DataT` | `uint64_t` | 64-bit | one limb — a true 64-bit value by default; a 32-bit value in a 64-bit slot under `-DBIGMATH_LIMB_64=0` |
| `BaseT` | `int64_t` | — | the radix; `CurrentBase = Base2_64 = 0` (sentinel) in the default build |
| `ULong` | `uint64_t` | 64-bit | accumulator where one 64-bit value suffices |
| `ULong128` | `unsigned __int128` | 128-bit | accumulator for 64×64 products and carry chains — the workhorse of the default limb width |
| `Long` | `int64_t` | 64-bit signed | signed accumulator (legacy 32-bit borrow paths; generic-base paths) |
| `Int` | `int32_t` | 32-bit signed | loop counters, signed limb diffs |
| `SizeT` | `uint32_t` | 32-bit | vector sizes, indices |

`Constants.h` also defines compile-time limb-width helpers that drive the carry/borrow/shift idioms throughout the library:

```cpp
#if BIGMATH_LIMB_64
    constexpr SizeT LimbBits = 64;
    constexpr ULong128 LimbBase = (ULong128)1 << 64;   // 2^64
    constexpr ULong LimbMask = 0xFFFFFFFFFFFFFFFFULL;
    constexpr BaseT CurrentBase = Base2_64;            // sentinel 0
#else
    constexpr SizeT LimbBits = 32;
    constexpr ULong LimbBase = 1ULL << 32;             // 2^32
    constexpr ULong LimbMask = 0xFFFFFFFFULL;
    constexpr BaseT CurrentBase = Base2_32;
#endif
```

`BigInteger::Base()` returns `CurrentBase` — `Base2_64` (the 0 sentinel) by default, `Base2_32` under `-DBIGMATH_LIMB_64=0`. Defined bases (most are I/O or test helpers):

```
Base2       = 2          ← used by tests, bit-level operations
Base10      = 10         ← parser inner-most digit chunk
Base100     = 100        ← legacy I/O helper
Base100M    = 100 000 000 (10⁸) ← legacy I/O helper
Base2_8     = 256        ← used by some tests
Base2_10    = 1024       ← used by some tests
Base2_16    = 65 536     ← Goldilocks-NTT internal coefficient base
Base2_31    = 2 147 483 648 (2³¹) ← used in normalization
Base2_32    = 4 294 967 296 (2³²) ← legacy production limb base (-DBIGMATH_LIMB_64=0)
Base2_64    = 0 (sentinel)        ← PRODUCTION limb base (default build)
```

The `Base2_64` sentinel deserves a sentence of its own: 2⁶⁴ does not fit in `int64_t`, so the constant is 0 — a value that can never be a real base. Dispatch code that merely compares (`if (base == Base2_64)`) works unchanged; arithmetic code that would have computed `% base` or `/ base` instead takes a dedicated branch using shift/mask on the full limb (or relies on natural 64-bit wraparound, e.g. `FastDivision`'s low-digit helper, where the `(ULong)base = 0` sentinel makes the historical mask expression degenerate and a branch is required).

The production code path is `Base2_64` everywhere outside the NTT internal transform (which by default uses three ~30-bit CRT primes with 32-bit coefficients — see [NTT input split](#ntt-input-split)) and decimal I/O (which parses in `Base10_18 = 10¹⁸` chunks and formats in `Base10_19 = 10¹⁹` chunks — both are chunking conventions, not storage bases).

---

## Defending the choice

### Why base 2^64 — and why it was base 2^32 for years

**The current answer is simple: half the limbs, and `__uint128_t` accumulators provide all the headroom the inner loops need.** Every loop over limbs — carry propagation, schoolbook partial products, normalization, reciprocal steps — runs half as many iterations as the 32-bit layout for the same number. The 2026-05 limb refactor (PRs #18–#30) made this the default after head-to-head benchmarks showed `BIGMATH_LIMB_64=1` matching or beating the 32-bit layout on every measured operation (see the bench note in `Constants.h`: small/mid multiplication −25 to −60 %, skewed division −36 to −59 %, parse −24 to −34 %, ToString −19 to −28 %; NTT-bound 1M-digit multiplication flat, because that path is coefficient-bound, not limb-bound). The measured wins are recorded in detail in [MULTIPLICATION.md](MULTIPLICATION.md) and [DIVISION.md](DIVISION.md).

**The historical rationale for base 2³² — kept here for context, and because it still governs the `-DBIGMATH_LIMB_64=0` fallback — was a single property:** every native multiplication of two limbs fits in a single 64-bit register.

In base 2³², `limb × limb < 2³² · 2³² = 2⁶⁴`. The product is exactly a `uint64_t`. Adding two such products plus an in-place limb plus a carry all fit in 64 bits with one bit of slack:

```
   prod_max = (2³² − 1)² = 2⁶⁴ − 2³³ + 1
   adding (r[i+j] + carry) where each ≤ 2³² − 1:
   total_max = (2⁶⁴ − 2³³ + 1) + 2(2³² − 1) = 2⁶⁴ − 1.
   ── fits in uint64_t exactly, no overflow, no extra carry plumbing.
```

That made the schoolbook inner loop three ARM64 instructions per partial product with single-register accumulators, and the worry at the time was that base 2⁶⁴ in pure C++ — `MUL` + `UMULH` per multiply, two-register accumulators, `ADCS` carry chains — would lose what it gained from halving the limb count, because compilers don't schedule `__uint128_t` arithmetic as well as GMP's hand-written assembly does.

**Why the 64-bit refactor won anyway.** The original analysis predicted base 2⁶⁴ in pure C++ would be net "about the same throughput, worse register pressure". Measurement said otherwise:

- Modern AppleClang/GCC turn `(ULong128)a * b + r + carry` into a tight `MUL`/`UMULH`/`ADDS`/`ADCS` sequence; the feared register-pressure penalty did not materialize on ARM64.
- Halving the limb count helps *every* loop in the library, not just the multiply leaf — carry propagation, comparison, copying, normalization, division fix-up steps, parse/format inner loops. Those wins compound.
- The 32-bit layout's slot headroom was only ever a substitute for a wider accumulator. `ULong128` provides the same headroom explicitly, without halving the information density of every load and store.
- The one place the limb width genuinely doesn't matter — the NTT-bound large-multiplication path — measured flat, exactly as predicted (the transform cost depends on coefficient count, not source limb width).

The 32-bit path is retained behind `-DBIGMATH_LIMB_64=0` for A/B testing and as a fallback for platforms where `__int128` codegen is poor.

### Why little-endian

The choice is **structural for arithmetic** (not a debate). In multi-precision arithmetic, carries propagate from low to high. Indexing a vector `a` with `a[0]` as the least significant limb means a single forward loop `for (i = 0; i < n; ++i)` handles carry propagation in addition, subtraction, and scalar multiplication.

Reverse this — store big-endian — and the same loop has to either decrement an index (`for (i = n-1; i >= 0; --i)`) or pre-reverse the buffer, neither of which the compiler vectorizes as cleanly. More importantly: extending a number (adding limbs at the high end during a carry that propagates past the current top) becomes `push_back(carry)` in little-endian and `insert(begin(), carry)` — an O(n) shift — in big-endian.

Little-endian also matches the on-disk and on-network conventions of most existing big-number libraries (GMP `mpn`, OpenSSL `BIGNUM` internals, Python's `_PyLong` digits, Java `BigInteger.mag[]` *internally is big-endian but exposed as such only because it pre-dates these tradeoffs*). Interop is easier when both sides agree.

A subtle benefit: **trimming leading zeros** (the canonical "is this number normalized" operation) is `while (vec.size() > 1 && vec.back() == 0) vec.pop_back();` — O(amortized 1) at the high end of the vector. The same operation in big-endian would erase elements from the front, an O(n) per-erase cost.

### The legacy 32-in-64 layout (historical)

> **Historical section.** Everything below describes the original design — 32-bit limb values in `uint64_t` slots — which was the production layout until the 2026-05 limb refactor and now lives behind `-DBIGMATH_LIMB_64=0`. It is kept because the rationale explains several code shapes that survive in the fallback paths.

The "2× memory for 32 bits of value" tradeoff was deliberate. Three reasons:

1. **Carry headroom.** A multiply-accumulate result `a[i] * b[j] + r[i+j] + carry` always fits in 64 bits with the 32-bit limb representation. Storing limbs in `uint32_t` would have required a separate `uint64_t` temporary for the accumulator at every loop iteration; keeping limbs in `uint64_t` slots let the accumulator live in the same register as the destination limb. In the current 64-bit-limb default, this role is played explicitly by `ULong128` accumulators instead — same headroom, expressed in the accumulator type rather than wasted storage bits.

2. **Alignment for the hybrid leaf.** The 64-bit hybrid Karatsuba leaf (see [Schoolbook multiplication inner loop](#schoolbook-multiplication-inner-loop)) packed pairs of 32-bit limbs into `uint64_t` for the schoolbook; the `vec[2k] | vec[2k+1] << 32` pack reads two 64-bit slots and produces one packed value with a single AArch64 `ORR`. This pack/unpack dance was the bridge that proved 64-bit inner arithmetic was a win — and the full refactor then removed the need for it by making 64-bit limbs the storage format.

3. **NTT coefficient extraction.** The Goldilocks NTT splits each 32-bit limb into two 16-bit coefficients via `(limb & 0xFFFF)` and `(limb >> 16)`; both extractions on a `uint64_t` slot are one instruction.

The memory argument cut the other way too: a 32-in-64 layout stores 2× the bytes per bit of number. Under the current default, the same `uint64_t` slot carries twice the value — so the 64-bit layout halves memory traffic per digit of number, which is part of why the small/mid-size wins materialized.

### Why a runtime `BaseT` parameter

Most modern big-number libraries hard-code their base at compile time. This library exposes `BaseT base` as a runtime parameter to every algorithm. Why?

**Historical**: the algorithm code was originally written to be base-agnostic for pedagogical clarity (the classical schoolbook algorithms in *TAOCP* are presented in terms of an abstract base B). The runtime parameter let the same code paths handle `Base10` for tutorial-level traces and the production base for real work.

**Practical**: a small set of code paths still need non-power-of-two bases. Specifically:

- `ClassicDivision::DivModTo(vec, divisor, base)` is called from `ToStringLinearAppend` with `base = CurrentBase` and `divisor = Base10_19`. The "base" parameter here is the **storage** base of the dividend vector, independent of the divisor.
- Test scaffolding builds small numbers in `Base10` or `Base16` for human readability before converting.
- Some legacy I/O helpers operate in `Base100M`.

The runtime parameter costs essentially nothing because of how the dispatching is structured: every hot inner loop branches on `base == Base2_64` (and, in code that retains the fallback, `base == Base2_32`) *once at the top of the function* and falls into a fully-specialized fast path. The branch is highly predictable (always the same in production), so the compiler keeps the fast path tight. Note again the sentinel discipline: the `Base2_64` comparison is safe everywhere, but a generic `% base` / `/ base` expression must never be reached with the sentinel — the specialized branch uses `LimbMask` / shift-by-64 / 64-bit wraparound instead.

In the rare loops where the base check was historically *inside* the inner loop (the Karatsuba helpers `AddPtr` / `AddToPtr` / `SubtractFromPtr`), the 2026-05 rewrite hoisted it out — see [MULTIPLICATION.md §Optimizations](MULTIPLICATION.md#optimizations-already-implemented). After the hoist, the runtime parameter has no measurable cost in any production path.

If the parameter ever did become a bottleneck, the path forward is a templated `BaseT` (compile-time constant) — feasible because `CurrentBase` is the only base any production caller passes. The templating would yield a non-production path with no specialization (the generic `% base` / `/ base` arithmetic), which is fine since those paths are only exercised by tests.

---

## How it works in practice

The base-and-endian choice cashes out as concrete instruction sequences. This section walks through the most performance-critical patterns as they exist in the default (`BIGMATH_LIMB_64=1`) build; the legacy 32-bit idioms are noted where they differ.

### Single-limb scalar multiply

`a` is a vector of 64-bit limbs. Multiply by a 64-bit scalar `b` in place (the `Base2_64` fast path in `ClassicMultiplication::MultiplyTo`):

```cpp
ULong carry = 0;
for (SizeT j = aStart; j <= aEnd; ++j) {
    ULong128 p = (ULong128)a[j] * b + carry;
    a[j]  = (DataT)p;          // low 64 bits become the new limb
    carry = (ULong)(p >> 64);  // high 64 bits become the next carry
}
if (carry) a.push_back((DataT)carry);
```

On ARM64 this compiles to `MUL` (low half) + `UMULH` (high half) + an add-with-carry per iteration. The carry is at most 2⁶⁴ − 1 (the high half of the previous product), so a single push suffices after the loop; the vector's exponential growth strategy makes the amortized cost of `push_back` O(1).

The legacy 32-bit path does the same thing with a `uint64_t` product — one `UMULL` per iteration, mask for the low 32 bits, shift for the carry — but runs **twice as many iterations** for the same number. That trade (cheaper iteration × 2× iteration count vs. richer iteration × 1× count) is the whole limb-width question in miniature, and measurement settled it in favor of 64-bit.

### Carry propagation in addition

The `Base2_64` branch in `Add` (`src/algorithms/Addition.cpp`):

```cpp
ULong128 carry = 0;
for (SizeT i = 0; i < n; ++i) {
    ULong128 sum = (ULong128)a[i] + b[i] + carry;
    r[i]  = (DataT)(sum & LimbMask);   // low 64 bits
    carry = sum >> 64;                 // 0 or 1
}
// propagate any final carry; a bare += can itself overflow the slot
```

`a[i] + b[i] + carry` can exceed 2⁶⁴, so the accumulator is `ULong128`; the high bits are exactly the carry-out, by construction. (An alternative formulation — plain `uint64_t` adds with `__builtin_add_overflow` — compiles to the same `ADCS` chain; the `ULong128` form was chosen for consistency with the multiply paths.)

The legacy 32-bit branch needs only a `uint64_t` accumulator (`a[i] + b[i] + carry ≤ 2³³ − 1`), which was one of the original arguments for the 32-in-64 layout. The cost of the wider accumulator in the 64-bit path is negligible next to halving the iteration count.

### Borrow propagation in subtraction

A signed `Long` can't hold a full 64-bit limb, so the `Base2_64` branch (`src/algorithms/Subtraction.cpp`) does unsigned arithmetic with explicit borrow detection:

```cpp
ULong borrow = 0;
for (SizeT i = 0; i < n; ++i) {
    ULong ai = a[i], bi = b[i];
    // Compute ai - bi - borrow with two-step borrow detection.
    ULong t1 = ai - borrow;
    ULong borrow1 = (ai < borrow) ? 1 : 0;
    ULong diff = t1 - bi;
    ULong borrow2 = (t1 < bi) ? 1 : 0;
    borrow = borrow1 + borrow2;
    r[i] = (DataT)diff;
}
```

Unsigned wraparound makes each step exact: when the subtraction underflows, the wrapped value is precisely `ai − bi + 2⁶⁴`, which is the correct limb, and the comparison records the borrow. The compiler lowers this to a `SUBS`/`SBCS`-style flag chain.

The legacy 32-bit branch is simpler — `a[i] − b[i] − borrow` fits in `int64_t` (`Long`), so a signed accumulator with an `if (diff < 0)` fixup suffices.

### Schoolbook multiplication inner loop

In the default build the schoolbook simply **is** 64-bit: limbs are already 64-bit values, partial products use `ULong128`, and dispatch (`include/biginteger/algorithms/Multiplication.h`) keeps operands on classic schoolbook up to 96 total limbs (a much larger classic window than the 32-bit layout had, because each 64-bit iteration does the work of four 32-bit ones):

```cpp
for (SizeT i = 0; i < nb; ++i) {
    ULong bi = b[i];
    if (bi == 0) continue;
    ULong carry = 0;
    for (SizeT j = 0; j < na; ++j) {
        ULong128 prod = (ULong128)a[j] * bi + r[i+j] + carry;
        r[i+j] = (ULong)prod;
        carry  = (ULong)(prod >> 64);
    }
    r[i + na] = carry;
}
```

**Historical note — the base-promotion trick.** Before the limb refactor, this exact loop existed as a *hybrid leaf* inside `KaratsubaMultiplication::MultiplyClassicPtr`: the `Base2_32` path packed pairs of 32-bit limbs into 64-bit values (`a64[k] = a[2k] | (a[2k+1] << 32)`), ran the schoolbook above on the packed representation, and unpacked the result back to 32-bit limbs. The inner multiplication ran in base 2⁶⁴ (recovering the per-multiply efficiency GMP gets natively) while storage stayed in base 2³². It was worth ~2× at the leaf — and it was the empirical proof that motivated making 64-bit limbs the storage format outright, which eliminated the pack/unpack overhead. The hybrid leaf survives in `KaratsubaMultiplication.h` as the `Base2_32` fallback path; its stack buffers cover 64 packed limbs (= 128 source limbs), comfortably above the Karatsuba leaf threshold of **32** (`BIGMATH_KARATSUBA_THRESHOLD`).

### NTT input split

**Default path: 3-prime CRT NTT with 32-bit coefficients.** `NTTMultiplicationCrt.h` runs the convolution modulo three ~30-bit NTT-friendly primes (2 013 265 921, 469 762 049, 1 811 939 329) and reconstructs via CRT. Each 64-bit limb splits into exactly **two 32-bit coefficients**:

```
   one 64-bit limb         →        two 32-bit CRT coefficients
                                    ┌──────────────────┬──────────────────┐
   limb (uint64, value      ──→     │   limb >> 32     │ limb & 0xFFFFFFFF│
   in [0, 2⁶⁴))                     └──────────────────┴──────────────────┘
                                          (high half)        (low half)
```

The 32-bit coefficient split halves the transform length relative to a 16-bit split, the per-prime butterflies run on 32-bit modular lanes (single `MUL` plus a simple reduction; Shoup-multiplied twiddles, NEON-vectorized on aarch64), and the coefficients pack two-per-cache-word. The CRT length ceiling is 2²⁶ ≈ 67 M coefficients ≈ 640 M-digit operand pairs. Dispatch enters this path for effectively every NTT-sized multiply (`BIGMATH_NTT_CRT_THRESHOLD = 256` total limbs, below the NTT entry point of 1280).

**Fallback path: single-prime Goldilocks with 16-bit coefficients** (`-DBIGMATH_NTT_CRT=0`, in `NTTMultiplication.h`). The Goldilocks prime `P = 2⁶⁴ − 2³² + 1` caps input coefficients at 16 bits for its convolution-accumulation bound. Under the default 64-bit limbs, each limb splits into **four 16-bit coefficients**:

```cpp
for (SizeT i = 0; i < a.size(); ++i) {
    SizeT j = i * 4;
    fa[j]     = a[i] & 0xFFFFULL;
    fa[j + 1] = (a[i] >> 16) & 0xFFFFULL;
    fa[j + 2] = (a[i] >> 32) & 0xFFFFULL;
    fa[j + 3] = (a[i] >> 48) & 0xFFFFULL;
}
```

(Under `-DBIGMATH_LIMB_64=0` the same path splits each 32-bit limb into two 16-bit coefficients.) The convolution overflow bound is ~2³¹ limbs (≈ 16 GB operands) — comfortably beyond practical sizes — which is why no Schönhage-Strassen layer exists in this library.

Either way, after the inverse transform(s), carry propagation runs over the convolution output at the coefficient width, and coefficient groups are reassembled into 64-bit limbs (`NttFinalizeBase2_64`; two 32-bit CRT coefficients or four 16-bit Goldilocks coefficients per output limb). The relationship is exact and lossless in both layouts — the transform cost is governed by the coefficient count, which is the same per byte of input regardless of limb width. This is precisely why the limb refactor measured flat on NTT-bound sizes.

### Decimal I/O via grouped-digit chunks

Decimal parsing reads input digits 18 at a time (`Base10_18 = 10¹⁸`), multiplying the accumulated number by 10¹⁸ and adding each chunk:

```cpp
// in ParseUnsignedLinear (src/common/Parser.cpp):
while (pos <= end) {
    ULong chunk = ParseChunk(num, pos, Base10_18_Zeroes);       // 18 ASCII digits → uint64
    pos += Base10_18_Zeroes;
    ClassicMultiplication::MultiplyTo(r, Base10_18, CurrentBase); // r ← r * 10¹⁸
    AddTo(r, chunk, CurrentBase);                                  // r ← r + chunk
}
```

`MultiplyTo(r, Base10_18, CurrentBase)` takes the `Base2_64` scalar-multiply fast path shown earlier: each 64-bit limb times the 60-bit constant accumulates in `ULong128`.

Decimal output goes the other way, in **19-digit chunks** (`Base10_19 = 10¹⁹`, the largest power of 10 that fits in a 64-bit limb):

```cpp
// in ToStringLinearAppend (src/common/Parser.cpp):
chunks.push_back((ULong)ClassicDivision::DivModTo(r, Base10_19, CurrentBase));
```

`DivModTo` divides `r` by 10¹⁹ in place and returns the remainder — the next 19 decimal digits to format. The per-limb step uses a Möller-Granlund "div2by1" reciprocal divide (`UMULH` + 128-bit add + ≤2 fixups), with the reciprocal precomputed once per call from the invariant 10¹⁹ divisor. Above the linear loop sits the divide-and-conquer `Pow10` chain (`BuildDecimalDcChain`) — see [STRING_CONVERSION.md](STRING_CONVERSION.md) for the full I/O pipeline.

Both directions exploit the fact that the chunk constants fit in a 64-bit register while the limb storage representation stays untouched. The base-2⁶⁴ storage and the base-10ⁿ I/O chunking are independent design choices that compose cleanly.

---

## Benchmark consequences

**Current numbers live in [BENCHMARK.md](../BENCHMARK.md)** (canonical run 2026-06-12 at v13.0), which reflects the full current configuration: 64-bit limbs, 3-prime CRT NTT with NEON Shoup butterflies, multithreading, and the 2026-06 dispatch retunes. Per-operation analysis is in [MULTIPLICATION.md](MULTIPLICATION.md) and [DIVISION.md](DIVISION.md). The tables below are **not** current.

> **Historical snapshot (pre-2026-05): 32-in-64 limbs, single-prime Goldilocks NTT, single-threaded.** Kept because the *structure* of the analysis — which costs are limb-bound and which are transform-bound — remains instructive.
>
> Schoolbook leaf (post-hybrid, pre-refactor):
>
> | op | BigMath ms | GMP ms | BM / GMP |
> |---|---:|---:|---:|
> | mul 1 000 × 1 000 digits | 0.003 | 0.001 | 3.5 × |
> | mul 5 000 × 5 000 digits | 0.043 | 0.012 | 3.5 × |
>
> NTT-bound (input-split is the binding constraint):
>
> | op | BigMath ms | GMP ms | BM / GMP |
> |---|---:|---:|---:|
> | mul 1 000 000 × 1 000 000 digits | 37.5 | 8.9 | 4.2 × |
>
> Skewed division (Newton iteration over NTT internal mults):
>
> | op | BigMath ms | GMP ms | BM / GMP |
> |---|---:|---:|---:|
> | div 200 000 / 50 000 digits | 22 | 1.7 | 14 × |
> | div 500 000 / 100 000 digits | 55 | 4.6 | 12 × |

Two things about how the snapshot aged are worth recording:

- The snapshot's claim that "switching to base 2⁶⁴ would leave the NTT-bound ratio unchanged" was **validated by the switch actually happening**: the limb refactor measured flat on the 1M-digit multiply, exactly as the analysis predicted (the transform is coefficient-bound). The NTT-bound gap was subsequently closed by different levers entirely — the CRT NTT, NEON Shoup butterflies, MFA, and threading.
- The snapshot's limb-bound rows (schoolbook leaf, skewed division) are precisely where the 64-bit refactor and the later division work landed their wins. See BENCHMARK.md for where each ratio stands today.

---

## Optimizations already implemented

A summary of representation-related optimizations that landed and stuck. (Each is covered in detail in the algorithm-specific documents.)

### Full 64-bit limb representation (the headline win — landed 2026-05, PRs #18–#30)

`DataT` holds true 64-bit values; every carry/borrow/product idiom got a `Base2_64` branch built on `ULong128` accumulators (or unsigned-wraparound borrow tracking where 128-bit isn't needed). Halved the limb count for every loop in the library. Measured wins (M1 Max, `bench_vs_gmp`): small/mid mul −25 to −60 %, skewed div −36 to −59 %, parse −24 to −34 %, ToString −19 to −28 %; NTT-bound sizes flat as predicted. Details in [MULTIPLICATION.md](MULTIPLICATION.md) and [DIVISION.md](DIVISION.md).

### 64-bit hybrid leaf (historical stepping stone)

`KaratsubaMultiplication::MultiplyClassicPtr`'s `Base2_32` path packs pairs of 32-bit limbs into 64-bit values, runs schoolbook in 64-bit space (`ULong128` accumulator), and unpacks — ~2× over scalar 32-bit schoolbook at the leaf (mul 5 000 × 5 000 went from 9.3× to 3.5× vs GMP at the time). This was the experiment that proved the full refactor's premise; under the default build the schoolbook is natively 64-bit and the hybrid survives only as the `-DBIGMATH_LIMB_64=0` fallback path.

### Hoisted base branch in inner loops

`AddPtr`, `AddToPtr`, `SubtractFromPtr` in `KaratsubaMultiplication.h` had the runtime base check pulled out of the inner loop and replaced with parallel specialized paths. The branch is highly predictable, but the loop body without it is one fewer instruction per iteration. Modest 0–4 % net wins. Discussed in [MULTIPLICATION.md §Optimizations](MULTIPLICATION.md#optimizations-already-implemented).

### NTT direct-fill of coefficient buffers

Both NTT paths build their transform input vectors directly from the limbs — two 32-bit coefficients per limb on the CRT path, four 16-bit on the Goldilocks fallback — skipping intermediate split-vector temporaries. The power-of-two limb format makes this a few shift/mask lines per limb. Saves allocation + copy overhead proportional to input size.

### `FastDivision` limb-specialized paths

Normalize, scalar divisor, and remainder unnormalize have `Base2_64` (and legacy `Base2_32`) code paths using `ULong128` accumulators and `UMULH`-style intrinsics, with the sentinel-aware low-digit/carry helpers branching once at the top. Covered in [DIVISION.md §FastDivision](DIVISION.md#fast-division-knuth-algorithm-d).

### Decimal D&C chain using `NewtonDivision::Divider`

`Parser.h` / `src/common/Parser.cpp` (`BuildDecimalDcChain`) builds a chain of `Divider` instances, each holding a `Pow10(d)` value (in `CurrentBase` limbs) and its precomputed reciprocal — the foundation of the divide-and-conquer ToString speedup. The underlying `ApproxReciprocal` / `DivideChunk` operations run on `ULong128` accumulators over 64-bit limbs.

---

## Future opportunities

Ranked by expected ROI per unit of effort.

### Templated `BaseT` (compile-time constant)

Replace the runtime `BaseT base` parameter with a template parameter `template <BaseT base>`. For production callers using `CurrentBase`, the compiler would inline the constant and eliminate any residual conditional branching that the hoisted-branch approach doesn't fully cover.

Effort: moderate (touches every algorithm signature). Risk: low — the only callers passing non-production bases are tests, which can be templated trivially or pinned to a runtime variant for their specific bases. Estimated win: 1–3 % across all base-sensitive paths. Net: marginal, but cumulative across the whole library.

### Landed: full 64-bit limb refactor (2026-05, PRs #18–#30)

This section used to be the largest "future opportunity" in the document; it shipped, and it is now the default build. The prediction table it carried (marginal gains at the Karatsuba leaf, zero at NTT-bound sizes, meaningful wins on Newton division and linear parse/ToString loops) held up qualitatively — NTT-bound multiplication was indeed flat, and the limb-bound paths improved by the margins recorded in `Constants.h` and the algorithm docs. See [Optimizations already implemented](#optimizations-already-implemented) above, and [MULTIPLICATION.md](MULTIPLICATION.md) / [DIVISION.md](DIVISION.md) for the measured deltas. The legacy layout remains behind `-DBIGMATH_LIMB_64=0` for A/B testing.

### Aligned-storage hint to the allocator

`std::vector<DataT>` doesn't guarantee 16-byte alignment of its data. Aligned loads are slightly faster on M1 and matter for NEON paths. A custom allocator or `alignas`-declared inline storage for small numbers could give 2–5 % on hot inner loops. Marginal, but cheap.

> **Historical note.** A "NEON-vectorized 32-bit limb arithmetic" item also used to live here (lane-parallel `UMLAL` partial products with deferred carry propagation, estimated 1.5–2× at the schoolbook leaf). It targeted the legacy 32-bit limb layout and is mooted by the 64-bit default; the NEON effort went into the CRT NTT butterflies instead (Shoup-multiplied twiddles — see MULTIPLICATION.md), which is where vector lanes actually pay off in the current architecture.

---

## Explored but rejected

### Base 10 / Base 10ⁿ as the storage representation

**Rejected.** Tempting because decimal I/O becomes trivial (no parse/format conversion overhead). Killed by every other operation: arithmetic in base 10 requires `% 10` and `/ 10` per limb step — both extremely slow compared to the bit-mask and shift that a power-of-two base uses. The classical schoolbook in pure base 10 would be 10–20× slower. Real implementations using decimal storage (e.g. `_decimal128` formats) are designed for fixed-precision applications where I/O cost dominates; arbitrary-precision libraries universally use power-of-two bases.

The library does internally use base 10¹⁸ / 10¹⁹ as **chunking conventions for I/O**, but not as a storage representation. The chunks are immediately multiplied into binary limbs during parsing.

### Base 2⁶³ (one bit reserved for sign)

**Rejected.** Some libraries reserve the top bit of each limb for the sign of the partial product, eliminating the separate sign field. Costs: every arithmetic operation must mask the sign bit, and the per-limb effective range shrinks. The savings (one boolean per `BigInteger`) are negligible compared to the loss of clean full-width arithmetic. The sign-as-separate-bool design is universally preferred in modern libraries.

### Signed-digit (non-adjacent form) representation

**Rejected.** Signed-digit representations (digits in `{−B/2, ..., +B/2}`) make addition carry-free at the cost of a non-canonical representation (multiple ways to represent the same number) and conversion overhead at boundaries (I/O, comparison). Useful in some cryptographic applications (Montgomery ladders, ECC scalar multiplication) but the library's use case (general-purpose big-number arithmetic) doesn't benefit. Conversion costs would dominate.

### Big-endian limb storage

**Rejected.** See [Why little-endian](#why-little-endian). The carry-propagation direction (low → high) matches the array iteration direction; trim-leading-zeros is O(1) at the high end of a `std::vector`. Big-endian inverts both, with no compensating benefit.

### Compile-time `BaseT` only (templated, no runtime parameter)

**Considered but not yet executed.** Discussed in [Future opportunities](#future-opportunities). The runtime parameter has no measurable cost in the current code (the `base == Base2_64` / `base == Base2_32` checks are hoisted out of inner loops), so the templating refactor is not currently prioritized. Could be done if testing infrastructure ever wants to template-instantiate algorithm code paths for fuzz testing of non-production bases.

### Variable-length limbs (e.g., 28-bit for digit-recurrence FFT)

**Not considered seriously.** Some FFT-multiplication schemes use 28-bit or 23-bit limbs to fit results in a `double` mantissa for IEEE 754 floating-point NTT. The integer-NTT paths used here don't benefit; the limb size is decoupled from the transform coefficient size by the input split (32-bit coefficients on the CRT path, 16-bit on the Goldilocks fallback).

### Native `__uint128_t` as the limb type

**Rejected.** `__uint128_t` is a GCC/Clang extension, not part of standard C++. Using it as `DataT` would lock out MSVC (where `__int128` doesn't exist) and complicate the `std::vector` allocator path. The `ULong128` accumulator is used pervasively but in narrow scopes (carry chains, NTT, division normalize) where the platform requirement is acceptable; pushing it to the limb type would force the whole library to be GCC/Clang-only at the storage level too, for a representation whose products would then need 256-bit accumulation.

The library is already effectively GCC/Clang-only because of `unsigned __int128` and `__builtin_*` usage — but that's a contained surface, and could in principle be replaced with portable C++23 `<stdckdint.h>` operations.

---

## References

### Multi-precision arithmetic foundations

- Knuth, D. E. *The Art of Computer Programming, Vol. 2: Seminumerical Algorithms*, §4.3 — the canonical treatment of positional notation, limb-based arithmetic, and base choice.
- [Brent, R. P. and Zimmermann, P. — *Modern Computer Arithmetic* (Cambridge, 2010)](https://members.loria.fr/PZimmermann/mca/pub226.html) — modern reference; Chapter 1 covers integer representation including the base/limb/endianness tradeoffs.

### Hardware multiplication primitives

- [ARM Architecture Reference Manual: `UMULL` / `UMULH` / `UMADDL`](https://developer.arm.com/documentation/ddi0596/2021-12/Base-Instructions/UMULL--Unsigned-Multiply-Long-) — `MUL`/`UMULH` underlie the 64-bit limb arithmetic on ARM64; `UMULL` underlies the legacy 32-bit paths.
- [Intel 64 and IA-32 Architectures Software Developer's Manual, Vol. 2: `MUL` / `MULX`](https://www.intel.com/content/www/us/en/developer/articles/technical/intel-sdm.html) — x86 multiplication primitives.

### Reference implementations and conventions

- [GMP — The GNU Multiple Precision Arithmetic Library](https://gmplib.org/) — uses base 2⁶⁴ (configurable to 2³² at build time on 32-bit platforms). Little-endian limbs. The same limb width this library now defaults to.
- [GMP manual: "Low-level Functions" (`mpn_*`)](https://gmplib.org/manual/Low_002dlevel-Functions) — documentation of the `mp_limb_t` and `mp_size_t` types and the calling conventions, including the little-endian convention.
- [OpenSSL `BIGNUM` internals (`bn.h`)](https://github.com/openssl/openssl/blob/master/include/openssl/bn.h) — `BN_ULONG` is either 32-bit or 64-bit depending on the platform; little-endian limbs.
- [Python `_PyLong` digit representation](https://github.com/python/cpython/blob/main/Include/cpython/longintrepr.h) — uses 30-bit "digits" stored in 32-bit words on 64-bit platforms, leaving 2 bits of carry headroom. Conceptually similar tradeoff to this library's legacy 32-in-64 layout.
- [Java `BigInteger.mag[]`](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/math/BigInteger.html) — `int[]` of 32-bit magnitudes, stored big-endian (`mag[0]` is most significant). The big-endian convention is a holdover from early Java's I/O orientation; modern libraries universally pick little-endian.

### Specific design references

- [Möller, N. and Granlund, T. — "Improved Division by Invariant Integers"](https://gmplib.org/~tege/division-paper.pdf) — single-limb reciprocal arithmetic where the limb size matters.
- [Polygon Plonky2 — Goldilocks field documentation](https://github.com/0xPolygonZero/plonky2) — production use of the Goldilocks prime `2⁶⁴ − 2³² + 1` and discussion of the 16-bit coefficient capacity that drives the fallback NTT's input split.

### Adjacent design spaces (not used here)

- [Booth encoding (signed-digit representation) — Wikipedia](https://en.wikipedia.org/wiki/Booth%27s_multiplication_algorithm) — signed-digit number systems, used in hardware multipliers and some cryptographic schemes.
- [Hensel coding and 2-adic representations — Wikipedia](https://en.wikipedia.org/wiki/P-adic_number#Hensel.27s_lemma) — alternative numeric representations explored in some specialized arithmetic libraries.

### This codebase

- `include/biginteger/common/Constants.h` — type aliases, base constants, `BIGMATH_LIMB_64` selector, and the `LimbBits` / `LimbBase` / `LimbMask` / `CurrentBase` helpers.
- `include/biginteger/BigInteger.h` — the `BigInteger` class definition; `Base()` returns `CurrentBase` (`Base2_64` sentinel by default, `Base2_32` under `-DBIGMATH_LIMB_64=0`).
- `include/biginteger/algorithms/Multiplication.h`, `include/biginteger/algorithms/Division.h` — top-level dispatchers; both take a `BaseT base` parameter.
- `include/biginteger/algorithms/multiplication/NTTMultiplicationCrt.h` — the default 3-prime CRT NTT with the 32-bit coefficient split (two per 64-bit limb).
- `include/biginteger/algorithms/multiplication/NTTMultiplication.h` — the single-prime Goldilocks fallback with the 16-bit input split (four per 64-bit limb; two per limb under the legacy layout).
- `include/biginteger/algorithms/multiplication/KaratsubaMultiplication.h` — Karatsuba with native 64-bit leaf; retains the historical pack/unpack hybrid leaf as the `Base2_32` fallback path.
- `include/biginteger/common/Parser.h` (implementation in `src/common/Parser.cpp`) — `Base10_18` parse chunking and `Base10_19` ToString chunking over `CurrentBase` storage.

### Companion documents

- [MULTIPLICATION.md](MULTIPLICATION.md) — multiplication algorithms built on this representation.
- [DIVISION.md](DIVISION.md) — division algorithms built on this representation.
