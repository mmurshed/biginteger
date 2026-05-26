# Number representation: base and endianness

A technical reference for the foundational choice underlying every arithmetic algorithm in this BigInteger library: **little-endian limbs in base 2³² stored as `uint64_t`**. This document defends the choice, shows how it works in practice, lays out the consequences across the codebase, summarizes optimizations enabled by the choice, and explores the alternatives that were considered and rejected.

---

## Table of contents

1. [Scope](#scope)
2. [The choice in one paragraph](#the-choice-in-one-paragraph)
3. [Type definitions](#type-definitions)
4. [Defending the choice](#defending-the-choice)
   - [Why base 2³² (not 2⁶⁴)](#why-base-2-not-264)
   - [Why little-endian](#why-little-endian)
   - [Why `uint64_t` storage for 32-bit values](#why-uint64_t-storage-for-32-bit-values)
   - [Why a runtime `BaseT` parameter](#why-a-runtime-baset-parameter)
5. [How it works in practice](#how-it-works-in-practice)
   - [Single-limb scalar multiply](#single-limb-scalar-multiply)
   - [Carry propagation in addition](#carry-propagation-in-addition)
   - [Borrow propagation in subtraction](#borrow-propagation-in-subtraction)
   - [Schoolbook multiplication inner loop](#schoolbook-multiplication-inner-loop)
   - [NTT input split](#ntt-input-split)
   - [Decimal I/O via base 10¹⁸ chunks](#decimal-io-via-base-1018-chunks)
6. [Benchmark consequences](#benchmark-consequences)
7. [Optimizations already implemented](#optimizations-already-implemented)
8. [Future opportunities](#future-opportunities)
9. [Explored but rejected](#explored-but-rejected)
10. [References](#references)

---

## Scope

This document explains why a single design decision — "limbs are 32-bit values in `uint64_t` slots, stored low-limb-first, with a runtime `BaseT` parameter" — propagates through the entire arithmetic stack and what the alternatives would cost.

Companion documents cover the algorithms built on this foundation:

- [MULTIPLICATION.md](MULTIPLICATION.md)
- [DIVISION.md](DIVISION.md)

This document is the prerequisite for both.

Assumed reader: a working C++ engineer familiar with `uint32_t` / `uint64_t` arithmetic, bit operations, and the basic shape of multi-precision algorithms.

---

## The choice in one paragraph

A `BigInteger` (`biginteger/BigInteger.h`) wraps a `std::vector<DataT>` where `DataT = uint64_t`, plus a sign bit. Each element of the vector holds a **32-bit limb value** (`0 ≤ limb < 2³²`). The upper 32 bits of each storage slot are headroom used during arithmetic to accumulate carries without spilling into a separate variable. Limbs are stored **little-endian**: `vec[0]` is the least-significant limb. The base is exposed as a runtime parameter `BaseT base` (in practice almost always `Base2_32 = 2³²`) so non-power-of-two bases (`Base10`, `Base100`, `Base100M`) can share the same algorithm code paths for I/O and testing.

```
                    most significant ──→
        ┌──────────┬──────────┬──────────┬──────────┐
   a =  │  a[3]    │  a[2]    │  a[1]    │  a[0]    │       sign = false  ⇒  a ≥ 0
        └──────────┴──────────┴──────────┴──────────┘
            ↑                              ↑
         high limb                      low limb (vec[0])

   each limb is uint64_t storage holding a value in [0, 2³²).
   the high 32 bits of each slot are headroom for carry accumulation.
```

---

## Type definitions

From `biginteger/common/Constants.h`:

| name | underlying | width | role |
|---|---|---|---|
| `DataT` | `uint64_t` | 64-bit storage, 32-bit value | one limb |
| `BaseT` | `int64_t` | constant `2³²` in production | the radix |
| `ULong` | `uint64_t` | 64-bit | accumulator for one product (fits 64×32 / 32×32 results) |
| `ULong128` | `unsigned __int128` | 128-bit | accumulator for 64-bit limb paths (NTT, hybrid basecase) |
| `Long` | `int64_t` | 64-bit signed | signed accumulator for subtraction/borrow |
| `Int` | `int32_t` | 32-bit signed | loop counters, signed limb diffs |
| `SizeT` | `uint32_t` | 32-bit | vector sizes, indices |

`BigInteger::Base()` returns `Base2_32`. Defined bases (most are I/O helpers):

```
Base2       = 2          ← used by tests, bit-level operations
Base10      = 10         ← parser inner-most digit chunk
Base100     = 100        ← legacy I/O helper
Base100M    = 100 000 000 (10⁸) ← legacy I/O helper
Base2_8     = 256        ← used by some tests
Base2_10    = 1024       ← used by some tests
Base2_16    = 65 536     ← NTT internal coefficient base
Base2_32    = 4 294 967 296 (2³²) ← PRODUCTION limb base
Base2_31    = 2 147 483 648 (2³¹) ← used in normalization
```

The production code path is `Base2_32` everywhere outside the NTT internal transform (which uses `Base2_16` for its convolution coefficients) and decimal I/O (which uses `Base10_18 = 10¹⁸`, a `uint64_t` chunk size for grouped decimal digits).

---

## Defending the choice

### Why base 2³² (not 2⁶⁴)

**The headline reason: every native multiplication of two limbs fits in a single 64-bit register.**

In base 2³², `limb × limb < 2³² · 2³² = 2⁶⁴`. The product is exactly a `uint64_t`. Adding two such products plus an in-place limb plus a carry all fit in 64 bits with one bit of slack:

```
   prod_max = (2³² − 1)² = 2⁶⁴ − 2³³ + 1
   adding (r[i+j] + carry) where each ≤ 2³² − 1:
   total_max = (2⁶⁴ − 2³³ + 1) + 2(2³² − 1) = 2⁶⁴ − 1.
   ── fits in uint64_t exactly, no overflow, no extra carry plumbing.
```

This is the central enabling property for compact and fast schoolbook multiplication inner loops:

```cpp
ULong prod = (ULong)a[j] * b[i] + r[i+j] + carry;
r[i+j] = (DataT)(prod & 0xFFFFFFFFULL);
carry  = prod >> 32;
```

Three instructions per partial product on ARM64: `UMULL` + `ADDS` + extract-high. No mispredicted branches, no register spilling beyond the carry chain.

In base 2⁶⁴, the same loop would need `__uint128_t` for `prod`. On ARM64 that's a `MUL` for the low half plus a `UMULH` for the high half (two instructions instead of one for the multiply), and the accumulator is two registers wide. The accumulator-plus-carry-plus-r[i+j] chain spans two registers and requires `ADCS` for carry-through-add. GMP handles this with hand-written assembly that interleaves multiple lanes and tunes instruction scheduling. Pure C++ compilers don't consistently produce equivalent code from `__uint128_t` operations.

So:

- **Base 2³² in pure C++** ≈ 1× of base-2⁶⁴-in-hand-tuned-asm performance for the inner partial product, *with* the simplification of single-register accumulators that the compiler reasons about well.
- **Base 2⁶⁴ in pure C++** ≈ 0.5× to 0.8× of the same baseline (compiler-dependent), recoverable only via assembly.

Combined with the further factor that base 2⁶⁴ means **half as many limbs to process**, the asymptotic comparison is:

```
   base 2³² pure C++:    one well-scheduled  UMULL per limb-pair
   base 2⁶⁴ pure C++:    two instructions (MUL+UMULH) per limb-pair,
                         but half as many limb-pairs
                         → net ~ same throughput, but worse register pressure
   base 2⁶⁴ asm (GMP):   one MUL+UMULH per limb-pair with optimal scheduling,
                         half as many limb-pairs
                         → net ~ 2× throughput vs base 2³² pure C++
```

The 2× factor is the structural gap to GMP that this codebase pays for staying in pure C++. Closing it requires either hand-written assembly or a fundamental shift to base 2⁶⁴ (analyzed in [Future opportunities](#future-opportunities)).

**A second reason: the 64-bit headroom enables the in-place scalar-multiply trick.**

Multiplying a vector of 32-bit limbs by a 32-bit scalar in place:

```cpp
for (j = 0; j < n; ++j) {
    a[j] = (a[j] & 0xFFFFFFFF) * scalar + carry;  // fits in upper bits of same uint64_t
    carry = a[j] >> 32;
    a[j] &= 0xFFFFFFFF;
}
```

The 32-bit headroom in each storage slot means the multiply result can live momentarily in the upper bits of the same slot before being separated into the new limb value and the carry-out. No temporary buffer required. This pattern shows up in `ClassicMultiplication::MultiplyTo` and the normalize step of `FastDivision`.

**A third reason: 32-bit-aligned bit packing for NTT.**

The NTT input split (see [NTT input split](#ntt-input-split) below) decomposes each 32-bit limb into exactly two 16-bit coefficients. The `2× 16-bit-per-limb` decomposition is alignment-trivial; the same algorithm at base 2⁶⁴ would split each limb into four 16-bit coefficients with the same coefficient count per byte of input, achieving nothing for the NTT but complicating the packing code.

### Why little-endian

The choice is **structural for arithmetic** (not a debate). In multi-precision arithmetic, carries propagate from low to high. Indexing a vector `a` with `a[0]` as the least significant limb means a single forward loop `for (i = 0; i < n; ++i)` handles carry propagation in addition, subtraction, and scalar multiplication.

Reverse this — store big-endian — and the same loop has to either decrement an index (`for (i = n-1; i >= 0; --i)`) or pre-reverse the buffer, neither of which the compiler vectorizes as cleanly. More importantly: extending a number (adding limbs at the high end during a carry that propagates past the current top) becomes `push_back(carry)` in little-endian and `insert(begin(), carry)` — an O(n) shift — in big-endian.

Little-endian also matches the on-disk and on-network conventions of most existing big-number libraries (GMP `mpn`, OpenSSL `BIGNUM` internals, Python's `_PyLong` digits, Java `BigInteger.mag[]` *internally is big-endian but exposed as such only because it pre-dates these tradeoffs*). Interop is easier when both sides agree.

A subtle benefit: **trimming leading zeros** (the canonical "is this number normalized" operation) is `while (vec.size() > 1 && vec.back() == 0) vec.pop_back();` — O(amortized 1) at the high end of the vector. The same operation in big-endian would erase elements from the front, an O(n) per-erase cost.

### Why `uint64_t` storage for 32-bit values

This is the "2× memory" tradeoff people notice first. It's deliberate. Three reasons:

1. **Carry headroom.** A multiply-accumulate result `a[i] * b[j] + r[i+j] + carry` always fits in 64 bits with the 32-bit limb representation. Storing limbs in `uint32_t` would require a separate `uint64_t` temporary for the accumulator at every loop iteration. Keeping limbs in `uint64_t` slots lets the accumulator live in the same register as the destination limb between iterations of the inner loop. The compiler's register allocator does noticeably better with this shape.

2. **Alignment for ULong128.** When the 64-bit hybrid Karatsuba leaf packs pairs of 32-bit limbs into `uint64_t` for the schoolbook (see [MULTIPLICATION.md §Classic schoolbook](MULTIPLICATION.md#classic-schoolbook)), the source layout `vec[2k] | vec[2k+1] << 32` reads two 64-bit slots and produces one 64-bit packed value with a single AArch64 `ORR` after the shift. If limbs were stored in `uint32_t`, the same pack would still be possible but each input read is now a 4-byte load — same data, no advantage, and the rest of the library's `ULong128` accumulator paths would need scattered casts.

3. **NTT coefficient extraction.** The Goldilocks NTT splits each 32-bit limb into two 16-bit coefficients via `(limb & 0xFFFF)` and `(limb >> 16)`. Both extractions on a `uint64_t` slot are one instruction. From `uint32_t`, the high-half extraction is still one instruction, but writes back into a `uint64_t` NTT working buffer require widening. The net is the same; the cleaner code path with `uint64_t` storage matters more.

The memory cost — 2× the storage of a hypothetical `uint32_t` representation — is rarely the binding constraint. A 1 000 000-digit BigInteger is ~417 KB of vector data either way; on M1 with 32 GB of RAM, the difference between 417 KB and 209 KB is invisible.

### Why a runtime `BaseT` parameter

Most modern big-number libraries hard-code their base at compile time. This library exposes `BaseT base` as a runtime parameter to every algorithm. Why?

**Historical**: the algorithm code was originally written to be base-agnostic for pedagogical clarity (the classical schoolbook algorithms in *TAOCP* are presented in terms of an abstract base B). The runtime parameter let the same code paths handle `Base10` for tutorial-level traces and `Base2_32` for production.

**Practical**: a small set of code paths still need non-power-of-two bases. Specifically:

- `ClassicDivision::DivModTo(vec, divisor, base)` is called from `ToStringLinearAppend` with `base = Base2_32` and `divisor = Base10_18`. The "base" parameter here is the **storage** base of the dividend vector (2³²), independent of the divisor.
- Test scaffolding builds small numbers in `Base10` or `Base16` for human readability before converting.
- Some legacy I/O helpers operate in `Base100M`.

The runtime parameter costs essentially nothing because of how the dispatching is structured: every hot inner loop branches on `base == Base2_32` *once at the top of the function* and falls into a fully-specialized fast path. The branch is highly predictable (always `true` in production), so the compiler keeps the fast path tight.

In the rare loops where the `base == Base2_32` check was historically *inside* the inner loop (the Karatsuba helpers `AddPtr` / `AddToPtr` / `SubtractFromPtr`), the 2026-05 rewrite hoisted it out — see [MULTIPLICATION.md §Optimizations](MULTIPLICATION.md#karatsuba-helper-rewrite-2026-05). After the hoist, the runtime parameter has no measurable cost in any production path.

If the parameter ever did become a bottleneck, the path forward is a templated `BaseT` (compile-time constant) — feasible because `Base2_32` is the only base any production caller passes. The templating would yield a non-`Base2_32` path with no specialization (currently the generic `% base` / `/ base` arithmetic), which is fine since those paths are only exercised by tests.

---

## How it works in practice

The base-and-endian choice cashes out as concrete instruction sequences. This section walks through the most performance-critical patterns.

### Single-limb scalar multiply

`a` is a vector of 32-bit limbs in base 2³². Multiply by scalar `b` (also 32-bit) in place:

```cpp
ULong carry = 0;
for (SizeT j = 0; j < n; ++j) {
    ULong prod = (ULong)a[j] * b + carry;
    a[j]  = (DataT)(prod & 0xFFFFFFFFULL);   // low 32 bits become new limb
    carry = prod >> 32;                       // high 32 bits become next carry
}
if (carry) a.push_back((DataT)carry);
```

Three ARM64 instructions per iteration: `UMULL` (32×32→64 multiply with add), an extract for the low half (free with the cast), and a shift for the carry. The `push_back` at the end handles the high-end growth case — only triggered when the original number's top limb times `b` overflowed into a new limb position. The vector's exponential growth strategy makes the amortized cost of `push_back` O(1).

In base 2⁶⁴ with `__uint128_t`, the same loop has `MUL + UMULH` for the multiply (two instructions instead of one for the multiply itself), and the accumulator is two registers wide. The carry propagation through the upper register requires `ADCS`. ~2× the instruction count per iteration.

### Carry propagation in addition

```cpp
ULong carry = 0;
for (SizeT i = 0; i < n; ++i) {
    ULong sum = (ULong)a[i] + b[i] + carry;
    r[i]  = (DataT)(sum & 0xFFFFFFFFULL);
    carry = sum >> 32;
}
if (carry) r.push_back((DataT)carry);
```

`a[i] + b[i]` is at most `2 · (2³² − 1) = 2³³ − 2`, plus a 1-bit carry, fits in `uint64_t` with room to spare. No overflow check needed. The high 32 bits of the sum are exactly the carry-out, by construction.

In base 2⁶⁴, `a[i] + b[i] + carry_in` overflows `uint64_t` for the worst case. The standard mitigation is either `__builtin_add_overflow` (one ARM64 `ADCS`-equivalent on supported platforms) or a `__uint128_t` accumulator. Both work; both add complexity that base 2³² avoids entirely.

### Borrow propagation in subtraction

```cpp
Long borrow = 0;
for (SizeT i = 0; i < n; ++i) {
    Long diff = (Long)a[i] - (Long)b[i] - borrow;
    if (diff < 0) { diff += (Long)1 << 32; borrow = 1; }
    else          { borrow = 0; }
    r[i] = (DataT)diff;
}
```

`a[i] − b[i]` is in `[−(2³² − 1), +(2³² − 1)]`, fits in `int64_t` (`Long`). The `if (diff < 0)` branch is well-predicted by modern branch predictors as long as the operand magnitudes are roughly similar.

In base 2⁶⁴, `a[i] − b[i]` for `b[i] > a[i]` would underflow `uint64_t`. The fix is `__builtin_sub_overflow` or a wider accumulator. Same overhead pattern as addition.

### Schoolbook multiplication inner loop

The fully specialized Base2_32 path inside `KaratsubaMultiplication::MultiplyClassicPtr` (after the 2026-05 hybrid rewrite — see [MULTIPLICATION.md §Classic schoolbook](MULTIPLICATION.md#classic-schoolbook)):

```cpp
// Pack pairs of 32-bit limbs into 64-bit values for inner-loop speed.
ULong a64[STACK_CAP];
ULong b64[STACK_CAP];
ULong r64[2 * STACK_CAP];

// Pack: a64[k] = a[2k] | (a[2k+1] << 32)
for (SizeT k = 0; k < na64; ++k) {
    ULong v = a[2*k];
    if (2*k + 1 < lenA) v |= (ULong)a[2*k + 1] << 32;
    a64[k] = v;
}
// (similarly for b64; zero r64)

// 64-bit schoolbook on the packed representation.
for (SizeT i = 0; i < nb64; ++i) {
    ULong bi = b64[i];
    if (bi == 0) continue;
    ULong carry = 0;
    for (SizeT j = 0; j < na64; ++j) {
        ULong128 prod = (ULong128)a64[j] * bi + r64[i+j] + carry;
        r64[i+j] = (ULong)prod;
        carry    = (ULong)(prod >> 64);
    }
    r64[i + na64] = carry;
}

// Unpack r64 back to 32-bit limbs.
for (SizeT k = 0; k < nr64; ++k) {
    if (2*k     < rLen) r[2*k]     = (DataT)(r64[k] & 0xFFFFFFFFULL);
    if (2*k + 1 < rLen) r[2*k + 1] = (DataT)(r64[k] >> 32);
}
```

This is the **base-promotion trick**: the inner multiplication runs in base 2⁶⁴ (recovering the per-multiply efficiency GMP gets natively), while the storage and the surrounding algorithms stay in base 2³² (preserving the headroom and the single-register-accumulator pattern that base 2³² enables everywhere else).

Throughput at the leaf:

```
   base 2³² scalar schoolbook:    1× (baseline)
   base 2⁶⁴ hybrid (this code):   ~2× (4× fewer multiplies, each ~2× cost)
   GMP hand-tuned ARM64 asm:      ~4× (same 4× fewer multiplies, optimal scheduling)
```

The hybrid recovers about half the gap to GMP at the leaf without touching the rest of the library. Stack buffers cover up to 64 packed limbs (= 128 source limbs), comfortably above the Karatsuba threshold of 48.

### NTT input split

The Goldilocks NTT uses prime `P = 2⁶⁴ − 2³² + 1`. Its convolution capacity constraint forces input coefficients to be at most 16 bits. The base-2³² limb representation maps to this with a straightforward bit-split:

```
   one 32-bit limb        →         two 16-bit NTT coefficients
                                    ┌──────────────────┬──────────────────┐
   limb (uint64 storage,    ──→     │   limb >> 16     │  limb & 0xFFFF   │
   value in [0, 2³²))               └──────────────────┴──────────────────┘
                                          (high half)        (low half)
```

In `NTTMultiplication.h::convolutionBase2_32`:

```cpp
vector<ULong> fa(n, 0);  // n is the NTT transform size (power of 2)
for (SizeT i = 0; i < A.size(); ++i) {
    SizeT j = i * 2;
    fa[j]     = A[i] & 0xFFFFULL;   // low 16-bit coefficient
    fa[j + 1] = A[i] >> 16;          // high 16-bit coefficient
}
```

This is *exactly* the operation the base choice was made for. With base 2⁶⁴ limbs, the same packing would produce four 16-bit coefficients per limb — same total coefficient count, more code, no algorithmic difference (the NTT working set is dominated by the coefficient count, not the source representation).

After the inverse transform, the carry propagation runs in base 2¹⁶ over the convolution output, then 16-bit-pair coefficients are re-assembled into 32-bit limbs:

```cpp
// carry propagation in base 2¹⁶
ULong carry = 0;
for (SizeT i = 0; i < c16.size(); ++i) {
    ULong total = c16[i] + carry;
    c16[i] = (DataT)(total & 0xFFFFULL);
    carry  = total >> 16;
}
// reassemble pairs into 32-bit limbs
for (SizeT i = 0; i < c16.size(); i += 2) {
    ULong low  = c16[i];
    ULong high = (i + 1 < c16.size()) ? c16[i + 1] : 0;
    c32.push_back((DataT)(low | (high << 16)));
}
```

The 16-bit ↔ 32-bit relationship is exact and lossless: every 32-bit limb is exactly two 16-bit coefficients, no padding, no edge cases. This wouldn't be true for base 2⁶⁴ (four 16-bit coefficients per limb, but consumers of the unpacked result still want 64-bit limbs, so the re-assembly is messier).

### Decimal I/O via base 10¹⁸ chunks

Decimal parsing reads input digits 18 at a time (because `10¹⁸ < 2⁶³`, the largest power of 10 that fits in a signed 64-bit chunk and still leaves room for a 32-bit multiplier without overflow):

```cpp
// in ParseUnsignedLinear:
while (pos <= end) {
    ULong chunk = ParseChunk(num, pos, Base10_18_Zeroes);  // 18 ASCII digits → uint64
    pos += Base10_18_Zeroes;
    ClassicMultiplication::MultiplyTo(r, Base10_18, Base2_32);  // r ← r * 10¹⁸
    AddTo(r, chunk, Base2_32);                                   // r ← r + chunk
}
```

`ClassicMultiplication::MultiplyTo(r, Base10_18, Base2_32)` multiplies each 32-bit limb by `10¹⁸` (a 60-bit constant). The result of one partial product is `(2³² − 1) · 10¹⁸ ≈ 2⁹² + small`, which fits in `ULong128`. The `Base2_32` parameter governs the storage base for `r`: low 32 bits become new limb, high bits become carry into the next position.

Decimal output is the inverse: `ClassicDivision::DivModTo(r, Base10_18, Base2_32)` divides `r` by `10¹⁸` in place and returns the remainder (the next 18 decimal digits to format). The per-limb step uses a Möller-Granlund "div2by1" reciprocal divide (`UMULH` + 128-bit add + ≤2 fixups), with the reciprocal precomputed once per call from the invariant `10¹⁸` divisor.

Both directions exploit the fact that `10¹⁸` fits in a 64-bit register while leaving the 32-bit limb storage representation untouched. The base 2³² storage and the base 10¹⁸ I/O chunking are independent design choices that happen to compose cleanly.

---

## Benchmark consequences

Most of the bench numbers in [MULTIPLICATION.md](MULTIPLICATION.md#benchmark-results-vs-gmp) and [DIVISION.md](DIVISION.md#benchmark-results-vs-gmp) trace back to the base choice. Three representative slices to illustrate the link:

### Scalar inner-loop throughput (the schoolbook leaf, post-hybrid)

| op | BigMath ms | GMP ms | BM / GMP |
|---|---:|---:|---:|
| mul 1 000 × 1 000 digits | 0.003 | 0.001 | 3.5 × |
| mul 5 000 × 5 000 digits | 0.043 | 0.012 | 3.5 × |

Both cases route to the Karatsuba leaf where the 64-bit hybrid runs. The 3.5× residual gap factors as ~2× (pure C++ vs hand-tuned ARM64 assembly in GMP's `mpn_mul_basecase`) × ~1.75× (vector overhead, dispatch, BigInteger object construction).

### NTT-bound throughput (the input-split is the binding constraint)

| op | BigMath ms | GMP ms | BM / GMP |
|---|---:|---:|---:|
| mul 1 000 000 × 1 000 000 digits | 37.5 | 8.9 | 4.2 × |

97% of this time is in NTT butterflies. The base-2³² limb storage **does not enter the cost** here — the binding constraint is the Goldilocks prime's 16-bit coefficient capacity, which is identical for any limb size ≥ 32 bits. Switching to base 2⁶⁴ would leave this 4.2× ratio unchanged. The gap is from compiler vs assembly in the butterfly inner loop.

### Skewed division (Newton iteration, dominated by NTT internal mults)

| op | BigMath ms | GMP ms | BM / GMP |
|---|---:|---:|---:|
| div 200 000 / 50 000 digits | 22 | 1.7 | 14 × |
| div 500 000 / 100 000 digits | 55 | 4.6 | 12 × |

Same NTT floor as large multiplication, but compounded by Newton iteration overhead (3 large mults per refinement, blockwise dispatch). Base choice has no effect at these sizes.

**Summary.** The base-2³² choice costs ~2× at the schoolbook leaf relative to GMP. The 64-bit hybrid Karatsuba leaf recovers about half of that. The remaining gap is pure C++ vs hand-tuned ARM64 assembly, not base-related. At NTT sizes, base is irrelevant.

---

## Optimizations already implemented

A summary of base-related optimizations that landed and stuck. (Each is covered in detail in the algorithm-specific documents.)

### 64-bit hybrid leaf (the headline win)

`KaratsubaMultiplication::MultiplyClassicPtr` Base2_32 path packs pairs of 32-bit limbs into 64-bit values, runs schoolbook in 64-bit limb space (`__uint128_t` accumulator), and unpacks. ~2× over scalar 32-bit schoolbook at the leaf, which translated to **mul 5 000 × 5 000 going from 9.3× to 3.5× vs GMP**. Discussed in [MULTIPLICATION.md §Classic schoolbook](MULTIPLICATION.md#classic-schoolbook).

### Hoisted `base == Base2_32` branch in inner loops

`AddPtr`, `AddToPtr`, `SubtractFromPtr` in `KaratsubaMultiplication.h` had the runtime base check pulled out of the inner loop and replaced with two parallel specialized paths. The branch is highly predictable, but the loop body without the branch is one fewer instruction per iteration. Modest 0–4% net wins. Discussed in [MULTIPLICATION.md §Karatsuba helper rewrite](MULTIPLICATION.md#karatsuba-helper-rewrite-2026-05).

### NTT direct-fill of 16-bit coefficient buffers

`NTTMultiplication::convolutionBase2_32` builds the NTT input vectors directly with the 16-bit split, skipping intermediate `a16` / `b16` temporary vectors. The base-2³² limb format makes this a one-line packing per limb. Saves allocation + copy overhead proportional to input size.

### `FastDivision` Base2_32 specializations

Normalize, scalar divisor, and remainder unnormalize have Base2_32-specific code paths using `ULong128` accumulators and `__builtin_umul128`-equivalent intrinsics. ~1.04–2× wins on those steps, depending on operand size. Covered in [DIVISION.md §FastDivision](DIVISION.md#fast-division-knuth-algorithm-d).

### Decimal D&C chain using `NewtonDivision::Divider` over base 2³²

`Parser.h::BuildDecimalDcChain` builds a chain of `Divider` instances, each holding a `Pow10(d)` value (in base 2³² limbs) and its precomputed reciprocal. This is the foundation of the 8.4× ToString speedup at 100k digits. The base 2³² limb representation is what makes the underlying `ApproxReciprocal` and `DivideChunk` operations efficient — each operates on `ULong128` accumulators over base 2³² limbs.

---

## Future opportunities

Ranked by expected ROI per unit of effort.

### Templated `BaseT` (compile-time constant)

Replace the runtime `BaseT base` parameter with a template parameter `template <BaseT base>`. For production callers using `Base2_32`, the compiler would inline the constant and eliminate any residual conditional branching that the hoisted-branch approach doesn't fully cover.

Effort: moderate (touches every algorithm signature). Risk: low — the only callers passing non-`Base2_32` are tests, which can be templated trivially or pinned to a runtime variant for their specific bases. Estimated win: 1–3% across all base-sensitive paths. Net: marginal, but cumulative across the whole library.

### Full 64-bit limb refactor

Discussed in detail in [MULTIPLICATION.md §Future opportunities](MULTIPLICATION.md#future-opportunities) and [DIVISION.md §Future opportunities](DIVISION.md#future-opportunities). For the base/endian story specifically:

- Switching from `uint64_t storing 32-bit values` to `uint64_t storing 64-bit values` would require replacing every `ULong` accumulator (1 product fits) with `ULong128` (1 product fits), every carry-handling pattern (32-bit add fits in 64-bit) with `__builtin_*_overflow` (64-bit add needs explicit overflow detection).
- The NTT input split would change from `2× 16-bit per 32-bit limb` to `4× 16-bit per 64-bit limb`, with no net change to NTT performance (coefficient count is the same).
- The 64-bit hybrid leaf could be removed (the schoolbook would just *be* in base 2⁶⁴), saving the pack/unpack overhead — small win (~1.1×).

Estimated cost: 1–2 weeks, 1000–2000 lines. Estimated wins:

| op | now | after | source of gain |
|---|---|---|---|
| mul 1 000 × 1 000 | 3.5 × | ~3 × | Karatsuba leaf already hybrid — marginal |
| mul 1 000 000 × 1 000 000 | 4.2 × | 4.2 × | **zero** — NTT input still 16-bit |
| div 200 000 / 50 000 | 14 × | ~7–9 × | Newton inner mults run in base 2⁶⁴ scalar arithmetic |
| parse 100 000 | 7.1 × | ~5 × | linear parser inner loop |
| ToString 100 000 | 17.7 × | ~12–14 × | linear formatter inner loop + Newton chain |

NTT-bound operations (large mult, large div) gain nothing — the Goldilocks prime caps coefficient width at 16 bits regardless of source limb size. The structural floor of 4× vs GMP at NTT sizes survives.

### NEON-vectorized 32-bit limb arithmetic

ARM NEON has 32×32→64 multiply-with-accumulate (`UMLAL` on 4-lane vectors, `UMLAL2` for the upper half). Vectorizing the schoolbook inner loop with 4 lanes could in principle give 4× throughput on the multiply, but the **carry propagation chain is fundamentally serial** — each partial product writes a limb that depends on the previous limb's carry-out. NEON helps only if the carries can be decoupled.

One scheme: compute partial products lane-parallel into a wider accumulator (e.g., 64-bit per lane), then propagate carries serially at the end. Estimated win: 1.5–2× on the schoolbook leaf, with careful tuning. Effort: ARM-specific intrinsics path, requires a scalar fallback, ~200 lines plus careful debugging.

Note: the 64-bit hybrid already captures most of the structural win at the leaf. NEON on top of the hybrid would only help if the leaf is still the dominant cost (which post-hybrid it isn't — see [MULTIPLICATION.md §Optimizations](MULTIPLICATION.md#optimizations-already-implemented) and the post-(b) profile analysis).

### Aligned-storage hint to the allocator

`std::vector<DataT>` doesn't guarantee 16-byte alignment of its data. Aligned loads are slightly faster on M1 and crucial for NEON. A custom allocator or `alignas`-declared inline storage for small numbers could give 2–5% on hot inner loops. Marginal, but cheap.

---

## Explored but rejected

### Base 10 / Base 10ⁿ as the storage representation

**Rejected.** Tempting because decimal I/O becomes trivial (no parse/format conversion overhead). Killed by every other operation: arithmetic in base 10 requires `% 10` and `/ 10` per limb step — both extremely slow compared to the bit-mask and shift that base 2³² uses. The classical schoolbook in pure base 10 would be 10–20× slower than the base 2³² implementation. Real implementations using decimal storage (e.g. `_decimal128` formats) are designed for fixed-precision applications where I/O cost dominates; arbitrary-precision libraries universally use power-of-two bases.

The library does internally use base 10¹⁸ as a **chunking convention for I/O**, but not as a storage representation. The chunks are immediately multiplied into base 2³² limbs during parsing.

### Base 2⁶³ (one bit reserved for sign)

**Rejected.** Some libraries reserve the top bit of each limb for the sign of the partial product, eliminating the separate sign field. Costs: every arithmetic operation must mask the sign bit, and the per-limb effective range shrinks. The savings (one boolean per `BigInteger`) are negligible compared to the loss of clean 64×64-or-32×32 arithmetic. The sign-as-separate-bool design is universally preferred in modern libraries.

### Signed-digit (non-adjacent form) representation

**Rejected.** Signed-digit representations (digits in `{−B/2, ..., +B/2}`) make addition carry-free at the cost of a non-canonical representation (multiple ways to represent the same number) and conversion overhead at boundaries (I/O, comparison). Useful in some cryptographic applications (Montgomery ladders, ECC scalar multiplication) but the library's use case (general-purpose big-number arithmetic) doesn't benefit. Conversion costs would dominate.

### Big-endian limb storage

**Rejected.** See [Why little-endian](#why-little-endian). The carry-propagation direction (low → high) matches the array iteration direction; trim-leading-zeros is O(1) at the high end of a `std::vector`. Big-endian inverts both, with no compensating benefit.

### Compile-time `BaseT` only (templated, no runtime parameter)

**Considered but not yet executed.** Discussed in [Future opportunities](#future-opportunities). The runtime parameter has no measurable cost in the current code (the `base == Base2_32` checks are hoisted out of inner loops), so the templating refactor is not currently prioritized. Could be done if testing infrastructure ever wants to template-instantiate algorithm code paths for fuzz testing of non-Base2_32 paths.

### Variable-length limbs (e.g., 28-bit for digit-recurrence FFT)

**Not considered seriously.** Some FFT-multiplication schemes use 28-bit or 23-bit limbs to fit results in a `double` mantissa for IEEE 754 floating-point NTT. The Goldilocks integer-NTT path used here doesn't benefit; the limb size is decoupled from the transform coefficient size by the 16-bit input split.

### Native `__uint128_t` as the limb type

**Rejected.** `__uint128_t` is a GCC/Clang extension, not part of standard C++. Using it as `DataT` would lock out MSVC (where `__int128` doesn't exist) and complicate the `std::vector` allocator path. The `ULong128` accumulator is used in narrow scopes (NTT, hybrid leaf, FastDivision normalize) where the platform requirement is acceptable, but pushing it to the limb type would force the whole library to be GCC/Clang-only.

The library is already effectively GCC/Clang-only because of `__builtin_*_overflow` usage in NTT — but that's a much smaller surface than the limb type, and could in principle be replaced with portable C++23 `<stdckdint.h>` operations.

---

## References

### Multi-precision arithmetic foundations

- Knuth, D. E. *The Art of Computer Programming, Vol. 2: Seminumerical Algorithms*, §4.3 — the canonical treatment of positional notation, limb-based arithmetic, and base choice.
- [Brent, R. P. and Zimmermann, P. — *Modern Computer Arithmetic* (Cambridge, 2010)](https://members.loria.fr/PZimmermann/mca/pub226.html) — modern reference; Chapter 1 covers integer representation including the base/limb/endianness tradeoffs.

### Hardware multiplication primitives

- [ARM Architecture Reference Manual: `UMULL` / `UMULH` / `UMADDL`](https://developer.arm.com/documentation/ddi0596/2021-12/Base-Instructions/UMULL--Unsigned-Multiply-Long-) — the instructions underlying base 2³² limb arithmetic on ARM64.
- [Intel 64 and IA-32 Architectures Software Developer's Manual, Vol. 2: `MUL` / `MULX`](https://www.intel.com/content/www/us/en/developer/articles/technical/intel-sdm.html) — x86 multiplication primitives.

### Reference implementations and conventions

- [GMP — The GNU Multiple Precision Arithmetic Library](https://gmplib.org/) — uses base 2⁶⁴ (configurable to 2³² at build time on 32-bit platforms). Little-endian limbs.
- [GMP manual: "Low-level Functions" (`mpn_*`)](https://gmplib.org/manual/Low_002dlevel-Functions) — documentation of the `mp_limb_t` and `mp_size_t` types and the calling conventions, including the little-endian convention.
- [OpenSSL `BIGNUM` internals (`bn.h`)](https://github.com/openssl/openssl/blob/master/include/openssl/bn.h) — `BN_ULONG` is either 32-bit or 64-bit depending on the platform; little-endian limbs.
- [Python `_PyLong` digit representation](https://github.com/python/cpython/blob/main/Include/cpython/longintrepr.h) — uses 30-bit "digits" stored in 32-bit words on 64-bit platforms, leaving 2 bits of carry headroom. Conceptually similar tradeoff to this library's 32-in-64 design.
- [Java `BigInteger.mag[]`](https://docs.oracle.com/en/java/javase/17/docs/api/java.base/java/math/BigInteger.html) — `int[]` of 32-bit magnitudes, stored big-endian (`mag[0]` is most significant). The big-endian convention is a holdover from early Java's I/O orientation; modern libraries universally pick little-endian.

### Specific design references

- [Möller, N. and Granlund, T. — "Improved Division by Invariant Integers"](https://gmplib.org/~tege/division-paper.pdf) — single-limb reciprocal arithmetic where the limb size matters.
- [Polygon Plonky2 — Goldilocks field documentation](https://github.com/0xPolygonZero/plonky2) — production use of the Goldilocks prime `2⁶⁴ − 2³² + 1` and discussion of the 16-bit coefficient capacity that drives the input split.

### Adjacent design spaces (not used here)

- [Booth encoding (signed-digit representation) — Wikipedia](https://en.wikipedia.org/wiki/Booth%27s_multiplication_algorithm) — signed-digit number systems, used in hardware multipliers and some cryptographic schemes.
- [Hensel coding and 2-adic representations — Wikipedia](https://en.wikipedia.org/wiki/P-adic_number#Hensel.27s_lemma) — alternative numeric representations explored in some specialized arithmetic libraries.

### This codebase

- `biginteger/common/Constants.h` — type aliases and base constants.
- `biginteger/BigInteger.h` — the `BigInteger` class definition; `Base()` returns `Base2_32`.
- `biginteger/algorithms/Multiplication.h`, `biginteger/algorithms/Division.h` — top-level dispatchers; both take `BaseT base` parameter.
- `biginteger/algorithms/multiplication/NTTMultiplication.h` — the 16-bit NTT input split for base 2³² limbs.
- `biginteger/algorithms/multiplication/KaratsubaMultiplication.h` — the 64-bit hybrid leaf that bridges base 2³² storage with base 2⁶⁴ inner-loop arithmetic.
- `biginteger/common/Parser.h` — `Base10_18` decimal I/O chunking over base 2³² storage.

### Companion documents

- [MULTIPLICATION.md](MULTIPLICATION.md) — multiplication algorithms built on this representation.
- [DIVISION.md](DIVISION.md) — division algorithms built on this representation.
