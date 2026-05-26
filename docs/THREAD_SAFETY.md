# Thread safety

## Summary

BigMath is **thread-safe for concurrent use of distinct objects from distinct threads**. Multiple threads may simultaneously construct, mutate, and destroy their own `BigInteger`/`BigDecimal` values without coordination. Multiple threads may also concurrently call `const` methods on a *shared* `BigInteger`, `BigDecimal`, or `NewtonDivision::Divider` instance.

Concurrent mutation of the *same* object from multiple threads is **not** supported (standard data-race semantics). Wrap shared mutable instances in your own mutex if needed.

## What's shared, what's not

### Per-thread caches (already isolated)

The library uses four `static thread_local` caches that are private to each thread. No coordination is needed; each thread pays a first-touch cost.

| cache | file | what it holds |
|---|---|---|
| `NTTCore::BuildPlan::cache` | `include/biginteger/algorithms/multiplication/NTTCore.h` | NTT plans (size → twiddles + bit-reversal table) |
| `NTTMultiplication::Multiply::fa,fb` | `include/biginteger/algorithms/multiplication/NTTMultiplication.h` | coefficient working buffers |
| `NTTSquare::Square::fa` | `include/biginteger/algorithms/multiplication/NTTSquare.h` | coefficient working buffer |
| `Pow10::cache` | `src/common/Parser.cpp` | memoized powers of 10 in current-base limbs |

These caches grow monotonically (NTT plan cache is keyed by transform size; Pow10 cache is keyed by exponent). They are never invalidated mid-process and never written from outside the owning thread.

For a thread-pool worker pattern, warm each pool thread by issuing one representative NTT-bound mult and one `Pow10(d)` call from each worker at startup. Otherwise the first call from each worker pays the cache-fill cost.

### Read-only namespace constants

Dispatch thresholds (`NTT_MULTIPLICATION_THRESHOLD`, `NEWTON_MEDIUM_B`, etc.) are `const SizeT` defined in `src/algorithms/*.cpp`. Initialized once at program start; never modified.

### `NewtonDivision::Divider` (and `ReciprocalDivision::Divider`)

The `Divider` class precomputes a reciprocal for a fixed divisor. Its `DivideAndRemainder` and `Divide` methods are `const`; the precomputed state is not mutated after construction. **One `Divider` instance can be safely shared across threads** as long as no thread is concurrently constructing it.

Typical pattern:

```cpp
// One Divider, many concurrent divides on different numerators.
ReciprocalDivision::Divider divider(b, BigInteger::Base());

#pragma omp parallel for
for (size_t i = 0; i < numerators.size(); ++i) {
  auto qr = divider.DivideAndRemainder(numerators[i]);  // safe
  ...
}
```

### `BigInteger` / `BigDecimal` instances

Each instance owns its `std::vector<DataT>` storage. Standard C++ rules apply: concurrent reads of a single instance are fine, concurrent writes (or read+write) are a data race. Operators (`+`, `-`, `*`, `/`, etc.) produce new instances and don't mutate operands, so chains like `a + b * c` from a single thread are always race-free regardless of whether `a`, `b`, `c` are shared with other threads (provided the others aren't writing).

## What is NOT thread-safe

- **Concurrent mutation of the same `BigInteger`/`BigDecimal` instance.** Includes `+=`, `-=`, `*=`, `/=`, `AddTo`, `SubtractFrom`, assignment.
- **Concurrent calls to `Divider` constructor on the same memory location.** Construction does heavy precompute work; finish construction before sharing.
- **Concurrent calls to `BigInteger::Base()` that race with a hypothetical future setter.** Currently `Base()` is a static `constexpr` and immutable; this caveat is forward-looking.

## Internal parallelism (`BIGMATH_USE_THREADS`, default on since 2026-05)

A small thread pool is linked in when `BIGMATH_USE_THREADS=1` (the default). Used by the CRT NTT path to dispatch 6 forward transforms + 3 inverse transforms as batched work units.

- **Pool size**: `min(hardware_concurrency(), BIGMATH_MAX_THREADS=8)` — on shared-L2 architectures (M1 Max etc.), going beyond 8 cores hits L2 cache pressure on NTT working sets (~512 KB at 32k-coeff transforms).
- **Linkage**: pool implementation lives in `src/common/Parallel.cpp`. Public headers stay free of `<thread>` so consumers don't pick up pthread unconditionally.
- **First-touch cost**: the `static thread_local` caches above remain per-thread. Each pool worker fills its own NTT plan / Pow10 caches on first use. For latency-sensitive workloads, warm the pool with one large `Multiply` from each worker at startup.
- **Caller participation**: the calling thread runs the first work chunk itself, so effective parallelism = pool size (not pool size + 1).
- **The user-facing thread safety guarantees above are unchanged.** Internal parallelism is an implementation detail of single operation calls, not a change in the concurrency model.

Opt-out: `-DBIGMATH_USE_THREADS=0` reverts to fully serial code paths and drops the pthread linkage. Useful for embedded targets or strict-header-only consumers.

See [`DIVISION.md` §Multithreaded NTT](DIVISION.md#multithreaded-ntt-2026-05-prs-32-38-default-since-39) for the design rationale and measured speedups.
