# Thread safety

## Summary

BigMath is **thread-safe for concurrent use of distinct objects from distinct threads**. Multiple threads may simultaneously construct, mutate, and destroy their own `BigInteger`/`BigDecimal` values without coordination. Multiple threads may also concurrently call `const` methods on a *shared* `BigInteger`, `BigDecimal`, or `NewtonDivision::Divider` instance.

Concurrent mutation of the *same* object from multiple threads is **not** supported (standard data-race semantics). Wrap shared mutable instances in your own mutex if needed.

## What's shared, what's not

### Per-thread caches (already isolated)

All of the library's internal caches and scratch buffers are `static thread_local` — private to each thread. No coordination is needed; each thread pays a first-touch cost. Grouped by subsystem:

| cache group | file | what it holds |
|---|---|---|
| NTT plan caches | `include/biginteger/algorithms/multiplication/NTTCore.h` (`BuildPlan::cache`), `include/biginteger/algorithms/multiplication/NTTMultiplicationCrt.h` (per-prime `Plan` cache, factorization cache) | NTT plans (size → twiddles + bit-reversal table) |
| NTT working buffers | `include/biginteger/algorithms/multiplication/NTTMultiplication.h` (`fa,fb`), `NTTSquare.h` (`fa`), `NTTMultiplicationCrt.h` (`fa1..fa3`/`fb1..fb3` coefficient buffers, pack buffers, MFA tile buffers and `mfaScratch`/`cycScratch`) | coefficient and tile working buffers |
| Newton division scratch | `include/biginteger/algorithms/division/NewtonDivision.h` (`Scratch()` → `ScratchBuffers`) | reusable temporaries for the reciprocal iteration and divide steps |
| `Pow10::cache` | `src/common/Parser.cpp` | memoized powers of 10 in current-base limbs |
| ToString divider-chain cache | `src/common/Parser.cpp` (`GetDecimalDcChain`) | per-size chains of `10^k` values + precomputed Newton reciprocals |
| `BigDecimal` `Pow10Bi` cache | `bigdecimal/BigDecimal.cpp` | `BigInteger` wrappers over cached powers of 10 |

The point is per-thread isolation, not the exact inventory — any new cache or scratch buffer added to the library follows the same `static thread_local` pattern. The map-like caches grow monotonically (NTT plan caches are keyed by transform size; Pow10 and divider-chain caches by digit count). They are never invalidated mid-process and never written from outside the owning thread.

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

A small thread pool is linked in when `BIGMATH_USE_THREADS=1` (the default). Its users:

- the **non-MFA CRT NTT path** dispatches 6 forward transforms + 3 inverse transforms as batched work units;
- the **fused MFA pipeline** (post-PR-#107, transform length ≥ 2^20) instead dispatches row-chunked `ParallelDo(6)` phases per fused stage, keeping 6–12 work units in flight;
- the **parse and ToString D&C fan-outs** (PRs #118/#119) dispatch one flat `ParallelDo` over 2³ subtree work items above ~100k digits.

`ParallelDo` is reentrant-safe: a thread-local nesting guard (`tl_chunkDepth` in `src/common/Parallel.cpp`) forces any `ParallelDo` issued from inside a chunk body to run inline serially, so subtree workers that reach Newton/NTT internals (which themselves call `ParallelDo`) cannot corrupt the pool's single work slot.

- **Pool size**: `min(hardware_concurrency(), BIGMATH_MAX_THREADS=8)` — on shared-L2 architectures (M1 Max etc.), going beyond 8 cores hits L2 cache pressure on NTT working sets (~512 KB at 32k-coeff transforms).
- **Linkage**: pool implementation lives in `src/common/Parallel.cpp`. Public headers stay free of `<thread>` so consumers don't pick up pthread unconditionally.
- **First-touch cost**: the `static thread_local` caches above remain per-thread. Each pool worker fills its own NTT plan / Pow10 caches on first use. For latency-sensitive workloads, warm the pool with one large `Multiply` from each worker at startup.
- **Caller participation**: the calling thread runs the first work chunk itself, so effective parallelism = pool size (not pool size + 1).
- **The user-facing thread safety guarantees above are unchanged.** Internal parallelism is an implementation detail of single operation calls, not a change in the concurrency model.

Opt-out: `-DBIGMATH_USE_THREADS=0` reverts to fully serial code paths and drops the pthread linkage. Useful for embedded targets or strict-header-only consumers.

See [`DIVISION.md` §Multithreaded NTT](DIVISION.md#multithreaded-ntt-2026-05-prs-32-38-default-since-39) for the design rationale and measured speedups.
