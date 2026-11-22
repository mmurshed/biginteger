# Thread safety

## Summary

BigMath is **thread-safe for concurrent use of distinct objects from distinct threads**. Multiple threads may simultaneously construct, mutate, and destroy their own `BigInteger`/`BigDecimal` values without coordination. Multiple threads may also concurrently call `const` methods on a *shared* `BigInteger`, `BigDecimal`, or `NewtonDivision::Divider` instance.

Concurrent mutation of the *same* object from multiple threads is **not** supported (standard data-race semantics). Wrap shared mutable instances in your own mutex if needed.

## What's shared, what's not

### Transform tables and scratch buffers

In threaded builds, callers share immutable CRT leaf plans and matrix roots for each prime and transform length. First use publishes the table under a per-length mutex. Each participating thread keeps an owning reference, so its warm lookup needs no lock. The table is released when its last owning thread exits. Mutable arithmetic buffers remain private to each thread.

| cache group | file | what it holds |
|---|---|---|
| NTT plan caches | `include/biginteger/algorithms/multiplication/NTTCore.h` (`BuildPlan::cache`), `include/biginteger/algorithms/multiplication/NTTMultiplicationCrt.h` (per-prime `Plan` cache, factorization cache) | NTT plans (size → twiddles + bit-reversal table) |
| NTT working buffers | `include/biginteger/algorithms/multiplication/NTTMultiplication.h` (`fa,fb`), `NTTSquare.h` (`fa`), `NTTMultiplicationCrt.h` (`fa1..fa3`/`fb1..fb3` coefficient buffers, pack buffers, MFA tile buffers and `mfaScratch`/`cycScratch`) | coefficient and tile working buffers |
| Newton division scratch | `include/biginteger/algorithms/division/NewtonDivision.h` (`Scratch()` → `ScratchBuffers`) | reusable temporaries for the reciprocal iteration and divide steps |
| `Pow10::cache` | `src/common/Parser.cpp` | memoized powers of 10 in current-base limbs |
| ToString divider-chain cache | `src/common/Parser.cpp` (`GetDecimalDcChain`) | per-size chains of `10^k` values + precomputed Newton reciprocals |
| `BigDecimal` `Pow10Bi` cache | `bigdecimal/BigDecimal.cpp` | `BigInteger` wrappers over cached powers of 10 |

Linear and cyclic CRT products on one thread reuse the same residue and transpose buffers. Each operation reconstructs its result into an owned limb vector before the next operation uses those buffers. Decimal powers and divider chains remain cached per thread.

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

The pool has one work slot, and exactly one caller owns it at a time (`dispatchMutex`, taken with `try_lock`). A second user thread that dispatches while the slot is busy runs its work **inline on its own thread** rather than waiting — so concurrent pool-backed operations from distinct user threads are safe, and under contention each extra caller simply degrades to the serial path. (Before 2026-11-20 there was no such ownership: a second dispatcher overwrote the job in flight, which lost chunks and corrupted the completion count — concurrent NTT-sized multiplies from two user threads hung. `tests/unit/test_parallel_pool.cpp` pins the fix.)

- **Dispatch model**: the caller and the workers claim slices of a job from a shared atomic counter tagged with the job's generation. The caller therefore never waits for a parked worker to wake up when it could run the slice itself, and a worker that arrives after a job finished cannot claim against the next one with stale state.
- **Lock-free job pickup** (2026-11-21): an awake worker takes a job without touching a mutex. The slot owner publishes in a fixed order — close the claim word under the new generation tag, write the job description (atomics), open the claim word at slice 0, publish the generation — and a worker reads the generation, then the description, then claims with that generation; a description that is stale or torn belongs to a generation whose claims have already been closed, so it is never acted on. The pool mutex is now used for parking only: a worker announces that it is about to park (it drops `awake` before its last look at the generation) and the dispatcher notifies only when it sees a worker missing. Before this, every worker took the pool mutex on every job to read four words; seven of them arriving at once serialized on it and the losers slept in the kernel, which cost ~13 µs per warm dispatch and made NTT timings depend on unrelated code layout (the same multiply ran 15–25% slower in one build than another). Measured and rejected: per-worker condition variables for the cold wake-up — the stagger between workers after a wake is kernel wake latency, not a lock chain.
- **Pool temperature**: `ParallelWarm()` reports whether enough workers are awake for a dispatch to be picked up immediately, and `ParallelPrewarm()` wakes parked workers without giving them work (the caller wakes one, which relays to the rest). Mid-sized multiplies and squares use these through `ParallelWarmOrPrewarm()` / `ParallelNoteColdDone()`: with a parked pool a small NTT costs about its serial time (~52 µs at 320–512 limbs, against 24 µs warm), so the dispatcher takes Karatsuba when cold and the NTT when warm. A thread that runs two such cold operations within `BIGMATH_POOL_SPIN_US` of each other prewarms the pool, so a burst switches to the NTT from its third operation; an isolated operation wakes nothing. The choice affects timing only — both algorithms are exact.
- **Spin before park**: an idle worker (and a caller waiting for its last slices) polls for `BIGMATH_POOL_SPIN_US` microseconds (default 100) before parking on a condition variable. Waking a parked thread costs 20–30 µs on macOS — more than a whole small NTT — while one multiply issues its dispatches tens of microseconds apart. Measured `ParallelDo(6, empty)`: 22.5 µs → 2.9 µs back-to-back. Build with `-DBIGMATH_POOL_SPIN_US=0` to park immediately (lowest idle CPU, old latency).

- **Pool size**: `min(hardware_concurrency(), BIGMATH_MAX_THREADS=8)` — on shared-L2 architectures (M1 Max etc.), going beyond 8 cores hits L2 cache pressure on NTT working sets (~512 KB at 32k-coeff transforms).
- **Linkage**: pool implementation lives in `src/common/Parallel.cpp`. Public headers stay free of `<thread>` so consumers don't pick up pthread unconditionally.
- **First-touch cost**. Pool workers share constructed CRT tables. A worker still allocates its own mutable workspace and decimal power cache as needed.
- **Caller participation**: the calling thread runs the first work chunk itself, so effective parallelism = pool size (not pool size + 1).
- **The user-facing thread safety guarantees above are unchanged.** Internal parallelism is an implementation detail of single operation calls, not a change in the concurrency model.

Opt-out: `-DBIGMATH_USE_THREADS=0` reverts to fully serial code paths and drops the pthread linkage. Useful for embedded targets or strict-header-only consumers.

See [`DIVISION.md` §Multithreaded NTT](DIVISION.md#multithreaded-ntt-2026-05-prs-32-38-default-since-39) for the design rationale and measured speedups.
