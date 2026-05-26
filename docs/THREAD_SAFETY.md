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

## Future: opt-in internal parallelism

A `BIGMATH_USE_THREADS` build flag is reserved for adding internal thread-pool-driven parallelism inside individual operations (e.g. parallel NTT butterflies, parallel `Forward(fa)`+`Forward(fb)`). When this lands:

- The flag will default to **off**. Single-threaded users see zero overhead.
- A small `BigMath::ThreadPool` will be linked in via `src/`. Public headers stay free of `<thread>` to avoid pulling pthread linkage into every consumer.
- Pool size defaults to `min(hardware_concurrency(), 8)` — beyond 8 cores on shared-L2 architectures (M1 Max etc.), NTT working-set L2 pressure dominates over compute parallelism.
- The thread_local caches above remain per-pool-worker. Pool warm-up at startup is recommended (issue one large mult per worker).
- The user-facing thread safety guarantees above do not change. Internal parallelism is an implementation detail of single operation calls, not a change in the concurrency model.

See [`DIVISION.md` §Improving skewed division beyond the current floor](DIVISION.md#improving-skewed-division-beyond-the-current-floor) for the design rationale and expected wins.
