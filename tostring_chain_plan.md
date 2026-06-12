# ToString Chain Optimization Plan

Status: PLANNED (written 2026-06-12, post PR #107–#112 MFA run).
This is the largest user-visible gap left vs GMP: ToString at 100k–2M
digits runs 1.85–4.0× GMP's time; at ≥5M digits the warm steady state is
already 1.17–1.28× (near parity) and the residual is cold chain build.

## Current state (post-#112 canonical run + A/B verification)

| digits | BigMath | GMP | ratio | note |
|---|---:|---:|---|---|
| 100k | 9.7 ms | 2.4 | 4.0× | worst point |
| 500k | 53.4 | 20.8 | 2.57× | |
| 1M | 105.6 | 49.0 | 2.16× | |
| 2M | 220.4 | 119.1 | 1.85× | |
| 5M | 510 warm / 879 cold | 398.6 | 1.28× / 2.2× | cold = first call at a size |
| 10M | 1 062 warm / 1 699 cold | 908.3 | 1.17× / 1.87× | |

Architecture (docs/STRING_CONVERSION.md): D&C formatter over a
`DecimalDcChain` — entries at `d = L/2, L/4, …` down to ~512 digits, each
holding `Pow10(d)` and a precomputed `NewtonDivision::Divider`
reciprocal. Chain is cached `thread_local` keyed by `topDigits`
(src/common/Parser.cpp:339), so repeat calls at the same size are warm;
the cold cost is one `Pow10` + one full Newton reciprocal build per
level. Per-level divmods then cost O(M(n)).

## Where the time goes (to be confirmed by profile — do this FIRST)

Hypothesized split, cold call at 1M digits:
1. Chain build: log₂(L/512) ≈ 11 Divider setups, each a full Newton
   reciprocal at its level's size. Dominant cold cost.
2. Recursive divmods: 2^k divmods at level k, each O(M(L/2^k)) — the
   steady-state cost, entirely SERIAL today (zero ParallelDo/ParallelFor
   in Parser.cpp).
3. Leaf linear formatting (Möller-Granlund div2by1 loop): already tight.

Profile protocol: `sample` on a cold 1M-digit ToString and a warm one;
attribute to BuildDecimalDcChain vs ToStringDivConquer vs leaf. The
levers below are ordered by expected (win ÷ risk) — re-order after the
profile.

## Levers

### T1. Parallelize the D&C recursion (high win at ≥500k digits, low risk)

After each divmod, the (q, r) subtrees are fully independent and their
output substrings land at disjoint, precomputed offsets (`padTo`
arithmetic already determines exact digit positions). Today the
recursion is serial on the caller's thread while the pool sits idle.

Constraint: the pool is non-reentrant (single work slot — see
docs/MULTIPLICATION.md threading notes), so nested ParallelDo from
inside a worker deadlocks. Shape the parallelism as ONE dispatch:
descend serially to depth k (2^k subtrees, k ≈ 3–4 → 8–16 subtrees of
similar size), collect them in a work list, `ParallelDo(2^k)` over the
list with each worker running its subtree serially into its own output
range, then the already-written shared `out` buffer needs no stitch.
The top-of-tree divmods before depth k stay serial — they are the big
O(M(L)) ones and can themselves use threaded multiplies internally
(they already do: Newton divmod → NTT → ParallelDo inside), so the
machine is busy throughout.

Expected: divmod tree below depth k is ~half the total divmod work;
8-way parallel ⇒ 1.3–1.7× on the warm path at ≥500k digits. This is
also the only lever that helps the WARM 5M–20M rows (1.17–1.28× → at or
below GMP).

Risk: low — output ranges are disjoint by construction; the workers
call Divider::DivideAndRemainderInto which must be checked for hidden
shared state (the chain entries are read-only after build; verify the
Divider is const-callable / has no internal scratch races — if it has
thread_local scratch, that is per-worker and safe).

### T2. Bottom-up reciprocal chain (kills most of the cold cost, medium risk)

Chain entries are exact halvings and `Pow10` is already built by
squaring (`Pow10(d) = Pow10(d/2)²`). Reciprocals compose the same way:
`1/10^d = (1/10^(d/2))² · 10^... ` — i.e. the level-d reciprocal is the
square of the level-d/2 reciprocal up to scaling, correct to the
smaller precision; ONE Newton refinement step lifts it to full
precision at level d. So build dividers bottom-up: full Newton build
only at the smallest level (~512 digits, trivial), then each larger
level costs one squaring + one refinement multiply instead of a full
reciprocal build. Halves-to-thirds the chain-build cost; the cold 100k
(4.0×) and 500k (2.57×) rows are mostly this.

Requires a `NewtonDivision::Divider` constructor that accepts a seed
reciprocal at half precision — check whether the invertappr-style
internal loop already supports a warm start (it iterates from a seed
anyway; expose the entry point).

### T3. Persist the chain across sizes (small win, trivial)

The cache is keyed by exact `topDigits`; a 1M-digit call then a
999 873-digit call builds two chains, though all lower levels could be
shared. Key the cache by level value (`d`) instead of by top size:
each `(d, value, divider)` triple is independent of the top. Saves
rebuilds in mixed-size workloads (calculator REPL, BigDecimal toString).
No effect on single-size benches — do not expect bench movement; this
is API-hygiene-level.

### T4. Leaf threshold re-sweep (cheap, do alongside T1/T2)

`BIGMATH_TOSTR_DC_THRESHOLD = 1024` and the chain floor (~512 digits)
were tuned before the MFA run cheapened every multiply ≥2^20
coefficients and before cyclic products reached the fused pipeline. The
leaf/linear crossover and `DecimalDcThreshold = 2048` (parse) deserve a
re-sweep — same playbook as the division-floor re-sweeps (PRs #92/#101,
and the #109 gate retune that this run showed pays compounding
dividends).

## Measurement methodology

- Warm/cold split is mandatory (PR #112 lesson: single-iter chain rows
  vary ±20–60% with process state). Per size: 1 cold call (fresh
  process), then best-of-3 warm — report both.
- Sizes: 100k / 500k / 1M / 2M / 5M / 10M / 20M digits, interleaved
  baseline-vs-change, quiet machine.
- GMP reference: `mpz_get_str` timed in-process (bench_vs_gmp rows).
- Verify output identity: byte-compare BigMath vs GMP strings at every
  measured size, every change (cheap, total).

## Correctness protocol

ToString output is end-to-end checkable — byte equality vs
`mpz_get_str` IS the oracle. Every lever change must:
1. Byte-compare vs GMP at 10k / 100k / 1M / 10M digits, random inputs,
   plus all-9s and power-of-10 boundary values (carry/pad edge cases —
   the `padTo` zero-padding logic is where T1's offset math can break).
2. Keep the existing parse/ToString round-trip tests green.
3. T1 gets a transient-break: shuffle the work-list order and confirm
   output corrupts (proves offsets, not luck, give the ordering); then
   restore.

## Sequencing

1. Profile (half a day) — confirms the T1/T2 split before any code.
2. T1 (parallel subtrees) — single PR, biggest steady-state lever.
3. T2 (bottom-up reciprocals) — single PR, biggest cold lever.
4. T4 threshold re-sweep — fold into whichever of T1/T2 lands last.
5. T3 only if a real mixed-size workload shows up.

Target: 100k digits 4.0× → ≤2×; 1M 2.16× → ≤1.3×; warm 5M–20M from
1.17–1.28× → parity or better.
