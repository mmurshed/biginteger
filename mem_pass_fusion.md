# MFA Memory-Pass Fusion Plan

Status: PLANNED (written 2026-06-12, post PR #82–#105 optimization run).
Owner note: this is the last identified structural performance lever. Every
band below 50M digits is at or better than GMP parity; the 50M–200M-digit
band (balanced mul 1.28–1.33×, div 1.23–1.37× vs GMP) is
**memory-bandwidth-bound**, and cutting DRAM passes is the only mechanism
that moves it.

## Evidence this band is bandwidth-bound

Profiled 2026-05-31 (M1 Max, see memory + `docs/MULTIPLICATION.md`):

- IPC 2.49 on the NTT phases with thread scaling dead (+0.5% from extra
  threads) — cores are stalled on memory, not compute.
- 4-process aggregate throughput 2.6× a single process — the bandwidth is
  there, one transform's serial pass chain can't use it.
- The NEON pass (PRs #96–#99) confirmed it from the other side: a 4×
  compute-kernel win delivered −26…40% below the bandwidth wall and only
  −10…13% above it.
- A large CRT multiply streams ~400 MB through DRAM across its pass chain at
  100M-digit operands (6 forward + 3 inverse transforms, each MFA stage a
  separate plane sweep, plus pointwise, Garner/carry, pack/unpack).

## Current pass inventory (one CRT multiply, per prime ×3, MFA regime)

For each operand (×2), forward transform = `ForwardMFA`:
1. `FusedForwardA` — gather rows + row FFT(n2) + cross-twiddle (fused, PR #78/#72)
2. `FusedForwardB` — column FFT(n1) + scatter (fused)
   (non-fused fallback path: step2 sweep + `Transpose` + step5 sweep + `Transpose`)

Then:
3. Pointwise multiply — separate full-plane sweep over all 3 primes
   (`fa[i] = F::Mul(fa[i], fb[i])`)
4. `InverseMFA` — `FusedInverseA` + `FusedInverseB` (mirrors forward)
5. `FinalizeProduct` — Garner CRT recombination + carry propagation, a
   separate sweep reading all three primes' planes and writing limbs
6. `PackOperand` (input side) — one write sweep per prime per operand

Counting plane-sized DRAM traversals per multiply (N = transform length,
4 bytes/coeff, 3 primes): pack 6, forward 2×3×(2 stages), pointwise 3×2
(read both + write), inverse 3×(2 stages), finalize 3 reads + 1 write.
≈ 30+ plane passes. Each fusion below removes 3–6 of them.

## Fusion opportunities, ordered by (win ÷ risk)

### F1. Pointwise into the last forward stage of operand B (high win, low risk)

`ForwardMFA(fb)` finishes with `FusedForwardB` writing fb's plane; the
pointwise sweep then re-reads fa and fb and rewrites fa. Instead: after
operand A's forward is complete, run operand B's `FusedForwardB` tile loop
with a fused `fa[i] = Mul(fa[i], fb_tile[i])` at scatter time — fb's plane is
never written to DRAM at all, and fa is updated in the same pass.
Saves ~3 full-plane writes + 6 reads per multiply (per prime: fb write,
fb read, fa read absorbed). Constraint: A's forward must complete before B's
last stage (already true — transforms run per-prime batched; needs the
per-prime ordering A-then-B within each worker, which the batched ParallelDo
can express by pairing fa/fb of the same prime in one work unit).

### F2. Garner + carry into the first... last inverse stage (high win, medium risk)

`FusedInverseB` scatters each tile to the output plane; `FinalizeProduct`
then re-reads all three primes' planes serially. Instead: keep the three
primes' inverse tiles for the same coefficient range resident (the inverse
stage B tile loop already walks coefficient-contiguous ranges), and run
Garner + carry per tile range, writing limbs directly. The carry chain is
the hard part: carries propagate left-to-right across tile boundaries, so
tiles must be finalized in ascending coefficient order with a carried-in
value per boundary (single-threaded stitch, or compute per-tile with
carry-out and do one cheap fix-up ripple pass — the ripple touches limbs,
not coefficient planes, and limbs are 1/2 the coefficient bytes).
Constraint: requires the three primes' inverse stages to be coordinated
per-range instead of per-prime-batched — this is the biggest structural
change and collides with the current cross-prime threading model. Prototype
single-threaded first; the band is bandwidth-bound, so losing thread overlap
on this stage may cost nothing.

### F3. PackOperand into FusedForwardA's gather (medium win, low risk)

Stage A already gathers input rows into tiles; today it gathers from the
packed coefficient plane, which `PackOperand` produced in a separate sweep
from the limb array. Pack per-tile inside the gather instead (limb → 3-prime
coefficients on the fly): saves 6 plane writes + 6 reads per multiply.
Constraint: the three primes pack the same limbs — per-prime packing reads
the limb array 3× (limb bytes = coeff bytes ×2... still net positive since
the coefficient plane round-trip disappears). Alternatively pack all three
primes in one fused gather when the per-range coordination from F2 exists.

### F4. Transpose blocking / elimination audit (low win, already partial)

PR #78 fused the transposes into stages A/B for the FUSE path; the non-FUSE
fallback and the recursive (`n2 > LEAF`) path still do standalone
`Transpose` sweeps. Either widen the fused path's applicability (currently
gated to `n1, n2 ≤ LEAF`) or accept and document.

## Measurement methodology

- Primary: `mul_xl_bench` limb sweep 2M–10M limbs (38M–193M digits),
  best-of-3 interleaved baseline-vs-fused, **quiet machine only** (load < 2
  — this band's numbers are meaningless under ambient load).
- Secondary: `/tmp`-style paired digit-driver runs at 50M/100M/200M ×
  {balanced, 10:1 skew} + div 200M×40M (inherits via Newton).
- Bandwidth confirmation: `powermetrics`/Instruments DRAM counters before vs
  after — the claim is fewer bytes moved, not faster code; verify the bytes.
- Expected ceiling: pass inventory says F1+F2+F3 remove ~12–15 of ~30 plane
  passes → up to ~1.3–1.6× in the MFA band if fully bandwidth-bound.
  Target: balanced 50M–200M from 1.28–1.33× to ≈1.0× (GMP parity), 200M×40M
  div from 1.37× toward ~1.1×.

## Correctness protocol (non-negotiable, learned from PR #103)

The MFA band sits above every routine harness — an ordering bug here ships
silently. Every fusion step must:
1. Keep `unit_tests` MfaTransform green (round-trip + MFA-vs-plain at 2^14–2^17,
   they exercise the fused tile paths directly).
2. Keep `mfa_roundtrip` (n=2^23, ctest) green.
3. Run the raw-limb GMP probe (mpz_import/export — never ToString, it routes
   through the band being tested) at tiers 2^24/2^25 before merging.
4. Any new pass ordering gets a transient-break validation: re-introduce the
   wrong order, confirm the tests fail, restore.
5. Watch for the fallback signature: a Newton divide pinning one core for
   minutes = corrupt multiply layer, not a perf mystery.

## Sequencing

1. F1 (pointwise fusion) — single PR, measurable alone.
2. F3 (pack fusion) — single PR.
3. F2 (Garner fusion) — prototype branch first (single-threaded per-range),
   measure, then decide on the threading redesign.
4. F4 audit last, only if F1–F3 leave the band short of target.
