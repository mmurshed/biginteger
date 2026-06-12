# Small-Skew Division Plan (10k–2M-digit dividends)

Status: DONE (S1 executed 2026-06-12; see outcome at bottom).
Honest expectations up front: this band is adjacent to documented dead
ends (sub-50k mul vs GMP's hand-tuned basecase), and GMP's lead here is
basecase quality, not algorithm choice. The realistic goal is to halve
the gap, not flip it. Read the cost/benefit note at the bottom before
starting.

## Current state (canonical run, PR #112)

| shape (digits) | limbs (÷19.27) | BigMath | GMP | ratio |
|---|---|---:|---:|---|
| 40k ÷ 10k | 2.1k ÷ 520 | 0.74 ms | 0.22 | 3.35× |
| 100k ÷ 10k | 5.2k ÷ 520 | 2.16 | 0.45 | 4.84× ← worst |
| 200k ÷ 50k | 10.4k ÷ 2.6k | 3.83 | 1.69 | 2.26× |
| 500k ÷ 100k | 26k ÷ 5.2k | 7.65 | 4.63 | 1.65× |
| 1M ÷ 200k | 52k ÷ 10.4k | 14.6 | 10.0 | 1.46× |
| 2M ÷ 500k | 104k ÷ 26k | 28.4 | 24.6 | 1.15× |
| 5M ÷ 1M | 260k ÷ 52k | 70.5 | 69.8 | 1.01× (parity from here up) |

Divisors in the losing rows are 520–26k limbs. Dispatch context
(DispatchThresholds.h): BZ from divisor 512; Newton bands from divisor
1024 (5/2), 4096 (8/5), 768 (8/1 high-skew); below those, FastDivision
(Knuth D with Möller-Granlund 3/2). The 520-limb divisor of the two
worst rows sits at the BZ floor and below every Newton band.

## Why the gap exists (hypotheses — profile FIRST, same rule as always)

1. **Stale floors.** Every Newton/BZ floor was tuned 2026-06-11/12
   BEFORE the MFA run made Newton's internal multiplies and cyclic
   products 1.2–3× cheaper (PRs #107–#111). The same staleness that hid
   the 2^24→2^20 MFA gate win (PR #109 — floors tuned against a slower
   engine) almost certainly affects `BIGMATH_NEWTON_MEDIUM_B = 1024`,
   `BIGMATH_NEWTON_RATIO2_B = 4096`, `BZ_DIVISOR_THRESHOLD = 512`, and
   `BIGMATH_CYCLIC_NTT_THRESHOLD = 1280`. A floor re-sweep is the
   cheapest possible experiment and has paid out twice already
   (PRs #92, #101).
2. **10:1 shape runs quotient-sized work in the wrong place.** 100k÷10k
   has a 4.7k-limb quotient against a 520-limb divisor; the
   QuotientSizedDivision band starts at b ≥ 8192 (thin-quotient) /
   24576. A small-divisor high-skew path that chunks the dividend
   against a cached divisor reciprocal (the ToString chain pattern —
   `NewtonDivision::Divider` reuse across chunks) may beat FastDivision
   well below the current floors now that reciprocal builds are
   cheaper.
3. **FastDivision basecase itself.** GMP's `mpn_sbpi1_div_qr` /
   divide-and-conquer ladder is hand-tuned assembly territory. If the
   profile says the time is in the 3/2 GM inner loop rather than in
   dispatch-band misses, stop — that is the documented-dead-end wall,
   and the remaining ratio is the price of portable C++.

## Levers, ordered

### S1. Post-MFA-run floor re-sweep (cheap, do first, likely the whole win)

Re-run `dispatch_tuner` / `division_shape_bench` across divisor sizes
256–8192 × ratios {2, 4, 8, 10, 16} against the current engine. Sweep
`NEWTON_MEDIUM_B`, `NEWTON_RATIO2_B`, `BZ_DIVISOR_THRESHOLD`,
`NEWTON_HIGH_SKEW` floor (768), and `CYCLIC_NTT_THRESHOLD` (1280)
downward. Warm-state, interleaved, quiet machine — the PR #107–#112
methodology. Expected from history: each prior re-sweep moved its band
1.5–2.5×; even half of that closes 200k÷50k to ~1.5× and 1M÷200k to
near parity.

### S2. Small-divisor high-skew Divider path (medium effort)

For `b < NEWTON floors` and `a ≥ 8b`: build one cached
`NewtonDivision::Divider` for the divisor (cost amortizes over a/b
chunks), then chunk the dividend top-down — the ToString-chain
consumption pattern, which the formatter already proves out at
520–5 200-limb divisor sizes. Prototype = wire
`Divider::DivideAndRemainderInto` into a loop in Division.h dispatch
under a new band; measure before any polish. If S1 already moved the
floors below 520 limbs, skip — S1 subsumes this.

### S3. Accept and document (the likely end state for ≤100k÷10k)

If the profile shows GM 3/2 inner-loop dominance after S1/S2, record
the residual as the portable-C++ basecase price next to the sub-50k
mul dead-end note, with the profile attached, and stop. Do NOT reach
for assembly or SIMD here without a separate cost/benefit discussion —
the band's absolute times are 0.7–15 ms.

## Measurement + correctness

- `division_shape_bench` + the raw-limb probe pattern: extend
  `mfa_mul_gmp_probe` CheckDiv with the band's shapes (520–26k-limb
  divisors) — quotient limb-compare vs `mpz_tdiv_q` is the oracle.
- `div_correctness` (identity q·b + r = a, r < b) stays the canonical
  harness; any new dispatch band must appear in its cross-check matrix.
- Floors changed by S1 get the knife-edge treatment from PR #92's
  lesson: test divisor sizes at band±1 limb and 2^k±1 around each new
  floor (the BZ non-pow2 pathology lives exactly there).

## Cost/benefit note (read before starting)

The whole band is single-digit milliseconds. ToString
(tostring_chain_plan.md) repays effort better: its worst row is 4× at
a user-facing operation, and its T1/T2 levers are structural rather
than tuning. Do S1 here (it is a half-day sweep with proven odds),
then decide between S2 and the ToString plan with the numbers in hand
— do not run both speculative tracks at once.

## Outcome (2026-06-12, executed)

S1 ran and found something better than stale floors: **BZ fell back to
FastDivision wholesale whenever the divisor size was odd** at entry or
at any recursion level (`BurnikelZieglerDivision.h`), i.e. on ~half of
all real divisor sizes BZ was silently O(n·Δ). This was the root cause
of the documented "2^k+1 family" pathology, and it had biased every
prior Newton/BZ floor sweep (smooth probe sizes measured real BZ; real
shapes often got Fast).

Fix: pad both operands with k bottom zero limbs so the divisor size is
`ceil(n/2^d)·2^d` (d = recursion depth); quotient unchanged, remainder
drops the k zero pad limbs. Overhead ≤ ~1/BZ_THRESHOLD. Measured (odd
sizes, dispatch before → after): 8193÷ at ratio 1.5: 43.5 → 4.9 ms;
12289 at 1.5: 95 → 8.1; 4097 at 1.5: 10.5 → 2.2; 2049 at 2.0:
5.3 → 1.8; 1037 at 2.0: 1.35 → 0.67.

Newton frontier re-swept against the fixed BZ (equal-rep interleaved
probe, `tests/performance/division_floor_probe.cpp`): now
(640, 8/1), (1024, 7/2), (1280, 14/5), (1792, 5/2), (4096, 8/5),
(24576, 4/3). Knife-edges at integer digit-ratios stay on the Newton
side (7/2 covers 4.0000±1, 14/5 covers 3.0000±1).

The plan's two worst rows (520-limb divisor) were already optimally
dispatched — BZ wins them, all our algorithms are ≥2.5× GMP there, and
the gap is basecase quality (hypothesis 3 — S3 accepted, documented
here). The MFA-staleness premise (hypothesis 1) was wrong for this
band: the MFA gate is 2^20 transform length, far above these sizes.
S2 is moot: blockwise Newton with a cached reciprocal IS the chunking
pattern, and it loses below 640 limbs.

Follow-up (same day, second PR): with padding making recursion depth
size-independent, the BZ basecase `BIGMATH_BZ_RECURSION_THRESHOLD`
swept 512 → 128 — another 1.2-1.5× across the whole BZ band (flat
64-128, worse at 48), which moved every Newton floor up again and
raised the balanced floor 24576 → 131072 (padded BZ wins ratio 1.4-1.5
through ~98k limbs; QSized main floor decoupled at 24576). Even the
"basecase wall" rows improved: 5200÷520 limbs 2.5 → 1.7 ms. Frontier:
(896, 8/1), (1280, 7/2), (1792, 14/5), (2560, 5/2), (4096, 2/1),
(8192, 8/5), (131072, 4/3).
