# EX-GC34L — J1's first-order gravitomagnetism error is not in the source or the readout; it is in the geometry (EXPLORATORY)

**PI CLOSED 2026-11-23** (Richard: "Close 34K and 34L", after ChatGPT's follow-up review `review/EX-GC34J_followup_chatgpt_review.md`).

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "yes, go ahead. write the plan and run it." This followed the
ChatGPT review of 34K (`review/EX-GC34K_chatgpt_review.md`).
- The plan `spec/EX-GC34L_plan.md` and the `explore/LOG.md` entry were pushed before any run (c9e8537).
- Code: `src/exgc34l.py`.
- Outputs: `audit` → `audit.json`; `bench` → `bench.json`; `fit` → `fit.json`, `kernel_C.npz`.

This is our construction, not LaViolette's Model G.

## Conclusion

**L-FAIL, then C-FIXED.** Correcting the source and making the readout gauge-invariant does not remove J1's
first-order error. Nor can any exact kernel at range 1 or 2 remove it.

**1. The source mismatch the review found is real, and the derived source removes it exactly.**
- The current source (S0) overlaps the re-slicing direction by kh/2 of the time-derivative term: 2.5% at kh = 0.05
  in the worst case. Its forced residual is up to 5e-5.
- The source derived from exact lattice continuity (S1) puts the momentum on links as a current, per axis
  j_a = ρ v_a e^{iq_a/2} (q_a/2)/sin(q_a/2). Its overlap is ≤ 2e-15 in all 36 case–spacing combinations.
- But E4 stays first order. With the invariant readout, E4 is 1.55e-2 at kh = 0.05, against 1.25e-2 with S0, and
  the ratio is 2.00. The worst case moves from the axis-parallel source to the generic oblique one.

**2. The current readout is gauge-dependent, but the first-order error is not a readout artifact.**
- **The current readout (R0) depends on the gauge.** Adding re-slicing content equal to the solution's norm changes
  its O(v) Riemann by 2.7% of the signal at kh = 0.05. That change is itself first order: 0.67%, 1.3%, 2.7% and 5.4%
  at the four spacings.
- **The corrected readout (R1) is exactly invariant.** It changes by ≤ 2e-13 under any of the 4 gauge directions.
- **R1 gives the same E4 as R0.** The solve now fixes both gauges (W†x = 0), and on such solutions R1 = R0. The
  numbers match to 5 digits.
- **The 34K solve gives the same E4 too.** As a separate check I used 34K's solve (re-slicing left free) and added
  an arbitrary re-slicing. R1 then gives 1.55211407e-2 at kh = 0.05, the same value.
- **The readout's inner product moves the coefficient only slightly.** The weighted variants R1w give 1.49–1.55e-2,
  a 4% change, still first order.

**3. The leading coefficient is essentially fixed across all exact kernels (Part C, registered).**
- I linearised the 324-component first-order coefficient vector (9 cases × 36 components, real and imaginary parts)
  in the full null space: 123 dimensions at range 1 and all 1263 at range 2.
- Minimising its maximum within ‖N y‖ ≤ 0.5‖J1 − T‖ only lowers it from 0.310 to 0.292 at range 1, and to 0.292 at
  range 2. The radius bound is active at both.
- Both optimised kernels develop a growing mode at kh = 3:
  - range 1: Im ω/k = 5.0 and 11.2;
  - range 2: Im ω/k = 1.9.
- Checked directly, both are still first order.
  - Range 1 has a ratio of 2.09.
  - Range 2's ratio of 3.1 is not an improvement: its error grows at every spacing, for example 3.2e-2 at
    kh = 0.05 against J1's 1.55e-2.
- **Verdict: C-FIXED** (linear min-max ≥ 0.1).

**What this means.** The first-order term belongs to J1's lattice geometry: the forward re-slicing gradient together
with the vertex lapse coupling of the 34D potential. It is not the moving source, the curvature readout, or the choice
among exact kinetic kernels.
- The 34K post hoc result showed that a centred gradient cannot satisfy the identity with this potential.
- So second order would need a change on the potential side, such as where the lapse couples or how it is
  discretised. That would make a different vacuum action, not J1.
- Option (a), a finer lattice, remains the working answer: E4 ≤ 1e-2 needs kh ≤ 0.03 with S0, or about kh ≤ 0.032
  with S1 (interpolated, kh/3.1).

**A correction to 34K.** 34K A2 found the worst component's coefficient spread by 91% across exact kernels. That used
the inconsistent source, with the re-slicing content left to the solver. Here, with a consistent source and an
invariant readout, the maximum coefficient can move by only about 6% within the same bound. So most of 34K's
spread was probably readout and gauge dependence, not real freedom in the kernel. This is an inference from the
linearised fit. I did not rerun 34K's random samples with S1/R1.

## Registered results

**Part B (S0/S1 × readouts; E4 = maximum over 9 cases of the O(v) Riemann error relative to max |GR|):**

| source × readout | kh 0.0125 | 0.025 | 0.05 | 0.1 | ratio | worst case |
|---|---|---|---|---|---|---|
| S0 × R0 (34K benchmark) | 3.13e-3 | 6.25e-3 | 1.251e-2 | 2.51e-2 | 2.00 | (1,0,0) parallel |
| S0 × R1, R1w | same | same | same | same | 2.00 | same |
| S0 × R2 | 3.79e-3 | 7.59e-3 | 1.52e-2 | 3.03e-2 | 2.00 | generic oblique |
| S1 × R0 | 3.88e-3 | 7.76e-3 | 1.552e-2 | 3.11e-2 | 2.00 | generic oblique |
| **S1 × R1 (registered)** | **3.88e-3** | **7.76e-3** | **1.552e-2** | **3.11e-2** | **2.00** | generic oblique |
| S1 × R1w (2,1,1) | 3.73e-3 | 7.45e-3 | 1.49e-2 | 2.98e-2 | 2.00 | generic oblique |
| S1 × R1w (1,1,2) | 3.88e-3 | 7.76e-3 | 1.552e-2 | 3.11e-2 | 2.00 | generic oblique |
| S1 × R2 | 3.13e-3 | 6.25e-3 | 1.251e-2 | 2.51e-2 | 2.00 | (1,0,0) parallel |

**Verdict: L-FAIL** (1.55e-2 > 1e-2; ratio 2.00 < 3). S1 × R2 is S0 × R0 shifted by half a cell at both ends, which
is why they agree.

**Part A audit at kh = 0.05** (`audit.json` has all spacings):

| case | S0 re-slicing overlap / ωκ | S1 overlap | S0 forced residual | S1 forced residual |
|---|---|---|---|---|
| (1,0,0) parallel | 2.5e-2 | 7e-16 | 2.5e-5 | 4.7e-6 |
| (1,0,0) oblique | 1.5e-2 | 4e-16 | 1.6e-5 | 4.4e-6 |
| (1,1,1) parallel | 1.4e-2 | 2e-15 | 1.5e-5 | 2.0e-6 |
| generic parallel | 1.7e-2 | 1e-15 | 1.7e-5 | 2.4e-6 |
| perpendicular cases | 0 to 3.2e-3 | ≤ 1e-16 | ≤ 3.2e-6 | ≤ 2e-12 |

- The overlap grows linearly with kh: 6.3e-3, 1.25e-2, 2.5e-2 and 5.0e-2 at the worst case.
- The S1 residual of about 5e-6 is the displacement overlap, which is O(v²) (2.5e-7 of |f|).

**Readout gauge sensitivity** (change in the O(v) Riemann tensor, as a fraction of the signal, when gauge content
equal to the solution's norm is added; displacement x, y, z, then re-slicing):

| readout | kh 0.0125 | 0.05 | 0.1 |
|---|---|---|---|
| R0 | 2e-8, 2e-5, 2e-5, **6.7e-3** | 4e-7, 3e-4, 3e-4, **2.7e-2** | 2e-6, 1e-3, 1e-3, **5.4e-2** |
| R1, R1w | ≤ 2e-13 | ≤ 1e-15 | ≤ 2e-13 |
| R2 | 2e-5, 2e-5, 2e-5, 7e-6 | 9e-5, 3e-4, 3e-4, 1e-4 | 2e-4, 1e-3, 1e-3, 5e-4 |

**Part C** (with S1 × R1; c = first-order coefficient vector at kh = 0.0125):

| | range 1 | range 2 |
|---|---|---|
| null-space dimension used | 123 (full) | 1263 (full) |
| max \|c\| for the base kernel | 0.310 | 0.310 |
| linear min-max at ‖Ny‖ = 0.5‖J1 − T‖ | 0.292 | 0.292 |
| E4 at kh 0.0125 / 0.025 / 0.05 / 0.1 | 3.7e-3 / 7.4e-3 / 1.55e-2 / 3.6e-2 | 4.0e-3 / 1.0e-2 / 3.2e-2 / 1.08e-1 |
| ratio | 2.09 | 3.11 (error grows everywhere) |
| G1 reduced set | **FAIL** (2 growing at kh = 3) | **FAIL** (1 growing at kh = 3) |
| verdict | **C-FIXED** | **C-FIXED** |

**My stated expectations:**
- Audit: S1 overlap ≤ 1e-12, right. R0 changes under re-slicing, right. R1 is invariant, right.
- R2 also changes under displacement: right, but so does R0, at a similar level of about 3e-4 of the signal at
  kh = 0.05.
- Part B: I leaned towards L-FAIL. That was right.
- Part C: I had no expectation.

## Recomputation and shared dependencies

- **S0 × R0 reproduces the 34K benchmark** with a different solve (both gauges penalised): 1.2510e-2 at kh = 0.05,
  equal to 34J/34K's 1.25097e-2 to 5 digits.
- **The registered S1 × R1 number was recomputed by a separate path.** That path used 34K's solve (displacement
  penalty only, minimum-norm re-slicing), added an arbitrary complex re-slicing component, and applied R1. It gave
  1.55211407e-2 at kh = 0.05 and 7.7599e-3 at kh = 0.025, identical to `bench.json`.
- **Shared dependencies:** `exgc34e.riemann`, `gr_solution`, `metric_from_x`, `static_kernel`, and J1's coefficients.

## Compute

About 2 CPU-minutes in total:
- audit: 2 s;
- bench: 3 s;
- fit: 97 s.

Platform: Linux cloud, Python 3.13, numpy 2.5, scipy 1.18.

## Deviations

- **The audit's re-slicing overlap is normalised by κ|k||v| instead of |ωκ|.** That equals |ωκ| for v ∥ q, and stays
  defined in the perpendicular cases where ω = 0. The first audit run divided by ω ≈ 0 there, giving meaningless
  numbers. I changed the normalisation and reran before writing anything up. No verdict depends on it.
- **Part C at range 2 used the full 1263-dimensional basis,** as the plan allows below 15 CPU-minutes.
- None otherwise.
