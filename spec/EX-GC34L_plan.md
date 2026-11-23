# EX-GC34L — is J1's first-order gravitomagnetism error in the source, the curvature readout, or the geometry? (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "yes, go ahead. write the plan and run it." This followed the
ChatGPT review of 34K (`review/EX-GC34K_chatgpt_review.md`) and my reply proposing this study.

This plan and the `explore/LOG.md` entry are pushed before any 34L computation. Cap: 2 core-hours on the cloud; I
expect under 0.3. Our construction, not LaViolette's Model G.

## Held fixed

- The vacuum action: J1 (34J, forward gradient, range 1, μ = −1), from `results/EX-GC34J/kernels.npz`.
- The 34D static kernel and the GR comparison (`exgc34e.gr_solution`, `exgc34e.riemann`).
- The 9 cases of 34K: q along (1,0,0), (1,1,1) and (0.8,0.5,0.33), each with v parallel, perpendicular and oblique.
- The spacings kh = 0.0125, 0.025, 0.05 and 0.1, and |v| = 1e-3. The O(v) part is taken from ±v, as before.

## Pieces

**Gauge directions (11 components per vertex).**
- Displacement: G(q, ω) (3 columns), as in `exgc34e.gauge_modes`.
- Re-slicing: z = (n = −iω, β = μ g(q), edges 0), where g is the forward gradient, as in 34I/34J `I0`.
- Together they form W = [G, z] (11 × 4).

**Sources:**
- **S0 (current):** f = κ (1, −v, 0): the momentum sits at the vertex. Its re-slicing overlap is about kh/2 of the
  time-derivative term.
- **S1 (derived from discrete continuity):**
  - The momentum is a current on the links: f_β,a = −κ v_a e^{iq_a/2} (q_a/2)/sin(q_a/2).
  - This is the per-axis current that satisfies the lattice continuity equation −iωρ + Σ_a ∇⁻_a j_a = 0 exactly,
    with ω = k·v and ∇⁻ the backward difference (the adjoint of g).
  - So z†f = 0 to rounding. It tends to S0 as q → 0.

**Readouts (lattice x → continuum metric h, then linearised Riemann with k_μ = (−ω, k)):**
- **R0 (current):** h00 = −2n, h0i = β_i at the vertex, and h_ij from the edges (34D map).
- **R1 (gauge-corrected, primary):** L = L0 + ΔL, where:
  - L0 is R0 written as a 16 × 11 matrix;
  - ΔL = (Gc A − L0 W) W⁺, with W⁺ = (W†W)⁻¹W†;
  - Gc(q, ω) is the continuum gauge generator h_μν += i(k_μ ξ_ν + k_ν ξ_μ);
  - A = Gc⁺ L0 W is the least-squares continuum gauge that matches each lattice gauge direction.

  With this choice L W = Gc A, so R1's Riemann is exactly unchanged by either lattice gauge.
- **R1w (robustness, reported):** R1 with W⁺ taken in a weighted inner product on x, with block weights (lapse,
  shift, edges) = (2, 1, 1) and (1, 1, 2).
- **R2 (reported):** R0 with β read half a cell back, β_a e^{−iq_a/2}. This is the 34K A3 readout −½.

**The solve.** As before, M x = f with M = K_J1(q, ω) + static, and ω = k·v. Here the penalty covers both gauges,
M + W W†, so the re-slicing content of x is fixed by W†x = 0 and no longer left to the minimum-norm choice of a
least-squares solver.

## Part A — audit of the current source and readout (review step 1, reported)

For S0 and S1, in all 9 cases at all 4 spacings:
- the overlap with each gauge direction, |w†f|, divided by |ω κ| for z and by |f| for G;
- the forced residual |M x − f|/|f|;
- the per-component O(v) curvature error, as in 34K.

**Readout gauge sensitivity.** For R0, R1 and R2, at each spacing for the worst case of 34K:
- add each of the 4 gauge directions to the solution, scaled to the norm of x;
- report the largest change in the O(v) Riemann tensor relative to max |R_GR,O(v)|.

## Part B — the corrected benchmark (review steps 2 and 3, registered)

**Primary combination:** S1 + R1. E4 is the maximum over the 9 cases of the O(v) Riemann error relative to max |GR|,
as in 34E/34K, plus the static error.

**Registered outcome for (S1, R1):**
- **L-PASS** if E4 ≤ 1e-2 at kh = 0.05 and the order ratio E4(0.05)/E4(0.025) ≥ 3. Then the first-order error came
  from the source and the readout, and J1 stands as the vacuum action.
- **L-FAIL** otherwise.

Also reported: the full 4 × 4 table of E4 and order ratios for S0/S1 × R0/R1/R2/R1w. This separates the source
effect from the readout effect.

## Part C — only if L-FAIL: can the leading coefficient be removed within the exact kernels? (review step 3, registered)

- **The kernels.** Exact kernels at range 1 and range 2: J1 + N y, where N is the null space of the identity
  (dimensions 123 and 1263).
- **Fitting the error.** With (S1, R1), I fit the leading O(kh) error vector over all 9 cases and 36 components
  (c = err/kh at kh = 0.0125) as a linear function of y:
  - finite differences on each null-space direction, using a basis of N and step 1e-3 × ‖J1 − T‖;
  - at range 2, a random 300-dimensional subspace (seed fixed), if the full basis exceeds 15 CPU-minutes.
- **Minimising.** I find the y that minimises max |c|, as a linear programme with ‖N y‖ ≤ 0.5 ‖J1 − T‖.
- **Checking.** For the resulting kernel I recompute E4 directly at all 4 spacings and check G1 on the reduced set
  (axes, signed diagonals, generic; kh = 0.05, 0.2, 1, 3).
- **Verdict:**
  - **C-REMOVABLE** if that kernel has E4 ≤ 1e-2 at kh = 0.05, an order ratio ≥ 3, and G1 holds.
  - **C-FIXED** if the minimised linear max |c| stays ≥ 0.1, i.e. the first-order coefficient cannot be brought below
    0.1 within the bound.
  - **C-PARTIAL** otherwise.

## My expectations, stated before running

- **Audit:** S1 has re-slicing overlap ≤ 1e-12. R0 and R2 change at first order under re-slicing (R2 also under
  displacement). R1 does not change under either gauge.
- **Part B:** I lean towards L-FAIL. In 34K the closest earlier analogue, source +½ with readout −½, stayed first
  order. But R1 is new, so this is a weak lean.
- **Part C:** no expectation.

## Scope

Linear, flat background, continuous time. A gauge-invariant readout is a choice. R1 is the minimal correction in the
Euclidean inner product, and R1w shows how much the result depends on that choice.

## Files

`src/exgc34l.py` writes to `results/EX-GC34L/`:
- `audit` → `audit.json` (A);
- `bench` → `bench.json` (B);
- `fit` → `fit.json` and `kernel_C.npz` (C, only if L-FAIL).

It reuses `src/exgc34e.py`, `exgc34i.py` and `exgc34k.py` unchanged.
