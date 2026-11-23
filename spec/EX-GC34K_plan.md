# EX-GC34K — where the first-order gravitomagnetism error of the re-slicing-exact kinetic term comes from, and whether another gradient removes it (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "You can close them and use any inputs you get from the review
for the next test decision. I don't have any bias toward the 3 options. Maybe you should do all 3."

This plan and the `explore/LOG.md` entry are pushed before any 34K computation. Cap: 2 core-hours on the cloud; I
expect under 1. Our construction, not LaViolette's Model G.

## Inputs

**34J (PI closed).** J1 is the forward-gradient, range-1, re-slicing-exact kinetic term with μ = −1. It is stable and
has two isotropic waves at c, but its E4 gravitomagnetism error is first order: about kh/4, or 1.25e-2 at kh = 0.05.

**ChatGPT review (`review/EX-GC34I-34J_chatgpt_review.md`):**
- (i) Keep J1 as the linear vacuum benchmark.
- (ii) Before searching wider stencils, find the case and curvature component that set the error. Then test whether
  its leading coefficient varies within the identity-obeying kernel space, with source and readout consistent.
- (iii) If it is fixed, test a staggered shift, and seed the time domain with constraint-compatible data from the
  outset.

**The three options from 34J, all done here as the PI suggested:**
- (a) accept first order with a finer lattice;
- (b) a staggered shift;
- (c) a centred gradient at a wider range.

**A note on option (b), from the algebra before any run.** Suppose β is stored on half-cells, β_a at v + e_a/2, and
interpolated to vertices by averaging. Then the displacement gauge stays exact if u is staggered the same way. The
re-slicing gradient becomes 2i sin(q_a/2)/h on the half-cell, and its vertex value cos(q_a/2) · 2i sin(q_a/2)/h is the
central difference i sin(q_a)/h. So, in the re-slicing identity, the averaged staggered shift is the same as option
(c). Parts B and C below therefore test one family of re-slicing gradients that covers both (b) and (c).

## Part A — diagnosis (review priority ii)

- **A1: which case and component.** For J1, at kh = 0.0125, 0.025, 0.05 and 0.1, I record the O(v) Riemann error
  per case (3 q directions × 3 velocity directions) and per component (36). I report:
  - the case and component that set the maximum;
  - its leading coefficient c₁ = err/kh, the limit from the two smallest kh.
- **A2: does c₁ vary within the identity-obeying space?** Exact solutions are J1 + N z, where N is the null space of
  the identity at that range.
  - Range 1: 20 random z, with ‖N z‖ = 0.1 and 0.3 × ‖J1 − T‖ (10 each).
  - Range 2: 10 random z, from J3's space, at 0.1 × ‖J3 − T‖.
  - For each kernel I compute c₁ of the worst case and component from A1, and the maximum over all cases.
  - **A-fixed** if the relative spread of c₁ (max − min over the mean, over all samples, J1, J2 and J3) is < 1%.
  - **A-variable** otherwise. Then I also report the smallest c₁ found, and whether that kernel keeps G1 on a
    reduced direction set (axes, signed diagonals, generic; kh = 0.05, 0.2, 1, 3).
- **A3: consistent placement (reported).** The 9 combinations of source and shift-readout half-cell shifts (0, ±½)
  for J1, as in the 34J pre-registration probe, with c₁ for each.

## Part B — existence of local solutions for a family of re-slicing gradients (options b and c)

- **The family.** One-parameter, per axis: g_γ,a = ((1 + γ) e^{iq_a} − (1 + 2γ) + γ e^{−iq_a})/h.
  - γ = 0 is forward, γ = −½ central, γ = −1 backward.
  - Tested: γ ∈ {0, −0.25, −0.5, −0.75, −1}.
- **Also:** the "Kuhn-edge" gradient, the least-squares gradient fitted to the 7 forward edge differences of ξ.
- **Ranges.** R = 1 and 2 for every member. R = 3 for γ = −0.5 and γ = −0.25 if they have no solution at R ≤ 2.
- **Method.** As in 34I and 34J:
  - real Hermitian coefficients;
  - exact on an FFT grid (half the grid, using q → −q symmetry);
  - solution exists if the consistency residual is ≤ 1e-10;
  - μ = −1;
  - closest to T.
- **Reported for every member and range:** existence, the residual, the solution-space dimension, and c₁ (E4 at
  kh = 0.025 and 0.05) for members that exist.

## Part C — full scoring

**Scored candidates:**
- **J1** (registered again, for the time domain and for option a);
- **every Part B member** with a solution and E4 order ratio E4(0.05)/E4(0.025) ≥ 3.

**Criteria:**
- **I0**, as in 34I.
- **G1, G2, G3 and G7**, as in 34G/34I, with the full direction and spacing set.
- **G5**, as in 34J, including E4 ≤ 1e-2 at kh = 0.05 and the order ratio ≥ 3.
- **G4″: synchronous-gauge time domain**, generic and (1,1,0) at kh = 0.05, 40 periods, steps 0.02/ω_max and half
  that. Seeds are constraint-compatible from the outset:
  - every non-zero eigenmode that satisfies the lapse and momentum constraints and their time derivatives (relative
    residual ≤ 1e-8);
  - plus the constraint-satisfying part of the zero cluster (|ω| ≤ 0.01 k);
  - each seeded at 1e-8.
- **G4″ passes if:**
  - constraint residuals ≤ 1e-9;
  - the wave's tidal amplitude changes by ≤ 1e-6;
  - each seed group's tidal change is ≤ 1e-6 of the wave's;
  - both steps agree.
- **Option (a), reported for J1:** E4 at kh = 0.03 and 0.04, and the largest tested kh with E4 ≤ 1e-2.

**A candidate passes** if it meets I0, G1–G5 (including E4 and its order) and G7. J1 is expected to fail E4. Its
G4″ is scored, so the time-domain verdict for J1 becomes a registered one.

**My expectations, stated before running:**
- A: I lean towards A-fixed, but the review is right that it is unproven.
- B: forward (γ = 0) and backward (γ = −1) have solutions. Central has none at R ≤ 2, and I guess none at R = 3
  either. For intermediate γ I don't know. If only one-sided members exist, they are all first order.
- C: J1 passes G4″ but fails E4.

## Scope

Linear, flat background, continuous time.

## Files

`src/exgc34k.py`:
- `diag` → `results/EX-GC34K/diag.json` (A1–A3);
- `family` → `family.json` + `kernels.npz` (B);
- `spec` → `candidates.json`;
- `sync` → `sync.json`;
- `e4` → `e4.json`;
- `check` → `check.json` (C).

It reuses `src/exgc34e.py`, `exgc34f.py`, `exgc34g.py`, `exgc34i.py` and `exgc34j.py` unchanged.
