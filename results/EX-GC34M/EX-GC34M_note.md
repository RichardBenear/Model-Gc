# EX-GC34M — a scalar packet's own sources conserve energy exactly on J1's lattice, but not momentum locally (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "I suggest that the next step is to do number 1) Start with a
scalar wave packet before returning to the knot. Derive its energy, momentum, and stress sources by varying one
discrete matter action with respect to J1's lapse, shift, and spatial variables. The first gate is whether those
sources satisfy the lattice conservation identities as the packet moves."
- Source: `review/EX-GC34J_followup_chatgpt_review.md`, recommendation 1.
- The plan `spec/EX-GC34M_plan.md` and the `explore/LOG.md` entry were pushed before any run (f252989).
- Code: `src/exgc34m.py`.
- Registered outputs: `m0` → `m0.json`; `run` → `identities.json`.
- Post hoc output: `posthoc` → `posthoc.json`.

This is our construction, not LaViolette's Model G.

## Conclusion

**The packet does not satisfy J1's lattice conservation identities: M1 passes, M2 fails.**

- **The sources come from one action.** The scalar action is varied with respect to J1's lapse, shift and edge
  lengths. A direct finite-difference check agrees with the analytic sources to 6e-11 (M0).
- **Energy is conserved exactly (M1 PASS).** The re-slicing identity holds to rounding: ≤ 1.1e-14 relative, at every
  spacing, time and direction.
- **Momentum is conserved only globally (M2 FAIL).**
  - Total lattice momentum is conserved exactly.
  - Locally, the displacement identity is violated at **first order** with an O(1) coefficient: 14–20% of the
    momentum-change scale at h = 0.25, and 47–97% at h = 1. It halves when h halves (ratios 1.4–2.3).
  - So the identity that J1's displacement symmetry requires is not met. M x = f is inconsistent at first order for
    this packet.
- **The integrated sources are right.** Total energy, total momentum and the integrated stress tensor all converge to
  the continuum values at second order (M3; for example dE 7e-3 → 1.9e-3 → 4.7e-4). The stress is correct on average
  and only misplaced locally.

**The cause is post hoc, but clear-cut.** It is the same half-cell conflict that 34J guessed at and 34K/34L circled
around, now seen in the matter sector:
- J1's re-slicing symmetry uses a forward gradient. To make energy exactly conserved, the shift must couple along the
  forward axis edge, so the momentum density sits half a cell along each axis.
- J1's displacement symmetry puts the shift at the vertex (β += u̇ there).
- Moving the momentum density back to the vertex turns the trade around:
  - local momentum then converges at roughly second order (ratios 3.1–3.7; 0.6% / 1.8% at h = 0.25);
  - energy breaks at first order instead (16% / 14% at h = 0.25).

  **One shift location cannot serve both symmetries.** This is the review's step 3 case: "redesign the lapse and
  shift couplings together."

## Registered results

**M0, derivation check:**

| h | finite difference vs Σ s·x (Richardson) | lapse / shift / edges separately |
|---|---|---|
| 1 | 5.9e-11 | 2.7e-10 / 8.7e-11 / 4.0e-9 |
| 0.5 | 2.7e-12 | 4.3e-10 / 2.2e-10 / 3.0e-9 |

The parts use ε = 1e-5 central differences. M0 **PASS**.

**M1 and M2** (m = 1, |k₀| = 0.6, σ = 3, L = 32; residual relative to the scale of ∂_t s; t = 0, 6, 12):

| direction | h | r_E (energy) | r_P (momentum), t = 0 / 6 / 12 |
|---|---|---|---|
| axis | 1 | ≤ 2.5e-15 | 0.97 / 0.53 / 0.47 |
| axis | 0.5 | ≤ 3.4e-15 | 0.42 / 0.35 / 0.34 |
| axis | 0.25 | ≤ 9.4e-15 | 0.21 / 0.17 / 0.19 |
| generic | 1 | ≤ 2.7e-15 | 0.66 / 0.61 / 0.61 |
| generic | 0.5 | ≤ 3.4e-15 | 0.35 / 0.29 / 0.33 |
| generic | 0.25 | ≤ 1.1e-14 | 0.17 / 0.14 / 0.17 |

- **M1 PASS.**
- **M2 FAIL.** The r_P ratio per halving is 1.37–2.30. That is first order, not the ≥ 3.5 needed for M2-SECOND.

**M3** (lattice vs continuum, integrated over the box; the same at all three times to 3 digits):

| direction | h | energy | momentum | stress tensor (max component) |
|---|---|---|---|---|
| axis | 1 / 0.5 / 0.25 | 7.5e-3 / 1.9e-3 / 4.7e-4 | 8.9e-2 / 2.3e-2 / 5.8e-3 | 2.9e-2 / 7.3e-3 / 1.8e-3 |
| generic | 1 / 0.5 / 0.25 | 5.7e-3 / 1.4e-3 / 3.6e-4 | 5.9e-2 / 1.5e-2 / 3.8e-3 | 1.3e-1 / 3.4e-2 / 8.6e-3 |

All converge at second order (×4 per halving).

**My stated expectations:**
- M0 passes: right.
- M1 exact: right.
- M2 SECOND: **wrong.** It is first order with an O(1) coefficient. I expected the mesh-motion error (second order)
  and missed the half-cell placement, which is first order.
- M3 second order: right.

## Post hoc (labelled, `posthoc.json`, t = 6)

| direction | h | total momentum conserved | r_P as derived | r_P, momentum density re-centred to the vertex | r_E as derived | r_E re-centred |
|---|---|---|---|---|---|---|
| axis | 1 / 0.5 / 0.25 | ≤ 1e-13 | 0.53 / 0.35 / 0.17 | 0.070 / 0.022 / 0.0060 | ≤ 9e-15 | 0.56 / 0.30 / 0.16 |
| generic | 1 / 0.5 / 0.25 | ≤ 4e-14 | 0.61 / 0.29 / 0.14 | 0.19 / 0.065 / 0.018 | ≤ 8e-15 | 0.52 / 0.28 / 0.14 |

"Re-centred" is the average of the momentum density at v and v − e_a. It is a diagnostic, not a different action.

## What this means

- **The missing link is real.** A matter field whose sources come from its own action, on J1's lattice, cannot meet
  both of J1's gauge identities beyond first order. Energy and momentum need the shift in different places.
- **This explains why 34L could not see it.** 34L's prescribed dust source had no stress at O(v), so its momentum
  identity was trivially satisfied.
- **It points to the same redesign as 34J's option (b), now with a concrete target.** That target is a lattice in
  which the displacement gauge also acts on half-cell shifts:
  - the shift lives on links;
  - u is staggered the same way;
  - the edge strain map is rebuilt to match.

  Then a matter shift coupling on links could be exact for both symmetries. That is a new vacuum action (not J1), and
  34K's centred-gradient obstruction does not apply to it as such.
- **Review step 2 (energy changes gravity) should wait.** Coupling this packet to J1 would mix the physical signal
  with a first-order inconsistency of the same size as 34L's error.

## Compute

About 75 CPU-seconds:
- m0: 3 s;
- run: 73 s;
- posthoc: about 35 s.

Platform: Linux cloud, Python 3.13, numpy 2.5.

## Deviations

- **M0 numerical method.** The first M0 code used plain central differences. Their O(ε²) truncation gave 1.25e-8 at
  ε = 1e-5, just over the 1e-8 threshold. The convergence (1.2e-6 at ε = 1e-4, 1.2e-8 at 1e-5) showed it was
  truncation, not the derivation. I added Richardson extrapolation (6e-11) before the M1–M3 run. The raw rows are
  kept in `m0.json`.
- None otherwise.
