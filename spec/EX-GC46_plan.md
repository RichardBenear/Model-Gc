# EX-GC46 — heavy pairs with all four upper bands, and the small-g limit

2026-11-27 | EXPLORATORY | Claude | Plan only. PI request in thread: "the 45 review is available. We should consider the test suggested there" (2026-11-27 04:06 UTC). PI go-ahead 2026-11-27 04:09 UTC ("go ahead and close 43, 44, 45 and write the plan for 46 and then execute the runs"), given before this plan was pushed. The plan is pushed before any run. Follows the [Codex review of GC45](../review/EX-GC45_codex_review.md), tests 1 and 2. Test 3, a finite-range interaction with a pre-fixed strain law, is deferred until these two show whether the contact result is robust. Cap 2 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

## Why

45 found, in a heavy-band-projected model, that a bound pair of heavy carriers has a binding contribution to its spatial charge of order its binding energy (η = δQ_a/E_b from −2.7 to +4.9). Codex's review accepts the broad result but points out:
- the "|η| > 0.1 everywhere" clause fails with all four upper bands (η ≈ −0.04 at the middle point), and the projected solver's error is largest where η is small;
- m₀ and g changed together, so the sign flip is not attributable to g;
- the binding part can be tiny next to the pair energy even when η > 1, so total Q_a/ΔE must be reported beside it;
- the smallest-g bound pairs were not measured, so "dust at finite mobility" is shown only for a single carrier.

η here is minus the logarithmic derivative of the binding energy with respect to the hopping scale, at fixed U. It is a strain sensitivity under a contact law. It is not an observational Nordtvedt parameter. It changes the active source through Q_a.

## Model

As 45: the 40A network with one heavy mode per site, two species, lapse-homogeneous weights, hop ∝ 1/a, on-site contact V = −U Σ n_h1 n_h2 (w = 0). The valence sea stays frozen (not changed by this study).

**Four-band solver.** Each carrier may occupy any of the four 7-band bands above μ (indices 3–6). At total momentum K the pair basis is (k, b₁, b₂), with D = e_b1(k+K) + e_b2(−k) and a rank-1 contact amplitude Φ(k, b₁, b₂) = conj(u_b1(k+K)_6)·conj(u_b2(−k)_6). The bound state solves U S(E)/N = 1, with S = Σ|Φ|²/(D − E). This is the 45 solver with a larger basis. The free-pair threshold is the minimum of D over all (k, b₁, b₂), so E_b is measured against the four-band threshold. Q_a = t ∂E/∂t at fixed U by central difference; δQ_a = Q_a(pair) − Q_a(threshold).

## Parts

### 46A — four-band solver at 45's nine points (Codex test 1)

- Points: (m₀, g) = (0.3, 0.2), (0.5, 0.1), (0.7, 0.05), at E_b/W_h = 5 / 15 / 30%, L = 32. W_h is the heavy-band width, as in 45.
- Two comparisons: (i) fixed U, taken from 45, to measure projection error; (ii) U re-tuned in the four-band model so E_b equals the same fraction of W_h, to compare across carrier models at matched E_b/W.
- Record: ΔE, E_b, δQ_a, η, total Q_a/ΔE, binding share E_b/ΔE.

### 46B — the small-g pair limit (Codex test 2)

- m₀ ∈ {0.5, 0.7}, g ∈ {0.0125, 0.025, 0.05, 0.1}, E_b/W_h = 0.15 and 0.30, L = 32, four-band solver with U tuned per point.
- Record, separately: total Q_a/ΔE; the binding part δQ_a, η, and E_b/ΔE; the pair's mobility, taken as the dispersion curvature t_pair = [E₀(K) − E₀(0)]/|K|², with K = (2π·2/32, 0, 0); and the single carrier's t_eff from 45B for comparison.
- Two m₀ values at the same g separate hybridisation from the on-site scale.

## Predictions (scored)

- **P1 (46A).**
  - At fixed U the four-band E_b differs from the projected value by less than 5%, and δQ_a changes by more than 5% at the (0.5, 0.1) point (as Codex found).
  - At matched E_b/W the signs of η follow 45: negative at (0.3, 0.2), positive at (0.7, 0.05) for all three fractions.
  - |η| at (0.5, 0.1), 15% is below 0.1.
  - The spread of η across the nine points exceeds 20%.
- **P2 (46B).**
  - Total pair Q_a/ΔE at g = 0.0125 is below 5e-3 at both m₀.
  - Binding share E_b/ΔE is below 1e-3 at g = 0.0125.
  - t_pair scales as g² within 20% (g = 0.0125 vs 0.025).
  - η has a finite limit as g → 0: |η(0.0125)/η(0.025) − 1| < 20% at each (m₀, fraction), because δQ_a and E_b both shrink as g².
- **P3 (descriptive).** The sign of η at fixed g for the two m₀ values, to show whether the sign belongs to g or to the on-site scale.

If η has no finite limit as g → 0, or changes sign inside the sweep, the note says so.

## Checks

- **C1.** Q_N = ΔE for pairs (uniform-N finite difference with U scaled), within 1e-8.
- **C2.** With the basis restricted to band 3, the four-band solver reproduces 45's projected solver within 1e-12. Hellmann–Feynman: U ∂E/∂U equals ⟨V⟩ within 1e-8.
- **C3.** A dense four-band pair matrix at L = 6 equals the rank-1 solver within 1e-10.
- **C4.** The ordering: E_b ≥ 0 and the bound energy lies below the four-band threshold at every point.

## Outputs

- `src/exgc46.py` and `results/EX-GC46/`.
- A separate recomputation, `src/exgc46_recompute.py`: an independently built 7-band H (opposite phase convention), a dense four-band pair matrix at L = 4 and 6, and η at one 46A point, with no use of the driver's solver.
- Not done here (deferred): finite-range interaction, and any change to the frozen valence sea.
