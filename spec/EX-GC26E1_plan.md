# EX-GC26E1 — single knot vs bound pair in a scalar gravity field (EXPLORATORY)

Claude, 2026-11-21. PI go-ahead (Richard, 2026-11-21): "I have attached ChatGPT's response to the outline. If you
agree with it, go ahead with the EX-GC26E1 test". Based on `spec/EX-GC26E_outline.md` as amended by ChatGPT's response:
B = 2 first, rotating and breathing states deferred to 26E3, both centroids retained, two resolutions give a
*candidate* violation at most. Base commit 1e40f94. Massless Skyrme reference model (c = 1), not Model Gc; the lattice
construction and coupling are ours. Cap 2 core-hours (cloud), one numerical thread per process, at most 4 workers.

## Model and coupling

- Inherited 26B2 eight-orientation mirror-symmetric stencil, RK4 with projection, dt = 0.05h, fixed vacuum walls,
  h = 0.2, L = 7. **Same box for both species.**
- Nordström coupling from the action with g = e^{2Φ}η. Only the quadratic term changes: at each site i, its kinetic
  and gradient parts are multiplied by W_i = e^{2Φ(x_i,t)}. The quartic (Skyrme) term is unchanged, because it is
  conformally invariant in 3+1D.
- Φ = λ(t) g x along the axis, with g = 0.0125 (F/E ≈ 0.025, as in 26B3). λ(t) is the 26B3 half-cosine ramp from 0
  to 1 over t ∈ [0, 1]; the Ẇ·v term from the ramp is included.
- **Non-energy control:** the 26B2/B3 squared-deficit coupling V = −f λ(t) Σ x (1−φ₀)², with **one** constant f for
  both species. f is chosen so that B = 1 has the same predicted force as under the energy coupling.

## Species

- **B = 1:** hedgehog relaxed from scratch.
- **B = 2:** torus relaxed from the rational map R(z) = z² with the same orientation convention (degree −2).
- Both use the 26B2 optimizer (L-BFGS, 300 iterations; up to 600 for B = 2). Report the residual, E₂, E₄, the
  geometric degree, Q = Σ(1−φ₀)², and the binding energy E_bind = 2E(B=1) − E(B=2). Also report the fraction of
  energy within 1.5 of the walls, as a box check.

## Predictions, written before any fall run

- Energy coupling. Lattice prediction for each species: a/g = −2E₂/E (static, first order in g).
  Continuum: a/g = −1 for both.
- Control coupling: a = −f·Q/E for each species, with the same f.

## Runs (T = 5)

1. B = 1 energy coupling
2. B = 2 energy coupling
3. B = 1 control
4. B = 2 control
5. B = 1 no field
6. B = 2 no field

The no-field runs are subtracted. Priority order: 1, 2, 6, 4, 3, 5. Stop before the cap and record any omissions.

## Readouts

- Fall acceleration from a quadratic fit after the ramp (t ≥ 1), with t ≥ 2 as a window check. Two centroids,
  reported separately and in agreement or not:
  - the **baryon-density** centroid (eight-orientation determinant density);
  - the **energy** centroid (coupled energy density).
- Momentum gained versus the predicted impulse.
- Geometric degree at the start and end, transverse momentum, energy conservation after the ramp, and the
  no-field drift.

## Primary comparison

R = (a/g)_{B=2} ÷ (a/g)_{B=1} under the energy coupling, against:

- the lattice prediction (2E₂/E)_{B=2} ÷ (2E₂/E)_{B=1};
- the continuum value 1.

Under the control coupling, check that the measured ratio follows (Q/E)_{B=2} ÷ (Q/E)_{B=1}.

Reference error baseline: 26B3, a few percent on force response, below 1% on direct mass. Interpretation:

- **Consistent with universality:** R within about 1% of the lattice prediction and of 1.
- **Candidate violation:** an R difference well beyond the baseline. A candidate needs 26E2 (h = 0.15) before any
  claim.
- **Control fails:** if the control does not follow its predicted Q/E ratio, the method cannot see violations and the
  energy-coupling result is reported as not interpretable.

## Validation (before production)

- With W ≡ 1, the new kernel matches the 26B2 kernel bitwise or to 1e-13.
- An independent NumPy Euler–Lagrange sum with W and Ẇ agrees with the kernel to 1e-10.
- Finite-difference directional gradient of the static energy with W.
- Exact energy conservation rate with a static W (zero to 1e-7).
- Mirror check transverse to the field direction.
- Vacuum zero.

Raw states stay local with SHA-256 hashes. Compact time series, code and figure are committed.
