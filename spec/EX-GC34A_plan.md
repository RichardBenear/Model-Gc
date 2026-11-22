# EX-GC34A — one shared link structure: do all fields share one light cone? (EXPLORATORY)

Claude, 2026-11-22. PI go-ahead (Richard, 2026-11-22): "Go ahead with the 34A plan." Cloud, cap 2 core-hours.
First stage of EX-GC34 (outline `docs/EX-GC34_outline.md`). Our construction, not LaViolette's Model G. No Skyrme
term in this stage: the knot itself moves onto the shared links in 34B.

## Question

Suppose every field can move only through one shared lattice structure: a site part (clock and volume) and nine link
weights per site. Do physically different fields then share one light cone, on flat, anisotropic and curved
backgrounds, without any per-field speed coefficient? And does the solver reproduce the textbook light-bending
multiples (0, 1, 2 times "Newtonian") on prescribed backgrounds? 34C needs that machinery to measure light bending
around a knot.

**This is expected to pass by construction for the linear fields.** Its job is to validate the code and the
metric-to-link mapping, to measure lattice and time-step errors, to test a nonlinear field, and to build and validate
the deflection readout before 34C. A GC31A-style control (a field with its own link weights) must fail.

## Shared structure

Spatial metric γ_ij and lapse N, mapped once onto the lattice (cubic, spacing h):

- site kinetic coefficient Ms = √γ / N; site potential coefficient Vs = N √γ;
- link weights w_d for the 3 axis and 6 face-diagonal directions, at link midpoints, with second moment
  Σ_d w_d d dᵀ = A = N √γ γ⁻¹. The 3 left-over weight combinations (null space) are fixed by β_ij = (A_ii + A_jj)/12,
  which gives the fourth-moment isotropic stencil (axis 1/3, diagonal 1/6) on flat space. The alternative β_ij =
  (A_ii + A_jj)/16 is the null-space test (same continuum metric, different lattice artefacts).

Every field's action is built from (Ms, Vs, w) only:
L = Σ_s h³ [½ Ms ḟ² − ½ Vs μ² f²] − Σ_s Σ_d h³ w_d (f(s+d) − f(s))² / (2h²).

## Fields

| id | field | why |
|---|---|---|
| F1 | massless scalar | "light" |
| F2 | massive scalar, μ = 1 | massive matter |
| F3 | two components with field-space metric G = [[1, 0.4], [0.4, 1]] in *both* kinetic and gradient terms, plus internal mass term ε(u₁ − u₂)², ε = 0.6 | internal structure and mixing; modes: symmetric (massless) and antisymmetric (μ = √(2ε/(1−g)) = √2). Contrast GC31A, where the mixing sat in the gradient term only |
| F4 | O(4) sigma model, elliptically polarised packet π = 0.5 Re[z ê], ê ∝ (1, 0.6i, 0.4) | nonlinear (the knot's field without the Skyrme term) |
| F5 | control: massless scalar with its *own* link weights, 1.10 × the shared ones | must fail (≈ 4.9% faster) |

Packets: Gaussian envelope σ = 2, carrier |k| = 3 (kh = 0.6), small amplitude 1e-3 for linear fields; initial
velocity from the local constant-coefficient lattice dispersion so the packet moves one way. Integrator: velocity
Verlet, dt = 0.25h; sigma model renormalised and projected each step.

## Predictions (separate code path: `src/exgc34a_rays.py`)

Ray tracing with H(x,k) from the exact lattice dispersion with local coefficients, including the exact leapfrog
time-step correction (sin(W dt/2) = ω dt/2). Rays sample the packet's Wigner function on a 4⁶ Gauss–Hermite grid,
weighted by energy (∝ ω²). The predicted readout is the weighted mean ray position, fitted exactly as the field
centroid is. The ray code re-implements the background formulas and weight split; it does not import the solver.
Flat-space method check before registration (F1 only, the only run made before this plan): measured 0.94067 vs
predicted 0.94112 (−0.05%), energy drift 2e-5.

## M1: uniform backgrounds (periodic box 80³, h = 0.2, T = 5, start at origin)

Backgrounds: B0 flat; B1 uniform anisotropic γ = Q diag(1.25, 1, 0.8) Qᵀ (Q = rotation 30° about z, then 20° about
x), β/12; B1′ same metric with β/16. Directions k̂: (1,0,0), (1,1,1)/√3, (0.8, 0.5, 0.33) normalised.
Readout: centroid velocity of energy density, linear fit t ∈ [1, 5]; F3 per mode.

**Pass, every background and direction:**
- F1, F2, F3 (both modes): |v − v_ray| / |v_ray| ≤ 0.3% and angle between v and v_ray ≤ 0.003 rad;
- F4: same with 1% and 0.01 rad (envelope readout of a nonlinear packet);
- F5: deviates from the shared-structure prediction by ≥ 3%, and agrees with its own-structure prediction within 0.3%;
- B1 vs B1′: each agrees with its own prediction; their difference is reported as the size of lattice (null-space)
  artefacts at kh = 0.6.

## M2: curved backgrounds, deflection (box x ∈ [−24, 24], y ∈ [−6, 14], z ∈ [−6, 6], h = 0.2)

Prescribed weak field Φ = −A/√(r² + a²), A = 0.05, a = 2:
- B2 lapse only: N = 1 + Φ, γ_ij = δ_ij (PPN γ = 0; "Newtonian", light deflection ½ of GR);
- B3 GR-like: N = 1 + Φ, γ_ij = e^{−2Φ} δ_ij (γ = 1);
- B4 Nordström: N = e^Φ, γ_ij = e^{2Φ} δ_ij (conformally flat; light deflection 0 in the continuum).

Packet from (−18, 6, 0) along +x, T = 36. Fields F1, F2, F3, F4 in each background (12 runs). Readout: θ from the
fitted centroid velocity in early (t ∈ [2, 8]) and late (t ∈ [28, 36]) windows; deflection α = θ_late − θ_early,
identical windows for rays.

**Pass:** for every field and background, |α − α_ray| ≤ 0.05 × |α_ray(B2, same mass)|. Massless modes are compared
with massless rays, massive modes with rays of their own μ. Reported: α(B3)/α(B2) (continuum weak-field expectation
2) and α(B4)/α(B2) (expectation 0), for massless fields.

## Separate recomputation (`src/exgc34a_check.py`)

From saved final fields (local, not committed; SHA-256 manifest committed): the propagation direction from the
energy-weighted mean group velocity of the final field's spatial Fourier spectrum (massless M2 runs and all M1
runs), compared with the centroid readout. Shared dependencies: the solver's fields and background builder.

## Stop and report

The study stops and reports if any linear field fails, since that would be an implementation error, not physics.
A nonlinear-field (F4) failure is reported as a result. No retuning of tolerances after the data.

## Compute

About 1.3 core-hours estimated (M2 dominates, the sigma model is about 4× a scalar); cap 2. Platform recorded.
Outputs: `results/EX-GC34A/` (records, note), figure in `paper/figures/`.
