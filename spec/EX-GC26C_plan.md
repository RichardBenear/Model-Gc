# EX-GC26C — knot sourcing its own scalar gravity field (EXPLORATORY)

Claude, 2026-11-21. PI go-ahead (Richard, 2026-11-21): "I sign off on GC26E2. Go ahead with the GC26C test." Base:
26E1/E2 code (closed). Massless Skyrme reference model with a dynamical Nordström scalar, c = 1. Not Model Gc; the
construction is ours. Cloud cap 2 core-hours, one thread per process.

## Question

Does a knot coupled to a dynamical scalar field G (here Φ, with G = e^{2Φ} acting as the metric factor) build its own
static 1/r field and stay stable? Is the far-field charge, its *active* gravitational mass, equal to what the common
action predicts, both for B = 1 and for the bound pair B = 2?

This is the 26E counterpart: there the knots responded to an imposed field (passive mass); here they source it
(active mass).

## Common action

L = L_m[φ; W = e^{2Φ}] + Σ_i h³ (Φ̇² − |∇_h Φ|²) / (8π G_N)

- L_m is the 26E1 matter action. Only the quadratic Skyrme term carries W.
- The field equation that follows is Φ̈ = ∇_h²Φ − 4π G_N ρ, with ρ_i = 2W_i (e2_i − k2_i)/12π², where e2 and k2 are
  the site gradient and kinetic quadratic densities. In the static case ∇²Φ = 4πG_N ρ: an attractive potential well,
  Φ → −G_N M_g / r, with active mass M_g = ∫ρ (lattice 2E₂, continuum E by von Laue).
- The matter feels W = e^{2Φ} and Ẇ = 2Φ̇W through the validated 26E1 kernel.
- Total conserved energy: H_m(W) + Σ h³ (Φ̇² + |∇Φ|²)/(8πG_N).

## Grids and run

| Item | Setting |
|---|---|
| Matter grid | h = 0.2, L = 7 (the closed 26E1 rest states, hash-verified) |
| Field grid | Φ on a larger box, L_Φ = 17, same h, matter grid embedded at its centre. Φ = 0 on its outer wall |
| Start | Φ = Φ̇ = 0 at t = 0, coupling on from t = 0, time-independent action |
| Wave front | The field grows as an outgoing front. In the continuum, Φ = −G_N M/r exactly behind it (r < t); the front reaches the outer wall only at t ≈ 17, so the wall plays no role up to T = 12 |
| Coupling | G_N = 0.01, so Φ at the core is about −0.02 (W ≈ 0.96): weak self-gravity |
| Time | T = 12, dt = 0.05h, RK4 for (φ, v, Φ, Φ̇) together |

**Runs:** B = 1 and B = 2. Plus, budget permitting, B = 1 with G_N = 0.02 as a linearity check.

## Predictions (stated before the runs, from the static records)

- M_g = 2E₂ at t ≈ 0: B = 1 1.1254, B = 2 2.1914; ratio 1.947.
- That ratio differs from the "binding does not source" value of 2. The continuum value (equal to the E ratio) is
  1.922.
- The far field r·Φ(r) is constant (= −G_N M_g) for r from about 3 to t − 2. In the time-dependent case M_g is the
  source at the retarded time.
- The knot stays put and keeps its degree. It relaxes slightly into its own well: the expected self-energy is of order
  G_N E²/R, about 1%.

## Readouts

| Readout | How measured |
|---|---|
| Active mass | Shell averages of Φ(r) at T. Fit −A/r on r ∈ [3, 10]; A/G_N against ∫ρ (current and retarded) and against 2E₂, E |
| Shape | Constancy of rΦ(r) in r. Angular spread on shells: monopole vs quadrupole, especially for the B = 2 torus |
| B = 2 vs B = 1 | Active mass ratio, against 1.947 (lattice), 1.922 (E ratio) and 2 (no binding) |
| Stability | Geometric degree at start and end, energy centroid drift, matter energy H_m(t) and field energy over time |
| Conservation | Total energy to ≲ 1e-6 relative while the front is inside the box |

## Validation (before production)

- Exact energy-rate zero for the coupled system on a small grid, analogous to the 26E1 ledger check.
- Discrete Laplacian checked against the point-source Green's function in the interior.
- Matter kernel unchanged (26E1 validation stands).

## Limits

- **Gauss's law:** the far-field charge equals the integrated source by construction. The test is that the *dynamics*
  realises it: stable, 1/r, energy-conserving, and with ρ set by the knot's own fields, including binding.
- **Scalar gravity:** no light bending (known).
- **One resolution.**
