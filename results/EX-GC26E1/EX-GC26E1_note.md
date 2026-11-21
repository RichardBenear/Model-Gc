# EX-GC26E1 — single knot vs bound pair in a scalar gravity field (EXPLORATORY)

Claude, 2026-11-21. PI go-ahead (Richard, 2026-11-21); **PI CLOSED 2026-11-21**. Plan `spec/EX-GC26E1_plan.md`, pushed before any run
(493cba7). Code `src/exgc26e1.py` and `src/exgc26e1_stencil.cpp`; analysis and figure `src/exgc26e1_analysis.py`;
records in `results/EX-GC26E1/`; figure `paper/figures/exgc26e1_fall.png`.

This uses the massless Skyrme reference model, not Model Gc. The lattice construction and the coupling are ours. The
stencil, integrator and optimizer are ChatGPT's 26B2/B3 code, unchanged.

## Conclusion

The bound pair (B = 2) and the single knot (B = 1) fall at the same rate to within about 2%. The ratio of their
fall rates is 0.995 to 1.020, depending on which centroid and fit window is used.

The deliberately non-energy control coupling, in the same runs, gives a 35–45% difference, close to its predicted
1.425. So the method clearly sees a violation when one is built in, and with the energy coupling there is none
beyond the lattice error.

What this does not yet settle is the specific question of whether **binding energy** gravitates. The binding energy is
only 3.9% of the pair's energy. If it did not gravitate, the ratio would be 1.041. The baryon-centroid values (0.995,
1.003) are about 4% from that. The energy-centroid value (1.020) lies halfway between 1 and 1.041. The finer grid
(26E2, h = 0.15) is needed to decide this.

## Setup

- **Grid:** h = 0.2, box ±7, the same box for both species.
- **Coupling:** Nordström metric e^{2Φ}η in the action. Only the quadratic (sigma-model) kinetic and gradient terms are
  weighted by W = e^{2Φ}; the Skyrme term is conformally invariant and unchanged.
- **Field:** Φ = λ(t)·g·x with g = 0.0125, ramped on over t ∈ [0, 1] with the 26B3 half-cosine; T = 5.
- **Control:** the 26B2/B3 squared-deficit coupling, with one constant f for both species. f is chosen so that B = 1
  has the same predicted force as under gravity.
- **Null subtraction:** each species has a no-field run, which is subtracted from its forced runs.

**Validation (`validation.json`):**

| Check | Result |
|---|---|
| W = 1 against the 26B2 kernel | equal (≤ 6e-17) |
| Independent NumPy Euler–Lagrange sum with W, Ẇ and the control force | 3e-14 |
| Directional derivative of the static energy | 1e-9 relative |
| Energy rate with a static W | 6e-11 |
| Energy ledger under a time-dependent W (dH/dt = −∂L/∂t) | 5e-5 relative, quadrature-limited |
| Mirror symmetry and vacuum | zero |

## Static states and predictions (stated before the fall runs)

| | E | E₂/E₄ | Q = Σ(1−φ₀)² | degree | lattice prediction a/g (energy) | prediction a/g (control) |
|---|---:|---:|---:|---:|---:|---:|
| B = 1 | 1.20240 | 0.880 | 0.0685 | −1 | −0.9360 | −0.9360 |
| B = 2 | 2.31063 | 0.902 | 0.1875 | −2 | −0.9484 | −1.3338 |

- Binding energy: 2E₁ − E₂ = 0.0942, or 3.9% (continuum about 4.3%).
- B = 2 relaxed from the rational map z² and converged (residual 2e-8). B = 1 stopped at the 300-iteration cap
  (residual 3e-6), as in 26B2.
- Energy within 1.5 of the walls: 0.14% (B = 1) and 0.12% (B = 2).
- The lattice prediction is a/g = −2E₂/E (the Derrick imbalance). The continuum value is −1 for both.

## Fall results (null subtracted)

| run | a/g, baryon centroid, t ≥ 1 / t ≥ 2 | a/g, energy centroid, t ≥ 1 / t ≥ 2 | measured ÷ prediction (baryon, t ≥ 1) | momentum gain ÷ impulse |
|---|---|---|---:|---:|
| energy, B = 1 | −0.920 / −0.937 | −0.887 / −0.890 | 0.983 | 0.964 |
| energy, B = 2 | −0.916 / −0.940 | −0.905 / −0.907 | 0.966 | 0.968 |
| control, B = 1 | −0.949 / −0.984 | −0.896 / −0.904 | 1.014 | 0.964 |
| control, B = 2 | −1.277 / −1.324 | −1.296 / −1.306 | 0.957 | 0.969 |

**Ratio B = 2 ÷ B = 1:**

| | baryon, t ≥ 1 | baryon, t ≥ 2 | energy, t ≥ 1 | energy, t ≥ 2 | lattice prediction | continuum |
|---|---:|---:|---:|---:|---:|---:|
| energy coupling | 0.995 | 1.003 | 1.020 | 1.020 | 1.013 | 1 |
| control coupling | 1.345 | 1.345 | 1.445 | 1.444 | 1.425 | n/a |

## Controls and quality

- **Null runs:** centroid drift ≤ 1e-8, momentum drift ≤ 2e-9 E. Subtraction makes no visible difference.
- **Topology:** every run keeps its degree (−1 or −2) at both ends in all three preimage checks.
- **Energy:** conserved after the ramp to ≤ 5e-14 (H, or H + V for the control). Transverse momentum is zero to
  rounding.
- **Momentum gain ÷ predicted impulse:** 0.964–0.969 for all four runs. This is the same factor as 26B3 at h = 0.2
  (0.965), so it is a common lattice effect, not a species effect.

## Caveats

- **The baryon centroid rings after the ramp.** Its fit residual is 3–7e-4 at t ≥ 1, against ≤ 1e-4 for the energy
  centroid, and its a/g moves by 2–4% between fit windows. The energy centroid is smooth and stable between windows,
  but it includes the tail. The two centroids differ by 3–6% for a single run. **The B = 2 ÷ B = 1 ratio is much less
  sensitive than that (0.995 to 1.020)**, because the common lattice effects cancel, as the comparative design
  intended.
- **Measured ÷ lattice prediction is 0.96–1.01.** That is within the 26B3 force-response baseline at h = 0.2 (F/a 4%
  high), but B = 2 comes out about 1.7% lower than B = 1. One resolution cannot say whether that is lattice error or
  physics.
- **This is a test body in an imposed field.** It does not test the knot sourcing the field (26C) or two-body forces
  (26F).
- **Axial direction only**, with T = 5. Displacements are 0.11 (energy coupling) to 0.16 (control): small compared
  with the knot size, but up to 0.8 lattice spacings.

## Next

26E2 should repeat the same four forced runs and two nulls at h = 0.15. If the ratio stays near 1, or moves toward 1,
while the lattice prediction 2E₂/E converges, binding energy gravitates. If it moves toward 1.041, it does not.

Estimated cost is about 3–4 core-hours (h⁻⁴ scaling from the 1.04 core-hours used here), so it needs a cap increase
or your PC. 26E3 (spinning and breathing) follows.

## Compute

1.04 core-hours under the 2 core-hour cap:

| Item | CPU time |
|---|---|
| Statics | 218 s |
| Six runs | 575–600 s each, 3527 s in total |
| Validation | under 1 s |

Cloud: Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1, Linux. Raw rest and end states are kept locally in
`local_raw/EX-GC26E1`, with SHA-256 hashes in the records.
