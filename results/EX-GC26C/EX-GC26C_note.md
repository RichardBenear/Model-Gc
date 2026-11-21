# EX-GC26C — knot sourcing its own scalar gravity field (EXPLORATORY)

Claude, 2026-11-21. PI go-ahead (Richard, 2026-11-21): "Go ahead with the GC26C test." Plan
`spec/EX-GC26C_plan.md`, pushed before the runs (2f8c24e). Code `src/exgc26c.py` (matter kernel from 26E1 unchanged);
analysis and figure `src/exgc26c_analysis.py`; figure `paper/figures/exgc26c_field.png`.

This uses the massless Skyrme reference model with a dynamical Nordström scalar, not Model Gc. The construction is
ours.

## Conclusion

Coupled through one common action, the knot builds its own gravitational field and stays perfectly stable.

- **The field is 1/r.** It grows as an outgoing front at the speed of light. Behind the front, r·Φ is constant to
  within 2% from r = 3 to 10.
- **Its strength (active mass) equals what the knot's energy predicts.** The field strength is 0.995 of the knot's
  integrated source. In the weak-coupling limit it is 2E₂ (the lattice form of the rest energy) to within about 0.4%.
- **The bound pair is a weaker source than two separate knots.** Its active mass is 1.936 times the single knot's, not
  2. Binding energy therefore sources gravity, as 26E showed it falls with gravity.
- **Energy is conserved and the system settles.** Total energy (knot plus field) is conserved to 3e-12. The knot keeps
  its topology and does not move, and the energy exchange settles by t ≈ 4.

## Setup

| Item | Setting |
|---|---|
| Action | L = L_m[φ; W = e^{2Φ}] + Σ h³(Φ̇² − \|∇Φ\|²)/(8πG_N) |
| Field equation | Φ̈ = ∇²Φ − 4πG_N ρ, with ρ = 2W(e2 − k2)/12π² |
| Matter grid | h = 0.2, L = 7, the closed 26E1 rest states (hash-checked) |
| Field grid | the same h, half-width 17, Φ = 0 on its wall. The front reaches the wall only at t ≈ 17, so the wall plays no role up to T = 12 |
| Start | Φ = Φ̇ = 0 at t = 0 |
| Integrator | RK4 for knot and field together, dt = 0.01 |

**Validation (`validation.json`):**

- The energy rate of the coupled system along the flow is zero to 1e-10 relative.
- The site source density matches the kernel's orientation sum exactly (4e-15).
- The discrete point-source field is −1/r plus a constant, as expected for a walled box.

## Results (T = 12; fit Φ = −A/r on shells r ∈ [3, 10])

| run | A/G_N | ÷ time-mean source ∫ρ | ÷ 2E₂ (static) | ÷ E | r·Φ spread | angular spread | knot energy → field | Φ at centre |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| B = 1, G_N = 0.01 | 1.1096 | 0.995 | 0.986 | 0.923 | 1.6% | ≤ 1.8% | 0.73% of E | −0.012 |
| B = 2, G_N = 0.01 | 2.1476 | 0.995 | 0.980 | 0.929 | 2.0% | ≤ 2.5% | 1.15% of E | −0.017 |
| B = 1, G_N = 0.02 | 1.0985 | 0.995 | 0.976 | 0.914 | 1.8% | ≤ 1.8% | 1.42% of E | −0.024 |

- **A fitted constant is negligible** (≤ 7e-6), as expected with no wall reached.
- **The ÷ E column is the lattice Derrick offset.** That is 2E₂/E = 0.936, as in 26E. It is not a field error.
- **Self-gravity weakens the source slightly.** In its own well W < 1, so the source falls by about 1% per 0.01 of G_N
  (1.125 → 1.110 at G_N = 0.01, and → 1.096 at 0.02). This is linear in G_N. Extrapolating the two couplings to
  G_N → 0 gives A/G_N ≈ 0.996 × 2E₂.
- **The B = 2 angular spread is slightly larger.** The torus is not spherical, and its quadrupole falls off as 1/r³.
- **The knot is stable.** The degree is unchanged in every run, the centroid moves ≤ 2e-7, and total energy is
  conserved to ≤ 3e-12.

**Active mass, B = 2 ÷ B = 1:**

| | value |
|---|---:|
| Measured from the field | 1.9356 |
| Source-integral ratio | 1.9370 |
| Lattice prediction 2E₂ ratio (G_N → 0) | 1.9472 |
| Energy ratio E₂/E₁ (continuum expectation) | 1.9217 |
| No-binding alternative | 2.0000 |

The measured value lies between the lattice prediction and the energy ratio and is 3.2% below 2. The extra 0.6% below
the lattice prediction is the stronger self-gravity of the heavier pair (Φ at its centre is −0.017 against −0.012).

## What this does and does not show

- **The far-field charge equals the integrated source by construction** (Gauss's law). What the run adds is that the
  coupled *dynamics* realises it: a stable knot, a clean 1/r field formed causally, exact energy bookkeeping between
  knot and field, and a source fixed by the knot's own fields, binding included. Together with 26E, active and passive
  gravitational mass both track the rest energy, for the single knot and the bound pair.
- **One resolution** (h = 0.2) and weak coupling only.
- **The knot is at rest.** Moving knots with their field (26D) and two-body attraction (26F) are not yet tested.
- **Scalar gravity limits remain**: no light bending.

## Compute

1.71 core-hours of the 2 core-hour cloud cap: three runs of about 2040 s each, plus validation. Platform: Python
3.13.18, NumPy 2.5.3, SciPy 1.18.1, Linux. End states (knot and field) are kept locally in `local_raw/EX-GC26C`, with
SHA-256 hashes in the records.

## Next (PI decision)

- **26D:** boost the coupled knot and its field, and check that it coasts without drag. The field must move with the
  knot, which tests the absence of a rest frame directly.
- **26F:** two knots attracting through the self-generated field.

Both need larger field boxes, so they are better suited to the PC.
