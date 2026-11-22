# EX-GC26E3 — does internal kinetic energy fall and resist force like rest energy? (EXPLORATORY)

Claude, 2026-11-22. PI go-ahead (Richard, 2026-11-22): "Go ahead with the GC26E3 test." Cloud, ≤ 2 core-hours.
Massless Skyrme reference model with the 26E1 Nordström coupling. Not Model Gc; the construction is ours.

## Question

26E1/E2 showed that binding energy (a static, negative energy) falls with the body. Here the question is internal
*motion*. Take a spinning knot, which carries extra kinetic energy. Does that energy count for its weight and its
inertia alike, so that the knot still falls at the universal rate?

- **Continuum expectation.** By von Laue's theorem, a time-averaged bound system's Nordström source ∫(−T) equals its
  total energy, so the fall is universal.
- **On the lattice.** The source is 2(E₂ − K₂), where K₂ is the quadratic-term kinetic energy. Its time average fixes
  the prediction.

## Setup

- **Start.** The closed 26E1 B = 1 rest state (h = 0.2, L = 7, hash-checked), set isospinning at angular frequency ω.
  The spin is a rotation in the (φ₁, φ₂) plane, which for a hedgehog is a spatial rotation about z. The fall is along x,
  perpendicular to the spin axis.
- **Internal motion includes breathing.** The spin-up is sudden, so the knot also breathes. Both motions are internal,
  and the predictions use time averages.
- **Rotational inertia.** K(ω) = 0.328 ω². The kinetic fractions are 3.1% (ω = 0.35), 6.4% (0.5) and 11.8% (0.7).
- **Field, ramp, T = 5, integrator, fits:** exactly as in 26E1. Code `src/exgc26e3.py`, reusing `exgc26e1.py` unchanged.
- **Baseline.** The static B = 1 fall from the committed 26E1 records (same code and settings).

## Stage 1: null runs (no field), ω = 0.35, 0.5, 0.7

These show whether the spinning knot is stable and fix the predictions, before any fall run. For each ω, from t ≥ 1
averages:

| quantity | formula |
|---|---|
| lattice energy-coupling ratio to the static fall, R_lat | [2⟨E₂ − K₂⟩/⟨H⟩] ÷ [2E₂/E]_static |
| alternative "internal motion has inertia but no weight", R_nw | E_static/⟨H⟩ |
| alternative "weight as the action gives, but no extra inertia", R_ni | [2⟨E₂ − K₂⟩/E_static] ÷ [2E₂/E]_static |
| non-energy control ratio, R_ctl | the 26E1 control coupling with the same constant |

The numbers, and the ω values used in stage 2, will be appended to this plan and committed before stage 2 runs. An ω
is dropped if its null run loses degree or energy conservation (> 1e-8 relative), or if the knot drifts by > 1e-3.

## Stage 2: fall runs

- Energy coupling at two ω values.
- The control coupling at the larger ω.

**Readout.** The fall rate a/g from the null-subtracted baryon centroid, quadratic fit over t ≥ 1 (also t ≥ 2 and the
energy centroid). Divided by the static 26E1 value, this gives R.

**Criteria**, fixed in stage 1 from the numbers: the energy coupling is supported at an ω if

- R is closer to R_lat than to either alternative, and
- the gap between R_lat and the nearer alternative is at least twice the 26E1 spread between readouts (about 1.5%).

**Cost:** about 0.17 core-hours per run (26E1 measured), six runs in total, about 1 core-hour of cloud compute.
