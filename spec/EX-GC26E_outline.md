# EX-GC26E — outline (DRAFT for PI discussion; not a registered plan)

Claude, 2026-11-21, after EX-GC26B3. Whoever runs 26E writes the pre-run plan (`spec/EX-GC26E_plan.md`) and the
explore/LOG.md entry from this outline. The massless Skyrme model is the reference model here, not Model Gc. The
lattice construction is ours.

## Question

Do different particle configurations fall at the same rate in a scalar (Nordström) gravitational field? Coupling to
energy predicts yes: by von Laue's theorem, every stable bound configuration has gravitational charge equal to its
total energy, binding and internal motion included. This is the 3D version of the EX-GC25 1+1D check, and the direct
answer to the Model Gc failures EX-GC01 to EX-GC04 (types and bound composites fell at different rates).

## Proposed order change (PI decision)

The original sequence was 26C (dynamic G field), 26D (boost with field), 26E (universality), 26F (two-body). I
suggest running 26E first, with a **fixed external field** (test-body limit). Reasons:

- It answers the main question without first building and qualifying the field dynamics.
- It reuses the qualified B2/B3 solver almost unchanged.
- 26C/D then add the self-generated 1/r field, and can reuse the 26E coupling code.

## Coupling (common action)

Nordström metric g = e^{2Φ}η inserted into the Skyrme action. In 3+1D the Skyrme (quartic) term is conformally
invariant, so only the quadratic term changes: both its kinetic and its gradient parts are multiplied by e^{2Φ(x)},
site by site. The quartic term is untouched. Use a weak uniform gradient Φ = g·x, ramped on smoothly as in 26B3, along
an axis and along the diagonal.

## Particle configurations ("species")

1. **B = 1 hedgehog.** The reference.
2. **B = 2 torus.** A bound pair with about 4–5% binding energy in the continuum. This is the EX-GC04 analogue: does
   binding energy fall with the body?
3. **B = 1 isospinning** at low spin. Internal kinetic energy, stationary rather than static.
4. **B = 1 breathing** (optional). Internal oscillation energy; the theorem holds only on time average.
5. **Positive control:** the same species with a deliberately non-energy coupling, e.g. Φ times ∫(1−φ₀)². The
   predicted fall rate is Q/E, which differs by species. This shows the test can detect a violation, as in EX-GC25.
   Note that a quartic-only coupling is *not* a valid control: by Derrick's theorem E₄ = E/2 for every static
   solution, so it would also look universal.

## The main known systematic: the lattice Derrick imbalance

On the lattice the coupling acts only on E₂, so the gravitational charge of a static state is 2E₂, not E. In the
continuum E₂ = E₄, so 2E₂ = E. On our lattice E₂/E₄ = 0.880 (h = 0.2) and 0.930 (h = 0.15), so:

- the predicted fall rate is a/g = 2E₂/E = 0.936 and 0.964;
- this 4–6% offset is larger than the 26B3 few-percent force systematics;
- it may differ between species.

So for every species, 26E should:

- compute the predicted lattice charge from the static (or time-averaged) state before running, as a cheap static
  number;
- test the measured fall against that prediction (dynamics check);
- compare species with each other, at h = 0.2 and 0.15 (physics check).

A real violation persists under refinement. A lattice artefact shrinks roughly as h².

## Readouts (same conventions for all species)

- Fall acceleration from the **baryon-density centroid**: topological, so insensitive to the tail (26B3 showed
  energy-centroid regional sensitivity). Energy centroid reported alongside.
- Momentum gain versus predicted impulse g·(charge), and P/v inertial mass, kept separate as in 26B3. Fit only after
  the ramp, with the same windows for every species.
- Geometric degree at start and end, the energy ledger including the source term, transverse momentum, and a no-field
  control per species.

## What counts as a result

Primary: the ratio (fall rate of species X) ÷ (fall rate of B = 1) at the same h, with its change from h = 0.2 to
0.15.

- **Universal:** within the 26B3 error baseline at both resolutions, and moving toward 1 under refinement.
- **Violation:** a difference well outside the baseline that does not shrink under refinement.
- **Control:** the positive control must show its predicted spread.

Anything in between is reported as unresolved. No precision gate. EXPLORATORY.

## Compute

B = 2 needs a larger box (L ≈ 9) and its own relaxation from rational-map initial data, plus validation. Estimate:
about 4 species × 2 resolutions × (field + control), roughly 6–12 core-hours. That exceeds the 2 core-hour cloud cap,
so it either needs the PI's explicit cap increase or should run on the PI's PC (about 33 to 35 of the 40 approved
core-hours remain). A cheap first pass: B = 1 and B = 2 at h = 0.2 only, plus the static 2E₂/E predictions for all
species.
