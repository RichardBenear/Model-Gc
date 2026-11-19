# EX-GC14 (EXPLORATORY): soliton route, Stage 3 – the knot as a Gc particle

Claude, 2026-11-19, at the PI's request ("Go ahead and make the knot a Gc particle"). Our construction, including the
coupling, which is a design choice.

## Answer
**Yes as a gravitational source: with the W core allowed to grow only where the frames are knotted, the knot carries
a stable Gc particle. The particle consumes G and produces the 1/r G field (within 0.1% of a lattice point sink).
Without a knot, the same particle dies. The knot keeps its degree throughout.** As a test mass (falling in a G
gradient) the result is not resolved. On this coarse lattice even a plain Gc particle barely moves, because the lattice
pins it.

## Coupling (design)
- **Knot → chemistry.** The W autocatalysis cYW² is multiplied by a gate g = clip(ρ_s/ρ_ref, 0, 1).
  - ρ = (1 − Re q)/2 is the knot density: 1 where the lift q is −1, at the knot's centre, and 0 in the vacuum.
  - ρ_s is ρ smoothed with a Gaussian of width 1.5, a "catalytic range".
  - ρ_ref = 0.017 is 0.4 of the smoothed peak of a relaxed knot.

  So W can grow only within about 2 length units of a knot.
- **Chemistry → knot.** Frame moves feel an extra energy −λ Σ W ρ (λ = 1), so the knot is pulled toward the W core.
- **Everything else is unchanged:**
  - Gc1 chemistry on the EX-GC09 lattice (h = 0.5, box ±6, deterministic);
  - frames with the EX-GC12 lift rule at T_f = 0.6, a sweep every 5 chemistry steps;
  - the initial knot is a hedgehog of radius 8 lattice units, which relaxes to a compact knot.

## Results (T = 60)
| Run | M (start → late) | Alive | Degree | Knot–core distance | Far deficit / (M × point sink) |
|---|---|---|---|---|---|
| knot + particle, seed 0 | 5.88 → 5.53 | yes | 1 throughout | 0.37 | 0.999 |
| knot + particle, seed 1 | 5.88 → 5.88 | yes | 1 | 0.14 | 1.000 |
| no knot + particle (control) | 5.88 → 0 | dies | 0 | — | — |
| knot + small W seed (0.5) | 0.40 → 0 | no particle forms | 1 | — | — |
| in gradient E = 0.02: plain Gc particle | 5.88 | yes | — | — | core moved −0.0095 (prediction μE·T = −0.031) |
| in gradient: knot particle, λ = 1 | 5.52 | yes | 1 | 0.37 | core moved +0.011 |
| in gradient: knot particle, λ = 0 | 5.86 | yes | 1 | 0.48 | core moved −0.0094 |

- **The knot holds the particle.** The core stays within 0.14–0.37 (a fraction of a lattice spacing) of the knot's
  centre. M depends on how well the knot sits on the core: when it is off-centre the gate is below 1 at the edge of the
  core, and M is 6% lower (seed 0).
- **A knot does not create a particle by itself.** A small W seed on a knot dies. As in plain Gc, the particle is one
  of two stable states and needs a finite kick to form. The knot is necessary but not sufficient.
- **Falling is not resolved.** The plain particle moved only 0.0095 against the 0.031 predicted from its continuum
  mobility (μ = −0.0255, EX008), and its motion was slowing down. That is lattice pinning at h = 0.5. The knot
  particles moved about as little, and with λ = 1 slightly the wrong way, so this lattice cannot say whether a knot
  particle falls like a plain one.

## What it means
- **One object now has all three properties:**
  - a conserved topological charge (degree 1);
  - spin-½-like and fermion-like loop classes (EX-GC13);
  - a G sink with a 1/r field, i.e. active gravitational mass.

  This holds under a chosen coupling.
- **Two consequences of the coupling:**
  - Particle number equals knot number, because a W core cannot survive without a knot.
  - Mass is set by the chemistry, not by the knot. A degree-2 knot would carry a particle whose mass we have not
    computed.
- **Still open:**
  1. Passive mass (falling) needs a finer lattice or a stronger gradient, ideally with tokens.
  2. Whether the coupling could come from one rule rather than two added terms.
  3. Two knot particles: do they attract as 1/d² and stay intact?

## Method and records
- `src/exgc14_knotgc.py` (`KnotField`: Metropolis with the extra −λWρ term; `run`), records `*.json`, and
  `src/exgc14_summary.py` → `summary.json`. The far field uses the lattice Green's function from
  `exgc09_summary.green_profile`.
- Platform: Linux cloud container, Python 3.13.18, numpy 2.5.3, scipy 1.18.1. Compute: about 0.25 core-hours.

## Deviations
1. **The first coupling killed the particle**, and those runs were discarded. Its settings were smoothing width 1 and
   ρ_ref taken from an unrelaxed hedgehog (0.124). The knot relaxed to a compact size, the gate fell to about 0.58,
   and the particle died at t ≈ 8. The catalytic range was widened (width 1.5, ρ_ref = 0.017, set from a relaxed knot)
   before the runs reported here.
2. The gradient runs could not test falling, for the pinning reason above. That was not anticipated in the log entry.
