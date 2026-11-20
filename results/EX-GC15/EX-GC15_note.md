# EX-GC15 (EXPLORATORY): does the knot particle fall?

Claude, 2026-11-20, at the PI's request ("go ahead and do the finer lattice / tokens run."). Follow-up to EX-GC14.
Our construction.

## Answer
**Yes. On the finer lattice (h = 0.25) the knot particle falls toward lower G at exactly the speed of a plain Gc
particle, and that speed matches the continuum prediction within 2%. The knot keeps up with its core. But this only
works when the knot is pulled toward the core (λ = 1). Without that pull, the knot shrinks to lattice size, its gate
weakens and the particle dies.**

| Run (h = 0.25, gradient E = 0.02, T = 60) | Late drift of the core (t ≥ 20) | M | Knot |
|---|---|---|---|
| plain Gc particle | −4.998 × 10⁻⁴ per time unit | 5.515 → 5.512 | none |
| knot particle, λ = 1 | −4.998 × 10⁻⁴ (identical to 4 digits) | 5.515 → 5.512 | degree 1 throughout; centre within 0.05 of the core |
| knot particle, λ = 0 | — | dies by t ≈ 6 | degree 1, wanders 0.2–0.4 |
| continuum prediction μE (EX008: μ = −0.0255) | −5.1 × 10⁻⁴ | | |

- **The finer lattice removes the pinning.** At h = 0.5 (EX-GC14) the plain particle moved only 0.0095 in 60 time
  units and was slowing. At h = 0.25 it moves steadily, 0.026 by t = 54, at the continuum speed within 2%.
- **The knot particle falls like the plain one.** Its core drift is identical to the plain particle's to four digits,
  and so is M. The knot sits on the core the whole time. Active and passive gravitational mass are therefore those of
  the Gc particle the knot carries.
- **Why λ = 0 dies.** A relaxed knot is a few lattice spacings across, so on the finer lattice it is half as big in
  physical units (about 0.5). Its smoothed density then sets the gate to only about 0.45, and the core starves. With
  λ = 1 the pull toward W swells the knot to about 1500 mismatched edges over the core, the gate stays at 1, and the
  particle lives.

  So in this coupling **the knot's size is set by the chemistry, through λ, not by the lattice.** That is the sensible
  way round, since the lattice spacing should not decide whether a particle exists, but it means λ is essential and not
  optional.

## What it means
- **One object now has every property tested so far:**
  - a conserved topological charge;
  - spin-½-like and fermion-like loop classes (EX-GC13);
  - a 1/r G sink (EX-GC14);
  - falling in a G gradient at the Gc rate, so active and passive mass agree as for plain Gc.

  This holds under our two-term coupling.
- **Universal free fall is not addressed.** Every knot particle here is the same Gc particle. EX-GC01–04 showed that
  different Gc particle types and bound pairs fall at different speeds. That open problem is unchanged.
- **The chemistry fully dresses the knot.** With λ = 1 the gate is 1 everywhere in the core, so the chemistry (M, the
  drift) is exactly that of a plain particle. The knot decides only whether the particle exists and where it is.

## Method
- `src/exgc14_knotgc.py`, with records written to this folder (`EXGC_OUT=EX-GC15`). Settings: h = 0.25, box ±6
  (49³), dt = 0.000375 (scaled by h²), a frame sweep every 20 steps (the same sweep rate per time unit as EX-GC14),
  coupling as EX-GC14 (gate width 1.5, ρ_ref = 0.017).
- The settled h = 0.25 particle is cached by `exgc10_tether.settled`, now with dt scaled by h². Its M is 5.515, between
  the radial 5.32 and h = 0.5's 5.88.
- The drift is a linear fit of the W-centre z position for t ≥ 20.
- Platform: Linux cloud container, Python 3.13.18, numpy 2.5.3, scipy 1.18.1. Compute: about 1.1 core-hours, plus a
  first attempt lost to a container restart.

## Deviations and limits
1. **No token runs.** At Ω = 3000 the core wanders about 0.4 by noise (EX-GC10), against an expected drift of about
   0.03, so a token test would need many replicas. This was stated in the log entry before running.
2. One seed per case, and a single gradient strength and direction.
3. The λ = 0 death was not predicted. The log entry expected the core might detach from the knot; instead the
   shrunken knot could not sustain it.
