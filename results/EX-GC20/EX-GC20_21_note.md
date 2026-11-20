# EX-GC20 / EX-GC21 (EXPLORATORY): wave-like G, and a faster knot

Our construction (Model Gc is ours, not LaViolette's). Pre-run entries: `explore/LOG.md` (EX-GC20, EX-GC21, and
EX-GC21b with the EX-GC20 follow-up runs, logged before they ran). Cloud compute: EX-GC20 1.76 core-hours (cap 2),
EX-GC21 0.79 (cap 1). Platform: Linux cloud container, Python 3.13, numpy 2.x. One run per case.

## Conclusion first
1. **Making G wave-like helps only a little.** A kicked particle still slows down. Its speed-loss rate falls from 0.023
   (diffusive G) to 0.016–0.017 (τ_J = 2 and 20), a 27% reduction. We predicted roughly 1/D_eff, i.e. 10–80 times
   smaller.
2. **The slowing is a real force from the particle's own G inflow.** With the recoil push switched off, the speed stays
   at exactly 0.05, so the flow model itself conserves momentum. The loss therefore comes from the push. As the particle
   moves, it keeps absorbing G that is at rest in the medium, and that pushes it backwards.

   **Hypothesis (not verified here):** this is accretion drag. A consuming particle moving through a medium that has a
   rest frame must keep taking in material at rest, and no change in how G propagates removes that. Only a medium with
   no rest frame would, as with Rueda–Haisch–Puthoff's Lorentz-invariant zero-point field. A testable consequence is
   that the drag should shrink with how strongly the particle depletes its surroundings. Our particles are very
   strong-field (G at the core is about 0.77).
3. **Fall in the wave-like medium.** For the first ~7 time units the particle accelerates exactly at push/mass (0.002).
   After that the drag builds up. K = 0.2 and K = 0.4 average 0.00171 and 0.00152, 12% apart (4% in diffusive EX-GC19).
   M rings by ±15% from waves bouncing off the box walls.
4. **The knot.** Four times more frame sweeps, and also λ = 4, did not help: the knot still lags by about 0.3 after the
   core has moved 0.39. **Carrying the knot with the flow** (shift it one site whenever the core's flow displacement
   reaches h) keeps it within one lattice step of the core, which is the best a lattice can do. That part works.

## Numbers
| run (h = 0.25, τ = 20, drag off, κ = 1) | result |
|---|---|
| kick 0.05, diffusive (τ_J 0.25, D_eff 1) | speed-loss rate 0.0232 (EX-GC19: 0.0219) |
| kick 0.05, τ_J 2 (D_eff 8) | 0.0164 |
| kick 0.05, τ_J 20 (D_eff 80, c = 2) | 0.0170; speed 0.037 at t = 19 |
| same, recoil push off | speed 0.0500 constant (control) |
| fall, τ_J 20, K 0.2 (external flux started at steady value) | mean acceleration 0.00171; M 5.5 → 6.2 ± 15% ringing |
| fall, τ_J 20, K 0.4 | 0.00152 |
| (first fall runs, flux started at 1/80 of steady: 0.00081, 0.00076; superseded) | |
| knot, 4× sweeps, λ 1 / λ 4 | knot z −0.11 / −0.17 when core is at −0.39 (t = 19) |
| knot carried by the flow | knot −0.28 (one shift of 0.25 behind core −0.39 by quantisation) |

Code: `src/exgc20_wave.py` (G flux on cell faces: ∂G/∂t = −∇·J − s, ∂J/∂t = −c²∇G − J/τ_J; recoil push 2J/G).
Figure: `paper/figures/exgc20_wave.png`.

## Deviations
- The first two wave-medium fall runs started the external flux at the diffusive value, 1/80 of its steady value. They
  were repeated with the correct start, logged before the rerun, and both sets are kept.
- The push-off control and the flow-carried knot were added after the first results. Both were logged before running.
- The knot gate was refreshed every 20 steps even when frames swept every 5.

## What it means for the plan
Propagation speed was not the bottleneck for inertia. The remaining obstacle is that G has a rest frame. Two ways
forward:
- **(a) Test the hypothesis cheaply.** Measure coasting time against particle strength (a weaker sink should coast
  longer, roughly as 1/depletion). Ordinary matter would then be in the extremely weak limit.
- **(b) A Lorentz-invariant G field.** This is a deeper model change toward field theory.
