# EX-GC04 (EXPLORATORY): does a bound pair of identical Gc particles fall like one particle?

Claude, 2026-11-18, at the PI's request ("Go ahead with the bound pair test"). Exploratory track: no claims on its
own. Model Gc is our construction, not LaViolette's.

## Answer
**No. A bound pair falls about 6% more slowly than a single particle.** Binding breaks universal free fall by an
amount set by how much each partner lowers the other's G, and the effect fades as 1/d at large separation.

- Two identical P1 particles (G feeds X, k = 0.2) have **two stable bound separations**: d* = 7.145 and d ≈ 10.11.
  The short-range Turing force pushes them apart inside each separation and pulls them together just outside it.
- In a weak external G gradient, the pair at d* = 7.145 falls at **0.939** of the single-particle speed. The same
  ratio, 0.9395, comes out on a coarser grid (h = 0.125 instead of 0.1). At d = 10.11 the ratio is 0.972.
- The pair consumes **0.954** of 2M at d* (0.967 at 10.11): each particle sits in its partner's G well and eats less.
  This "mass defect" is about 1.5 times smaller than the fall-speed deficit, so the pair is not simply a lighter
  particle that falls normally.
- **Mechanism (checked).** Consumption grows with the local G: dM/dG∞ = 4.58 (single particle at G∞ = 0.941). In a
  gradient, the upper partner sits in higher G, consumes more, and so pulls the lower partner up harder than the
  lower one pulls it down. The pair is dragged back toward higher G. To leading order this subtracts
  (dM/dG)/(8πd) = 0.182/d from the fall ratio. A fit of a/d + b/d² to the d = 14 and 20 points gives a = 0.177,
  within 3% of 0.182. At d = 7–10 the near-field force adds more.
- The ruled-out alternative: lowering the background alone does not slow the fall. A single particle at
  G∞ = 0.941 (the G its partner leaves it in at d*) has μ = −0.02565, the same as at G∞ = 1 (−0.02561).

## What it means
In Newtonian gravity a bound system falls exactly like its parts. That is the weak equivalence principle, and the
Nordtvedt effect would be its failure for self-gravitating bodies. In Gc, the same thing that makes a particle fall
breaks it: consumption that depends on the local G. A particle drifts because it consumes more on its high-G side.
The upper member of a pair consumes more for the same reason, and that pulls the pair back up. So within Gc1 this is
not a tuning problem. It is built into the pairing of "consumption ∝ G" with "drift from asymmetric consumption".

For comparison, lunar laser ranging gives a Nordtvedt parameter η (fall-speed change per unit binding fraction) of
order 1e-3 or less (Williams, Turyshev & Boggs 2004). If the consumption deficit (4.6%) plays the role of the binding
fraction, the fall deficit at d* (6.1%) gives η ≈ 1.3. As the model stands, a composite of Gc particles would violate
the equivalence principle by about its own binding fraction. This is speculative as physics; the numbers are
established for this model.

## Numbers (all from results/EX-GC04/summary.json)
| Case | Fall speed / single | Consumption / 2M | Pair speed at E = 0 |
|---|---|---|---|
| Bound, d = 7.145 (h = 0.1) | 0.9394 | 0.9544 | 2.3e-7 |
| Bound, d = 7.145 (h = 0.125) | 0.9395 | 0.9546 | — |
| d = 7.13, 7.19 (either side of d*) | 0.9398, 0.9381 | 0.954, 0.955 | ∓3.8e-5, +1.1e-4 |
| Bound, d = 10.1 | 0.9717 | 0.9671 | 6.6e-6 |
| Unbound, d = 14 | 0.9840 | 0.9767 | (still attracting) |
| Unbound, d = 20 | 0.9895 | 0.9841 | (still attracting) |

Single particle on the same solver: μ = −0.02561 (h = 0.1) and −0.02565 (h = 0.125), against −0.02554 from the radial
l = 1 calculation (0.3% and 0.4%). A uniform gradient does not stretch the identical pair: the relative-speed
response is below 1e-7 of the mean.

## Method
- Solver `AxiN` (`src/exgc04_boundpair.py`): full axisymmetric (ρ, z) domain, N particles on the axis, one co-moving
  speed and one phase condition per particle (Voronoi tanh switches), block elimination on one sparse LU. The outer G
  boundary is G∞ + E z with E = ±1e-4, and μ_pair is the central difference of the mean speed.
- Pair speed v(d): AxiN at E = 0, from the radial profile guess and by continuation in d (`src/exgc04_continue.py`).
- Summary and figure: `src/exgc04_summary.py` → `summary.json`, `paper/figures/EX-GC04_boundpair.png`.
- Platform: Linux cloud container, Python 3.13.18, numpy 2.5.3, scipy 1.18.1. Compute: about 1 core-hour, inside the
  2 core-hour cap.

## Deviations from the log entry, and limits
1. **The v(d) scan used the full-domain solver, not the mirror solver.** The mirror solver did not converge at d = 8
   from the profile guess. Its two converged points (d = 7.0 and 7.5) agree in sign with AxiN.
2. **The separation was not iterated until both speeds agreed.** For identical particles at a bound separation, the
   mean speed is antisymmetric in E, and its change with d is second order in E. The mean speed at d = 7.13, 7.145 and
   7.19 (0.938–0.940) confirms the result does not depend on d to that precision.
3. **No steady pair state was found for 7.62 < d < 8.25.** Continuation from both sides fails there, and one d = 8.0
   solve from the profile guess converged to an outlier (v = 6.9e-3, ten times its neighbours), which is excluded.
   On both sides of the window the pair attracts, so the window does not change the bound states. Whether the pair
   is time-dependent there was not checked.
4. **Stability is along the axis only.** The bound states are stable to changes in separation. Sideways and shape
   modes, and the full 3D (non-axisymmetric) problem, were not tested. Gradients are along the pair axis only.
5. The far-field estimate a = (dM/dG)/(8π) treats each partner as a point source with the far-field μ. It is a check
   on the mechanism, not a fitted law.

## Possible next steps
- A Gc variant whose drift does not come from G-dependent consumption (so dM/dG and μ come apart), and whether
  composites then fall universally.
- The perpendicular gradient (needs 3D or a different symmetry) and sideways stability of the bound pair.
