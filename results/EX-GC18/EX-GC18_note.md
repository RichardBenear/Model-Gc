# EX-GC18 / EX-GC18b (EXPLORATORY): recoil rule on the knot particle

Our construction (Model Gc is ours, not LaViolette's). Pre-run entries: `explore/LOG.md` (EX-GC18, then EX-GC18b
logged after the first runs). Cloud compute 1.76 core-hours (cap 2). Platform: Linux cloud container, Python 3.13,
numpy 2.x.

## Conclusion first
**This step failed.** Neither way of applying recoil makes the knot particle move as one object:
- **Kicking the knot alone tears the knot off the core.** The chemistry core stays where it is, held in place by its own
  G deficit. The knot gets about 0.8–1.0 away before the pull toward W balances the kicks. The stronger particle
  (K = 0.4) starts to die once its knot has moved off.
- **Pushing the absorber W locally does not move the particle either.** The push points toward lower G, and inside
  the core that means toward the centre, so the core squeezes itself (M −42% and −65%). The chemistry then simply
  re-anchors W: the core sits displaced by about 0.1 and does not travel.
- **Only a rigid shift of the whole particle moves it** at the recoil speed (about 0.048 at K = 0.2, 0.053 at K = 0.4).
  But that imposes rigidity by hand.

Even then, the two particle types differ by about 10%. Here, unlike on the token lattice (EX-GC17), the wake does not
compensate. The cause is that the core depletes G deeply (G at the core is about 0.77), and ∇ln G is steeper there for
the stronger consumer.

**What this teaches us.** A Gc particle is a pattern held together by slow X/Y/G chemistry. A push on any one part, the
knot or the absorber, is absorbed by that pattern rather than moving it. For recoil to move the whole particle, the
pattern needs a way to carry momentum as a whole. Overdamped reaction–diffusion has no such mechanism. This is the same
gap as the missing inertia (EX-GC17), so steps 1 and 2 of the plan are really one problem: **the Gc fields need a
momentum (or velocity) variable.**

## Runs (h = 0.25, box ±6, E = 0.02, T = 30, one seed; recoil v = −2 D_G ∇ln G, D_G = 1)
| run | core z at t = 29 | knot z | M(29)/M(0) | mean recoil v_z |
|---|---|---|---|---|
| knot kicked, K 0.2 | −0.04 | −0.82 | 0.983 | −0.0491 |
| knot kicked, K 0.4 | −0.23 (dying) | −0.97 | 0.759 | −0.0540 |
| knot kicked, K 0.2, E = 0 (control) | 0.00 | 0.00 | 1.000 | 0 |
| rigid shift, K 0.2 (no knot) | −1.26 (5 shifts) | — | 0.979 | −0.0483 |
| rigid shift, K 0.4 (no knot) | −1.52 (6 shifts) | — | 0.977 | −0.0535 |
| W pushed (field), knot K 0.2 | −0.05 | −0.04 | 0.585 | −0.046 |
| W pushed (field), knot K 0.4 | −0.11 | −0.13 | 0.350 | −0.048 |
| W pushed (field), plain K 0.2 / 0.4 | identical to the knot runs | | | |

The knot's degree stayed at 1 in every run. The K = 0.4 particle was settled freshly on the lattice (M = 8.17). For
comparison, the chemical drift without recoil is −5 × 10⁻⁴ (EX-GC15), about 100 times smaller than the recoil speed.

Code: `src/exgc18_knotrecoil.py`. Figure: `paper/figures/exgc18_knotrecoil.png`. No separate recomputation path was
needed: the failures are visible directly in the trajectories.

## Deviations
- EX-GC18b (field-level recoil) was a second design toward the same approved goal. It was logged after EX-GC18's first
  results and before its own runs.
- The rigid-shift K = 0.4 run was added after the first batch.
- Only one seed per run, and a single gradient.

## Possible next steps
- Give the Gc fields a velocity or momentum variable (for example a telegraph-type, second-order-in-time core or an
  advecting flow field). That would test whether recoil can then move the whole soliton, and the same change would
  supply inertia.
- Alternatively, accept rigid composites as a modelling choice: identity and position carried by the knot, with the
  chemistry re-centred on it. This is simpler but imposed.
