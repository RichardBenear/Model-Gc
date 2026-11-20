# EX-GC19 (EXPLORATORY): a flow (momentum) field for the Gc particle

Our construction (Model Gc is ours, not LaViolette's). Pre-run entry: `explore/LOG.md` (2026-11-20). Cloud compute
1.81 core-hours (cap 2). Platform: Linux cloud container, Python 3.13, numpy 2.x. One run per case.

## Conclusion first
1. **With the flow field, the particle moves as one body.** The whole pattern (X, Y, W) rides a shared velocity, is
   accelerated by the recoil push, and survives (M settles about 8% lower). This is what failed in EX-GC18.
2. **It accelerates rather than drifting.** With drag off it falls with steadily growing speed, about 0.002 per time
   unit, close to the push/mass prediction of 0.0025. With drag on it approaches a terminal speed.
3. **Acceleration is nearly universal.** The K = 0.2 and K = 0.4 particles accelerate within 4% of each other (0.00196
   against 0.00189 averaged over t ≤ 24; 0.0022 against 0.0023 early on). The rigid-shift speeds in EX-GC18 differed by
   10%.
4. **It does not coast forever, even with drag switched off.** After a kick, the speed decays with a time constant of
   about 46. The particle's own G deficit lags behind it (G diffuses), and the push toward that lagging deficit acts as
   drag. This wake drag is intrinsic to a diffusing medium, which has a rest frame. **True inertia would need G to
   propagate as waves rather than diffuse.**
5. **Cohesion is needed.** Without it ("zero viscosity"), the inward self-push crushes the core (M −70% by t = 24), as
   predicted.
6. **The knot cannot keep up.** The core reaches a speed of about 0.047, but the knot moves only 0.13 against the
   core's 0.61. The frame dynamics has a speed limit well below this. The particle is still alive at t = 24, since the
   gate is wide, but it would lose its knot.

## Model
These are additions to the deterministic Gc1 chemistry at h = 0.25, box ±6:
- X, Y and W deviations from background are carried by the flow: −∇·((f − f₀)V). G is the medium and is not carried.
- Inside the core (weight s/(s + ε), with s = KGW), mass density ∝ s:
  `dV/dt = [−2 D_G ∇ln G − drag·V]/τ`, with τ = 20 (mass/drag).
- Cohesion: −κ(V − ⟨V⟩), where ⟨V⟩ is the s-weighted Gaussian average over width 1.5, and κ = 1 (or 0).
- V is updated every 20 chemistry steps.

Code: `src/exgc19_flow.py`. Figure: `paper/figures/exgc19_flow.png`.

## Numbers (t = 24)
| run | core z | core speed | M/M₀ | comment |
|---|---|---|---|---|
| (a) fall, drag on | −0.44 | 0.0289 | 0.92 | model with total decay 0.073 predicts 0.0278 |
| (b) fall, drag off | −0.61 | 0.0471 | 0.91 | ½at² with a = 0.0021 gives 0.60 |
| (c) kick 0.05, drag off, E = 0 | −0.93 | 0.0296 | 0.92 | decay rate 0.0219 (time 46): wake drag |
| (d) kick 0.05, drag on, E = 0 | −0.57 | 0.0087 | 0.93 | decay 0.0729 = bare 0.05 + wake 0.023 |
| (e) fall, drag off, κ = 0 | −0.30 | 0.0015 | 0.30 | self-compression |
| (f) fall, drag off, K = 0.4 | −0.63 | 0.0454 | 0.86 | 4% below (b) |
| (g) knot, drag off | −0.61 (knot −0.13) | 0.0471 | 0.91 | knot lags by 0.48 |

The wake-drag rate from (c), 0.022, adds to the bare drag 1/τ = 0.05 and accounts for (d) (0.073) and (a) within 4%.
That is the internal consistency check.

## Deviations
- κ ("cohesion") stands in for viscosity. The log entry explains why: Laplacian viscosity cannot resist uniform
  compression.
- T = 25 instead of 30, to stay within the 2 core-hour cap.
- Only one run per case.

## What next
- **Make G propagate as a wave** (a telegraph or wave equation with a finite speed c). A moving particle's G field
  would then travel with it, and the wake drag should vanish for uniform motion, which is what inertia requires. This
  is also the "diffusive propagation" gap noted in docs/PHYSICS_MAPPING.md.
- Speed up the knot's frame dynamics (more sweeps per time unit, or a pull proportional to the gradient of W), or
  accept that the knot moves only slowly.
- An orbit test (two particles, sideways kick) makes sense once coasting works.

## Erratum (2026-11-20, after EX-GC22)
EX-GC22 (ChatGPT; reviewed in `review/EX-GC22_claude_review.md`) showed that most of the "coasting loss" reported here
was a defect of this flow model, not a physical drag. V was not carried with the moving particle (no (V·∇)V term). The
particle's own inward push builds a small radial V in the core, and as the pattern moves through it, the
consumption-weighted core velocity falls even with zero net force. The "own-wake drag" interpretation (EX-GC19) and the
"accretion drag" hypothesis (EX-GC20) are not supported. The push-off control could not distinguish the two
explanations, because it also removes the radial V. The numbers above are unchanged as records.
