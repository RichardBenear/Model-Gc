# EX-GC05 (EXPLORATORY), Stage A: is there a rotating Gc1 particle?

Claude, 2026-11-18, at the PI's request ("go ahead with number 1"): route 1 toward a 720° (spinor-like) particle,
meaning a rotating particle whose internal X/Y oscillation is locked 2:1 to its rotation. Stage A asks whether any
non-round mode of the particle can start oscillating, which is what a rotating pattern grows from. Exploratory track;
Model Gc is our construction, not LaViolette's.

## Answer
**No. In Gc1 the particle never starts to rotate.** Every core mode with angular index l = 0–3 stays real (no
oscillation) and decaying, across all three knobs tried:

- **Y response time τ_Y (the Brusselator oscillation knob), from 1 to 8.5.** Core modes stay real. The first
  oscillation to appear is the ether itself: at τ_Y = 8.64 the whole background X/Y medium starts oscillating, with
  angular frequency ω ≈ 5.3. By τ_Y = 9 every l grows at about +0.3 to +0.5 with ω ≈ 5, which is a bulk, not a
  particle, instability. The particle offers no window in which it oscillates while the ether is still steady.
- **G feeds Y instead of X (k = 0.2 and 0.4).** Feeding Y in the core should, in principle, make the core more
  oscillation-prone. It does not change the picture: all core modes stay real up to τ_Y = 8.4.
- **Slow core W (τ_W from 10 to 40, through the EX-GC02 drift threshold at 28.3).** The non-round core modes slow down
  exactly as 1/τ_W (τ_W × rate is constant at about 0.80 for l = 2), but they stay real and never cross zero.

**Why (established in part).** The Brusselator's local oscillation threshold does fall inside the particle: at
r ≈ 1 the local X/Y chemistry would oscillate at τ_Y > 4.8, compared with 8.64 in the ether. But that region is a thin
shell, and diffusion couples it to stable surroundings. The only trace of it is a near-collision of two real l = 2
eigenvalues at τ_Y ≈ 4 (−0.2469 and −0.2480), which separate again without becoming an oscillating pair. The core
species W obeys a single cubic equation, so it cannot oscillate by itself.

## What it means for the 720° route
A 2:1 lock needs an oscillating or rotating particle first, and Gc1 has neither. The X/Y chemistry oscillates only
as a whole medium. Route 1 would need a new ingredient in the core: for example, a second core species acting as a
slow inhibitor of W, so that the core has its own clock. In other three-species reaction–diffusion systems this gives
breathing (l = 0) and rotating or oscillating non-round spots (Bode et al. 2002). Routes 2 (an orientation field) and
3 (ribbon or braid structures in the hypergraph) need a new ingredient anyway. So this result does not favour route 1
over them.

A related idea is open: a self-propelled particle (beyond the EX-GC02 drift threshold) can turn onto a circular path
when its l = 2 deformation relaxes slowly enough (Ohta & Ohkuma, PRL 102, 154101, 2009). Here l = 2 relaxes at
0.80/τ_W, which is slow, so circular motion is possible. Finding it needs a nonlinear (3D time-dependent) calculation,
not this linear scan. That would be orbital motion of the whole particle, not internal rotation.

## Method
- `src/exgc05_rotation.py`: leading eigenvalues of J_l ψ = λ T ψ (T = diagonal response times; the steady particle is
  unchanged), found by shift-invert at shifts 0, 0.5i, …, 7i, 20 eigenvalues per shift. Each mode is tagged by its
  core share, the fraction of |ψ|² in X, Y and W within r < 10. The slow G-diffusion ladder and the ether modes have
  small shares. `src/exgc05_tauW.py`: the same for τ_W. Summary and figure: `src/exgc05_summary.py` → `summary.json`,
  `paper/figures/EX-GC05_rotation_search.png`.
- Particle: the EX-GC01 'ref' radial particle (G feeds X, k = 0.2, R = 120, dr = 0.05). The G-feeds-Y particles
  were continued from it with the EX-GC03 builder.
- Platform: Linux cloud container, Python 3.13.18, numpy 2.5.3, scipy 1.18.1. Compute: about 1.2 core-hours. A
  container restart lost one τ_W run, which was repeated.

## Deviations from the log entry
1. **B was not varied.** In the ether, B and τ_Y move the Hopf point the same way (threshold τ_Y = A²/(B − 1)), and τ_Y
   was scanned through it.
2. **Two knobs were added:** the G-feeds-Y coupling and a slow core W. Both test the same Stage A question.
3. **Core-share tagging is approximate.** Some ladder modes with share ≈ 0.5–0.6 appear among the "core" modes in
   `summary.json`, for example −0.0206 for l = 1. Their real eigenvalues do not change the conclusion, because no
   mode of any share is oscillatory before the ether Hopf.
4. **Stage B was not started**, as the log entry specified, because Stage A found no oscillatory l ≥ 1 instability.
