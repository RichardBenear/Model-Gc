# EX-GC06 (EXPLORATORY): a Gc core with its own clock

Claude, 2026-11-18, at the PI's request ("ok, go ahead with that test"). Exploratory track. "Gc3" (Gc1 plus a slow
core inhibitor) is our construction, not LaViolette's.

## Answer
**Yes, the core gets a clock, but it is a breathing clock, not a rotation.** A slow inhibitor V of the core species W
makes the particle pulse in size (l = 0 oscillation) before any other mode moves. Non-round modes behave as follows:
- l = 1 (a wobble of the core about the centre) starts oscillating later, at a much longer delay;
- l = 2 and l = 3 never oscillate in the range scanned.

| Brake κ | Steady consumption M | Breathing onset τ_V | ω at onset (period) | l = 1 oscillation onset |
|---|---|---|---|---|
| 0 (Gc1) | 5.32 | none | — | none |
| 0.5 | 3.26 | none up to 50 | — | none |
| 0.75 | 2.41 | none up to 50 (damped breathing, ω ≈ 0.2–0.4) | — | none |
| 0.9 | 1.89 | ≈ 5.8 | 0.31 (≈ 20) | none up to 50 (−0.02 at 50) |
| 1.0 | 1.51 | ≈ 3.0 | 0.36 (≈ 18) | ≈ 19, ω ≈ 0.2 |
| ≥ 1.1 | no particle | — | — | — |

The breathing frequency falls as τ_V grows. At κ = 1 and τ_V ≥ 20 the breathing
pair turns into a real growing mode, a collapse or blow-up rather than an oscillation, so the useful window is near
onset. With V spreading as fast as X (D_V = 1 instead of 0.1), nothing oscillates: the brake must stay local to the
core.

## What it means
- **An internal clock exists in a Gc-type particle** (item 1 of the earlier list). Its consumption, the particle's
  gravitational mass, would oscillate about a mean. Close to onset that oscillation is small, and far away the
  diffusing G field sees mainly the average.
- **No spin-like rotation yet.** l = 2 (a rotating "peanut") never oscillates. l = 1 does oscillate past
  τ_V ≈ 19 at κ = 1. In 3D an oscillating l = 1 mode can be a back-and-forth wobble or a circling of the core around
  the centre (m = ±1, two senses). The circling version would be a rotating particle with a breathing clock already
  running, which is the ingredient list for the 2:1 test. Which one appears, and whether the two lock, is a nonlinear
  question that needs a 3D time simulation. It also has to be done with breathing already unstable there.
- The clock frequency is set by the brake parameters (κ, τ_V), not by M, so there is no built-in de Broglie-type
  "frequency ∝ mass" relation. Across κ = 0.9 → 1.0 at onset, M falls 20% while ω rises 15%. That is the opposite
  sign, though it mixes two different onset points.

## Method
- `src/exgc06_clock.py`: 5-species radial model (X, Y, W, G, V):
  W_t = D_W∇²W − W + cYW² − W³ − κVW, τ_V V_t = D_V∇²V + W − V; X, Y, G as Gc1 (G feeds X, k = 0.2). V = 0 in the
  ether, so the far field and the G conservation are unchanged. The steady particle was continued in κ from the
  EX-GC01 reference (V = W at the start). Eigenvalues of J_l ψ = λ T ψ, T = diag(1, 1, 1, 1, τ_V), for l = 0–3 by
  shift-invert at 9 shifts (including positive ones). A mode counts as a core mode if ≥ 80% of |ψ|² is in X, Y, W, V
  within r < 10. Steady particles are saved as `particle_*.npz` (small).
- `src/exgc06_summary.py` → `summary.json`, `paper/figures/EX-GC06_clock.png`.
- Platform: Linux cloud container, Python 3.13.18, numpy 2.5.3, scipy 1.18.1. Compute: about 1.5 core-hours.

## Deviations and limits
1. **κ = 2 has no particle**, so κ = 0.75, 0.9, 1.1, 1.25 and 1.5 were added. Particles exist up to κ = 1.0; the
   continuation fails from 1.1 on.
2. **The first pass missed fast-growing real modes**, because its shifts were all on the imaginary axis. It was redone
   with positive shifts, and only the redone results are kept.
3. **μ (drift mobility) of the clocked particle was not computed**, although the log entry asked for it. The existing
   mobility code is written for 4 species.
4. **Linear only.** Whether breathing saturates into a steady pulsation (supercritical) or blows up, and what the
   oscillating l = 1 mode turns into, needs time simulations.
5. τ_V was sampled at 2, 3, 4, 5, 7, 10, 20, 50, so the onsets in the table are linear interpolations.
