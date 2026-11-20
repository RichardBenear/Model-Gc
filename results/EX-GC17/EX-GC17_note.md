# EX-GC17 (EXPLORATORY): drag and acceleration response of the recoil consumers

Our construction (Model Gc is ours, not LaViolette's). Pre-run entry: `explore/LOG.md` (2026-11-20). Prompted by the
PI's question of whether meeting tokens head-on acts like inertia, as in Rueda–Haisch–Puthoff. Cloud compute 1.00
core-hour. Platform: Linux cloud container, 4 cores.

## Conclusion first
1. **Head-on tokens give drag, not inertia.** A particle moving through the token gas feels a force proportional to its
   velocity. Per unit of inertia that force is γ ≈ m₀/6, within 1–2% at m₀ = 50 and 2–8% higher at m₀ = 500. Because
   drag and the gravity push both come from the same absorption, they scale together.
2. **Terminal fall speed is universal.** In an overdamped world (which Gc particles are), a particle falls at push ÷ drag.
   For all five EX-GC16 types this comes out within ±2% of 2E. That holds even at strong consumption, where the push
   alone differs by 2–7% between types. In the mean-field calculation, self-depletion cancels between push and drag to
   ≤1%.
3. **The response to acceleration has the wrong sign for inertia.** With an oscillating flow, the part of the force
   that follows acceleration pushes *along* the acceleration (κ > 0). That is anti-inertia: the depleted wake behind the
   particle lags behind. It is small, κ/γ ≈ 0.01–0.07 time units, and it scales like m₀²/Ω, so it vanishes in the weak
   limit.

So this medium gives Aristotelian free fall (a universal terminal speed), not Newtonian inertia. A coasting particle
would stop. Real inertia must come from somewhere else: either the particle's own internal dynamics (compare EX-GC02),
or a medium with no rest frame, which is the ingredient Rueda–Haisch–Puthoff rely on (a Lorentz-invariant zero-point
spectrum).

## Set-up
The lattice and consumer rule are the same as EX-GC16 (L = 10, Ω = 1000, dt = 0.1), with E = 0.
- The particle is held fixed while the gas drifts past it along z. Hop rates are 1 ± u/2 along ±z. This is equivalent to
  the particle moving at v = −u.
- **Steady flow** measures γ = force / u.
- **Oscillating flow**, u = u₀ sin ωt, is fitted as F = −γ v + κ a. Inertia would show up as κ < 0.

A separate path, `src/exgc17_summary.py`, solves the continuous-time linear mean field. It computes the steady field,
then the first-order response n₁(ω) = (iω − A₀)⁻¹(A₁n₀ + b₁), and linearises the arrival-imbalance force. It shares
only the unit geometry with the simulator.

## Numbers (per unit of inertia; one seed for steady runs, two for oscillating)
Drag γ/(m₀/6) in steady flow (m₀ = 50 at u = 0.1; m₀ = 500 at u = 0.3):
| type | m₀ 50 | m₀ 500 | mean field m₀ 500 |
|---|---|---|---|
| single | 0.987 | 1.038 | 1.042 |
| 2× rate, one site | 1.011 | 1.078 | 1.087 |
| spread over 7 sites | 0.995 | 1.022 | 1.021 |
| two units stacked | 1.024 | 1.035 | 1.037 |
| two units side by side | 1.002 | 1.042 | 1.049 |

The m₀ = 500 single unit at u = 0.05, 0.1, 0.2 and 0.3 gives 1.019, 1.043, 1.034 and 1.038: linear in u within noise.

Terminal fall speed, |EX-GC16 fall| / γ, divided by 2E:
| type | m₀ 50 | m₀ 500 |
|---|---|---|
| single | 1.010 | 0.985 |
| 2× rate, one site | 1.007 | 0.989 |
| spread over 7 sites | 1.009 | 1.009 |
| two units stacked | 0.979 | 1.014 |
| two units side by side | 1.001 | 1.013 |

For comparison, the push alone at m₀ = 500 is 1.023 to 1.067. The mean-field terminal speeds are 1.000 to 1.010.

Acceleration coefficient κ per unit (m₀ = 500, u₀ = 0.5):
| type | ω = 1 measured | mean field | ω = 0.2 measured | mean field |
|---|---|---|---|---|
| single | 0.83 ± 0.17 | 1.02 | 3.6 ± 2.0 | 3.04 |
| 2× rate, one site | 2.14 ± 0.06 | 2.13 | 6.0 ± 0.9 | 6.34 |
| spread over 7 sites | 0.51 ± 0.06 | 0.56 | 1.09 ± 0.08 | 2.34 |
| two units side by side | 1.17 ± 0.08 | 1.15 | 3.05 ± 0.40 | 4.31 |

At m₀ = 50 and ω = 1, κ = 0.05 and 0.02, against predictions of 0.01 and 0.02 (one seed each).

Every measured κ is positive. At ω = 1 they agree with the mean field. At ω = 0.2 the spread unit (and, at 3σ, the pair)
falls below the mean field. A likely reason is that the mean field splits a spread unit's sink equally over its sites,
whereas the rule draws in proportion to arrivals, which shift during the cycle. This is not checked further.

Figure: `paper/figures/exgc17_drag.png`.

## Deviations from the log entry
- The oscillation amplitude was u₀ = 0.5, not 0.1, to beat the noise. Linearity in u was checked for the steady drag.
- Steady runs at u = 0.3 were added for all types (lower noise).
- The κ errors come from only two seeds, so they are rough.
- No oscillating runs were done for the stacked pair.
