# EX-GC16 (EXPLORATORY): recoil rule for universal free fall

Our construction (Model Gc is ours, not LaViolette's). Pre-run entry: `explore/LOG.md` (2026-11-20). Cloud compute
1.20 core-hours (cap 2). Platform: Linux cloud container, 4 cores, Python 3 / numpy 2.x.

## Conclusion first
1. **Free fall is universal in the weak regime.** Five different particle types (one unit; double rate at one site;
   one unit spread over 7 sites; two units stacked; two units side by side) all fall at m₀E/3 per unit of inertia,
   within 0.2–0.5% (m₀ = 50, local G depletion about 1%).
2. **Two consumers attract with an inverse-square force.** With no gradient, two units pull toward each other, with
   equal and opposite forces. The measured pull matches the lattice mean-field prediction within 1–3% at d = 2–4 and
   within 9–11% (1σ) at d = 6–8. That prediction approaches the continuum m₀²/(12πΩd²) at large d. The "active mass"
   and the "passive mass" are the same number, m₀.
3. **Strong consumption breaks universality by a few percent.** At m₀ = 500 the local G dip is about 13%. The particles
   then fall 2–7% faster than m₀E/3, and the more concentrated the consumption, the faster they fall (double rate at
   one site falls fastest, the single unit slowest). The mean-field calculation predicts this ordering. This works
   like a self-energy (Nordtvedt-type) violation of the equivalence principle. Its size scales with the local depletion
   m₀/Ω, so it vanishes when the token background is dense.

## The rule
Here is how the test lattice is set up:
- G tokens random-walk on a cubic lattice: h = 1, D = 1, dt = 0.1, Ω = 1000 tokens per site at G = 1.
- A reservoir at the faces holds G = 1 + E z.
- Each consumer unit absorbs G at a fixed rate m₀, zero order in G. Each step it takes Poisson(m₀ dt) tokens from those
  arriving at its site(s) that step.
- Each absorbed token pushes the particle one unit along the token's direction of travel.
- Particles are held fixed, and the push per unit time is recorded. Fall rate = force / B.

The push is measured two ways: as the realised sum, and as its expectation given that step's arrival counts. The
expectation has lower noise, and both agree.

Why it works: tokens arriving from the low-G side are scarcer, so the net push points toward lower G. That push equals
m₀ × (fractional imbalance of arrivals), which is m₀ × 2E/6 at G = 1. The force is proportional to m₀, and so is the
"inertia" (number of units, each with the same m₀). Their ratio does not depend on the particle's makeup.

## Numbers (fall_z; ideal −m₀E/3; 2 seeds × T = 1950 each)
| type | m₀ = 50 | m₀ = 500 | mean-field m₀ = 500 |
|---|---|---|---|
| ideal | −0.3333 | −3.333 | |
| single | −0.3323 ± 0.0035 | −3.408 ± 0.032 | −3.472 |
| 2× rate, one site | −0.3392 ± 0.0191 | −3.556 ± 0.024 | −3.622 |
| spread over 7 sites | −0.3346 ± 0.0018 | −3.438 ± 0.046 | −3.435 |
| two units stacked (d 4) | −0.3342 ± 0.0036 | −3.499 ± 0.010 | −3.482 |
| two units side by side (d 4) | −0.3343 ± 0.0022 | −3.516 ± 0.019 | −3.497 |

E = 0 controls (single seeds, m₀ = 500): the single unit shows only noise (|F| ≤ 0.2). The stacked pair pulls together
with +0.42 against a prediction of +0.50, from one seed with noise of about 0.1.

Pair pull without a gradient (m₀ = 500, L = 12, 2 seeds):
| d | measured | lattice mean-field | ratio |
|---|---|---|---|
| 2 | 2.553 ± 0.066 | 2.592 | 0.985 |
| 3 | 1.015 ± 0.106 | 0.993 | 1.022 |
| 4 | 0.488 ± 0.002 | 0.501 | 0.973 |
| 6 | 0.223 ± 0.028 | 0.206 | 1.085 |
| 8 | 0.129 ± 0.015 | 0.116 | 1.11 |

In the mean-field calculation, the "mutual" part (with the box-wall effect removed) times d² falls from 10.4 at d = 2 to
6.8 at d = 8, approaching the continuum 6.63. The steeper fall-off at small d is a lattice effect.

Separate path: `src/exgc16_summary.py` solves the steady lattice diffusion equation with the sinks and the face
values, and computes the expected arrival imbalance. It shares only the unit geometry with the simulator.
Figure: `paper/figures/exgc16_recoil.png`.

## Deviations from the log entry
- The bound pairs used separation d = 4, not 3.
- The pair test used d = 2, 3, 4, 6, 8 on a larger box (L = 12, 25³), not d = 3–8 on 21³. This kept the box walls further
  away. For d = 6 and 8 we used T = 6000.
- Only one seed was used for each E = 0 control.

## Limits and what this does not show
- Inertia = B is assumed, not emergent. Particles are held fixed, so this measures force per unit of assumed inertia.
- **Drag is untested.** A Le Sage-type push predicts that a moving absorber meets more tokens head-on and feels a drag
  proportional to velocity. That is the classic objection to Le Sage gravity, and it is the next thing to test.
- The consumer units are bare rules, not yet the EX-GC14/15 knot particle. Joining the two (a knot whose core absorbs
  at a fixed rate and recoils) is the step after that.
