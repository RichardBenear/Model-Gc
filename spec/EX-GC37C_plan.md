# EX-GC37C — does the event-counting mechanism survive a moving mass? (preferred-frame test)

2026-11-26 | EXPLORATORY | Claude | PI go-ahead in thread ("close 37A/B. Go ahead with the next steps", 2026-11-26 01:31 UTC). Written before any run. Default cap 2 core-hours; expected under 0.1. This is our construction, not LaViolette's Model G.

## Question

In 37B, a firing-rate deficit set by a carrier field gives the weak-field metric. That was shown for a mass at rest. In GR, a mass moving uniformly carries the same field with it:
- Lorentz-contracted along the motion;
- the same ahead and behind;
- no self-force.

Which kind of carrier field reproduces that?
- **(a)** A **diffusing** carrier, as in Model Gc and 37B B2.
- **(b)** A **wave-like** field with speed 1, sourced by the same consumption: (∂_t² − ∇²)δN = source. This is the rest-frame-free field of EX-GC25.

## Method

- Comoving steady state of a Gaussian sink (core 1.5, total consumption fixed so that GM = 0.16 at rest) moving at v along x.
- Solved in Fourier space on a 128³ periodic lattice with the Y1 symbol K_a = 2 sin(k_a/2):
  - (a) diffusion: (−D K² + i v k_x) ĉ = ŝ, with D = 1;
  - (b) wave: (v² k̃_x² − K²) ĉ = ŝ, where k̃_x = K_x (the lattice time derivative is the same centred difference), for v < 1.
- Velocities: v = 0, 0.05, 0.1, 0.2, 0.4.
- Readouts:
  - **R1, fore/aft ratio:** the deficit 10 cells ahead divided by the deficit 10 cells behind (GR: 1);
  - **R2, side/axial ratio:** the deficit 10 cells to the side divided by the deficit 10 cells ahead, compared with GR's boosted 1/r, which gives γ = 1/√(1 − v²);
  - **R3, self-gradient:** the gradient of δN along x, averaged over the core with the sink profile as weight. A particle pulled toward lower N then feels a self-force F ∝ −∇N (GR, and any rest-frame-free field: 0). This is the drag test.
  - **R4, far-field reach:** the deficit at 40 cells ahead and behind, relative to the v = 0 field.

## Predictions (scored)

- **P1:** the diffusing carrier has a preferred frame. At v = 0.1 (Péclet number vd/D = 1 at d = 10), R1 differs from 1 by more than 10%, and R3 is nonzero and grows linearly with v. Expectation: yes. This is EX-GC25's objection made quantitative.
- **P2:** the wave field has none. For every v ≤ 0.4:
  - R1 = 1 within 1%;
  - R2 = γ within 2% (lattice);
  - |R3| ≤ 1e-3 of the diffusive value at the same v.

  Expectation: yes.

## What it means

- **If P1 and P2 both hold:** the counting mechanism survives motion only if the firing-rate field is a wave field, i.e. relativistic, rather than a consumed, diffusing substance.
- **That turns the open question** "why a coincidence?" into a sharper one: what in the substrate makes N obey a wave equation?
- **One natural candidate, an inference that is not tested here:** if changes in the firing rate themselves spread by two-site coincidences, they travel at the same speed as light.

## Outputs

`src/exgc37c.py`, `results/EX-GC37C/`. The separate check is a v = 0 limit, which must reproduce the 37B B2 fit, and an analytic continuum comparison for the wave case.

**Amendment A1 (before any run).** Y1 and the 35I/35K event hosts have continuous time, so the comoving time derivative of a moving profile is −v∂_x with the exact Fourier symbol ik_x, not the lattice symbol K_x. The wave case is therefore (K² − v²k_x²)ĉ = −ŝ/D. The grid is 160³, not 128³, to keep the 40-cell far-field readout away from the periodic images.
