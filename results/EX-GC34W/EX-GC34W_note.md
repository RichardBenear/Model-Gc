# EX-GC34W — cubic Y1 vacuum action compatibility (EXPLORATORY)

2026-11-24 | ChatGPT | PI go-ahead: “go ahead with the experiment”. [Plan](../../spec/EX-GC34W_plan.md) and the log were registered before the run (78ccbae and 1dd2dc2). PI CLOSED 2026-11-24 (Richard: “go ahead and close V and W”).

**W0 PASS, W1 FAIL.** Varying the cubic action with respect to Y1's ten metric variables produces a symmetric source. The tested action still violates the linear gauge compatibility condition on free vacuum waves at finite spacing. This occurs with the scalar field absent, so adding scalar matter cannot fix this candidate's vacuum defect. W2 therefore stopped before mixed-matter tests and reciprocal evolution.

## Candidate and method

The gravity candidate is the cubic part of the Einstein–Hilbert connection-product expression after dropping the continuum boundary term, discretized with Y1's staggered derivative and vertex averaging of individual factors. The cubic terms come from the first perturbation of the inverse metric times volume and the second perturbation of the connection. The existing quadratic Y1 action is kept fixed. The scalar coupling was held fixed and is inactive in the vacuum channels.

The two input tensor waves satisfy the full free equations to within 2.4e-16 relative. The source coefficients are defined by inserting each of the ten independent metric variations in the third leg. On the free input equations, the nonlinear variation of the quadratic action vanishes, so the cubic coefficient must vanish when the third leg is a linear coordinate transformation. All four gauge variations were checked. An independent polarization of the polynomial evaluated directly with a gauge leg agrees to within 7.5e-15 in the same normalization.

| Grid spacing | Max normalized contraction, repeated wave | Max normalized contraction, nonparallel waves |
|---:|---:|---:|
| Continuum | ≤ 2.14e-16 | ≤ 2.14e-16 |
| 0.5 | 0.041178 | 0.003936 |
| 0.25 | 0.009958 | 0.000995 |
| 0.125 | 0.002468 | 0.000249 |

The largest repeated-wave failure is the time re-slicing variation for two equal axis waves. For the largest nonparallel failure, the violated conditions are the transverse coordinate displacements in the two directions with nonzero diagonal input components. Both defects decrease about fourfold per halving, but fail the exact gate of 1e-10. No source was projected. The pure-vacuum failure narrows the design problem to the cubic gravity action before any question about the nonlinear scalar coupling arises.

## Implication and scope

This rules out this vertex-averaged cubic action. Endpoint averages supply an exact two-factor product identity on a link, but averaging every factor to a vertex in this action does not preserve the coupled staggered placements required for the cubic gauge identity. A possible follow-up is to write the discrete cubic gauge relation as linear constraints on short-range action coefficients, enforce the desired continuum coefficient, and test whether the constraints have a solution before attempting matter or time evolution. This calculation does not rule out a different local stencil or a different construction.

Run `python src/exgc34w.py`. Python 3.13.2 and NumPy 2.2.4 on Windows, 0.125 CPU seconds. The separate contraction uses the same action factors and free wave solutions, with a different extraction of the three-leg coefficient. No registered criterion changed.

[Raw data](results.json) · [Runner](../../src/exgc34w.py)
