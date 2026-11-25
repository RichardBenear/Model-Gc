# EX-GC35B — GC26 Skyrmions as weak sources on Y1 (EXPLORATORY)

2026-11-24 | ChatGPT | [Registered plan](../../spec/EX-GC35B_plan.md) · [Direction review](../../review/EX-GC34Q-Z_claude_review.md). **PI CLOSED 2026-11-24.**

**Registered topology, full-stress and feedback gates failed.** A finer-grid post hoc run recovered both knots and measured their leading response in the energy-only control. The full energy-plus-stress source remained incompatible with the static Y1 equations at the registered tolerance.

## Registered run and seed diagnostic

The implementation used a Gaussian radial B1 start in the first run, departing from the archived GC26 radial-profile prescription. The subsequent seed check used the archived GC26A sample values with linear interpolation and an inverse-square tail, and repeated the finer registered grid. In both runs B1 relaxed to vacuum at h=.3. Both initial coarse cases also relaxed to vacuum. The B2 case at h=.3 retained a localized configuration (oriented density integral −1.660, E=2.2341). Its full source residual was 4.45%, while its energy-only tidal error was 0.122%. The B1 vacuum output has no particle mass or binding interpretation. The archived seed check rules out the choice of Gaussian alone as the cause of B1 loss at the registered finer spacing.

[Registered output](results.json) · [Seed check](corrected_seed.json)

## Post hoc h=.2 recovery

Both configurations were regenerated on 61³ points with fixed vacuum on the outer faces. The accelerated objective was checked against the original autodifferentiated energy and gradient on a small grid. Both optimizations reached the 350-iteration cap, so the table characterizes approximate rest profiles. Degree below is the oriented baryon-density integral with the GC26 sign convention.

| Quantity | B1 | B2 |
|---|---:|---:|
| Degree integral | −0.930118 | −1.878636 |
| Integrated rest energy | 1.203673 | 2.312754 |
| Maximum optimizer gradient component | 2.12e−6 | 2.48e−7 |
| Full static forced residual | 2.09% | 1.55% |
| Energy-only static forced residual | 4.55e−16 | 4.35e−16 |
| Energy-only low-wavenumber tidal error vs matched continuum | 0.0563% | 0.0535% |
| Leading tangent acceleration norm | 0.1790 | 0.2897 |
| Relative change on halving the centered-difference strength | 5.55e−10 | 4.66e−10 |

The static shift source is zero by construction for rest profiles. The metric is the minimum-norm static solution for the nonzero Fourier modes retained in the benchmark. The full source includes the discrete spatial stress from metric differentiation of the Skyrme energy. Its residual measures the component outside the range of the static Y1 kernel and fails the 0.1% compatibility gate. The residual combines mapping, discretization and finite-box effects. This experiment does not isolate those contributions.

The B2 energy lies 3.93% below two copies of B1. The energy-only monopole ratio is 1.921414, equal to the energy ratio at roundoff after dividing by the measured form factor. This is a normalization control on energy sourcing, including the binding deficit in these profiles. It is not an independent measurement of active mass. Feedback uses the energy-only metric, checks the change in initial tangent acceleration at rest, and includes the variation of the kinetic inertia matrix. The peak imposed metric perturbations in the centered check were 5.46e−6 and 7.91e−6. These tests give no full-source reciprocal, translation, free-fall or binding-equivalence measurement.

The combined CPU time recorded across the three runs was about 0.156 core-hours. Reproduce with `python src/exgc35b.py`, `python src/exgc35b_diagnostic.py --spacing .3`, and `python src/exgc35b_refinement.py --spacing .2 --post-hoc`. The last uses Numba for the gradient. Python 3.13.2, NumPy 2.2.4, PyTorch 2.7.1 CPU. [Refinement output](post_hoc.json)
