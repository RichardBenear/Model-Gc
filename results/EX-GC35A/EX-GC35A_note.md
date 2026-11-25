# EX-GC35A — finite-time reciprocal scalar feedback on approximate Y1 (EXPLORATORY)

2026-11-24 | ChatGPT | [Registered plan](../../spec/EX-GC35A_plan.md) · [Direction review](../../review/EX-GC34Q-Z_claude_review.md). **PI CLOSED 2026-11-24.**

**A0 PASS, A1 PASS, A2 PASS.** The scalar packet and Y1 geometry were evolved through the first order at which the perturbed scalar feeds back into the geometry. The finite-time test used a standing massive scalar packet in a periodic box, synchronous gauge and perturbative strength 10⁻³. It used the inherited flat scalar sources and Y1 kernel, the scalar response operator from 34R, and a symmetrized position-referred, densitized quadratic effective source. The scalar's leading field was propagated analytically on the flat background. RK4 evolved its correction, both metric orders and the compensator corrections. All residuals use the norm of the original constraint source at the sampled time. The next-order columns include the strength factor 10⁻³.

| Run | Grid | Δt | Peak linear residual | Next residual at t=1 | At t=2 | At t=4 | Peak next residual |
|---|---:|---:|---:|---:|---:|---:|---:|
| Coarse | 24³ | 0.05 | 4.83e-7 | 7.99e-6 | 2.74e-5 | 5.70e-4 | 5.70e-4 |
| Fine | 48³ | 0.05 | 5.13e-7 | 2.08e-6 | 7.64e-6 | 1.72e-4 | 1.72e-4 |
| Half timestep | 24³ | 0.025 | 3.02e-8 | 7.99e-6 | 2.74e-5 | 5.70e-4 | 5.70e-4 |
| Continuum spectrum | 24³ | 0.05 | 5.24e-7 | 5.93e-10 | 1.43e-9 | 4.79e-9 | 4.79e-9 |

The continuum peak next-order residual before strength scaling was 4.79e-6. Halving the spacing reduced the final lattice residual by 3.310×, equivalent to an observed two-point order of 1.73 in this run. Halving the coarse timestep changed that residual by 1.09e-6 fractionally. The fine residual grows within the finite window as initialized constraints respond to the evolving source. These data give the size and time history over t≤4 and do not determine whether much longer evolutions develop secular drift.

Next-order curvature errors relative to the matched continuum run were 5.16% coarse and 1.35% fine at t=2, and 19.1% coarse and 5.27% fine at t=4, using the inherited tensor norm. The change due to scalar feedback alone had nonzero fine-grid tidal norms 0.08123 and 0.1350 at those times, before strength scaling. That contribution is the difference from holding the scalar correction out of the gravity source, by linear superposition at this order. Directly re-sourcing the linear kernel with only the scalar correction gave 2.67e-3 for the final fine constraint residual, versus 1.72e-4 for the approximate completed source. Peak first-order metric components times the perturbative strength were 3.07e-3 coarse and 3.21e-3 fine.

## Periodic energy and flux

An auxiliary stress cancels the uniform leading matter source. Its energy and current corrections obey the expanded conservation equation with zero initial corrections and prescribed zero spatial stress correction. Modes with zero spatial wavenumber are excluded from the metric solve. The mean next-order effective energy at t=4 was 5.43e-3 on the fine grid. This remaining homogeneous source is retained in the data and lies outside the nonzero-mode check.

The table integrates the mean-subtracted effective energy inside a box-centered cube and the corresponding outward flux, found by summing the spectral divergence within the cube. The packet is initially offset from the box center. The cube uses a strict interior mask, so the coarse and fine grids sample slightly different effective volumes. The coarse continuum and coarse lattice use the same mask. Values are signed and are coefficients before perturbative strength scaling.

| Run | Enclosed energy t=0 | Enclosed energy t=4 | Time-integrated outgoing flux | Δ energy + flux |
|---|---:|---:|---:|---:|
| Coarse | -1.5639 | -3.7111 | 2.2047 | 5.75e-2 |
| Fine | -1.7804 | -4.3297 | 2.5628 | 1.35e-2 |
| Half timestep | -1.5639 | -3.7111 | 2.2045 | 5.73e-2 |
| Continuum spectrum | -1.5598 | -3.7505 | 2.1909 | 2.64e-4 |

This ledger is a finite-box analogue of the energy and flux check in 34U. It measures a gauge- and prescription-dependent effective source in a compensated periodic setting. An isolated ADM mass and an evolving homogeneous metric require a different boundary-value benchmark.

Recorded CPU time 0.279 core-hours for the production run, within the two-core-hour cap. Reproduce with `python src/exgc35a.py` using Python 3.13.2 and NumPy 2.2.4. [Time series and controls](results.json)
