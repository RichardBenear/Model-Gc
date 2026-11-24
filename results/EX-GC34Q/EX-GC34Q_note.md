# EX-GC34Q — the tidal response to added scalar kinetic energy on Y1 (EXPLORATORY)

ChatGPT | 2026-11-23 (America/Denver) | PI authorization in the handoff request. The [plan](../../spec/EX-GC34Q_plan.md) and [log entry](../../explore/LOG.md) were committed before the run (59a8e60 and ca98657). Verdict below is the registered exploratory result. PI CLOSED 2026-11-23 (Richard: “we can close GC34Q”).

## Result

**PASS on all four registered gates.** Two standing preparations share the same initial scalar field and have initial time derivatives zero and one times that field in code units. The independent flat Hamiltonian gives an added kinetic energy of **39.05823652** in both directions and on both grids and their continuum references. Initial total momentum is zero to roundoff. The stress and local momentum are included in the difference source. The measured difference in the tidal tensor agrees with the response computed in continuum linearized GR. The fine-grid error is **0.031–0.771%** and improves by **3.95–4.04×** on halving the spacing (3.954–4.034 before rounding).

| Gate | Result | Verdict |
|---|---|---|
| Q0, initial conditions and coverage | Identical profiles within each pair. Initial momentum relative to the energy increment ≤ 6.4e-17. The increment matches the continuum within 1.6e-15 relative. Individual flat energies drift ≤ 1.8e-15. The registered cutoff 14 omitted up to 6.9e-6 of squared source norm, so the next registered cutoff, 18, was used. Omitted fraction ≤ 1.31e-12 in the audit. | PASS |
| Q1, constraints and timestep | Initial residual ≤ 3.9e-16. Maximum at dt = 0.05 ≤ 7.73e-7. In the off-axis coarse run, halving dt reduced 7.57e-7 to 4.73e-8 (15.99×). The final timestep effect was 1.34e-5 of the spatial discrepancy. | PASS |
| Q2, differential tidal curvature | All six fine-grid comparisons below 2%, with refinement ratios above 3. | PASS |
| Q3, alternative numerical path | Direct Fourier sums, per-mode kernel assembly, least-squares initialization and DOP853 agree with the batched RK4 calculation within 3.47e-7 relative for the three strongest modes. Two are a conjugate pair. | PASS |

## Differential tidal error against continuum linearized GR

| Direction | Time | h = 0.5 | h = 0.25 | Ratio |
|---|---:|---:|---:|---:|
| Axis | 0 | 0.1384% | 0.0343% | 4.034 |
| Axis | 2 | 0.8014% | 0.2027% | 3.954 |
| Axis | 4 | 3.0880% | 0.7709% | 4.006 |
| Off-axis | 0 | 0.1240% | 0.0308% | 4.023 |
| Off-axis | 2 | 0.5615% | 0.1419% | 3.956 |
| Off-axis | 4 | 2.3852% | 0.5939% | 4.016 |

The weighted norm of the continuum tidal difference ranged from 0.0102 to 0.0300 in code units. The comparison uses the same periodic boundary conditions, omits the homogeneous mode on both sides, and starts with matched constraint data. Evolving the difference source is equivalent to differencing the fields because the gravity equations and their initializer are linear.

## Scope and reproducibility

This establishes the expected tidal response to added initial kinetic energy in the tested scalar model, with all components of its stress-energy retained. It is a one-way weak-field test with an assigned matter–geometry coupling. The scalar evolves on the flat grid. Gravitational feedback and the origin of Y1's geometry remain open. The homogeneous gravitational response in a periodic box was outside the observable.

Run `python src/exgc34q.py`. The main path reuses GC34P's packet profile, source assembly, solver and tidal norm and GC34O's matter sources. The numerical cross-check shares the profile, matter sources and Y1 tensor formula. It changes the assembly, transform, initialization and integrator. No registered criterion was changed. The run used Python 3.13.2 and NumPy 2.2.4 on Windows, with one numerical thread and 1464 CPU seconds (0.407 core-hours) against the 2 core-hour cap.

[Full numerical record](results.json) · [Runner](../../src/exgc34q.py) · [Claude's handoff](../../docs/EX-GC34Q_handoff.md)
