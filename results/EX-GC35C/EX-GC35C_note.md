# EX-GC35C — longer reciprocal scalar evolution and constraint growth on Y1

2026-11-24 | EXPLORATORY | ChatGPT | [Registered plan](../../spec/EX-GC35C_plan.md) · [Claude review](../../review/EX-GC35AB_claude_review.md). PI closure pending.

**The perturbative scalar benchmark failed the long-time constraint gate.** In the t≤20 window the fine residual first exceeded 10⁻³ at the recorded t=5.75 sample. It reached a maximum of 3.68e-1. The vacuum spectrum rules out a GC34E-like exponential mode in the tested synchronous linear generator. The measured error grows in the sourced problem. This test does not isolate an individual nonlinear term as its cause.

## Reproduction and controls

The physical 24³ and 48³ lattices use the original GC35A Fourier supports and kernels. The fine run computes nonlinear products on a 24³ quadrature grid with physical .25-spacing translations implemented by Fourier phases. The flat packet, both metric orders and the scalar correction retain their original evolution. The maximum relative deviation from either published t=4 residual was 2.04e-14, passing the 10⁻⁵ reproduction gate. Leading constraints peaked at 7.00e-6 coarse and 7.51e-6 fine. At t=8 halving the coarse timestep changed the next-order normalized residual by a fraction of 4.55e-7. The compensated continuum control was evolved to t=20 with the same spectral cutoff and continuum symbols.

| t | Coarse normalized | Fine normalized | Coarse without compensator corrections | Continuum normalized | Fine unscaled absolute | Fine scaled second/first metric norm |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 5.70e-4 | 1.72e-4 | 5.83e-4 | 4.79e-9 | 1.73e-2 | 5.41e-4 |
| 8 | 2.36e-2 | 6.43e-3 | 2.34e-2 | 3.15e-8 | 3.07e-1 | 1.02e-3 |
| 12 | 1.74e-1 | 4.27e-2 | 1.72e-1 | 6.93e-8 | 1.68e+0 | 1.87e-3 |
| 16 | 4.29e-1 | 1.04e-1 | 4.24e-1 | 7.86e-8 | 5.57e+0 | 3.34e-3 |
| 20 | 1.56e+0 | 3.68e-1 | 1.55e+0 | 1.57e-7 | 1.42e+1 | 7.17e-3 |

Normalized residual multiplies the second-order defect by 10⁻³ and divides by the instantaneous original constraint-source norm. The second-to-first metric ratio uses both perturbative strength factors. The fine ratio peaked at 7.17e-3, below unity, so growth broke the constraint tolerance before the metric perturbative hierarchy broke down by this norm. Turning off the compensator corrections changes only the constraint sources in this implementation. The spatial mean and homogeneous metric remain omitted in both runs. The off control diagnoses sensitivity to the compensation prescription and does not model a second physical background.

At t=20 switching off the compensator corrections changed the coarse residual magnitude by 1.12% relative. Both compensated and uncompensated lattice runs show strong growth. The continuum control remains much smaller throughout the tested interval.

## Fitted growth and mode attribution

Each fit regresses the log residual against t or log t over the stated interval. RMSE values are in log-residual units. The last column fits the unnormalized constraint norm so variation in the leading source cannot by itself explain the growth. Fits describe the sampled intervals, which include oscillations, and do not identify an asymptotic law.

| Interval | Grid | Exponential slope | Exp. RMSE | Power exponent | Power RMSE | Absolute exponential slope |
|---:|---|---:|---:|---:|---:|---:|
| 2–4 | coarse | 1.57 | 0.141 | 4.58 | 0.0721 | 1.32 |
| 2–4 | fine | 1.60 | 0.155 | 4.68 | 0.0756 | 1.35 |
| 4–8 | coarse | 0.941 | 0.125 | 5.51 | 0.0448 | 0.740 |
| 4–8 | fine | 0.915 | 0.126 | 5.36 | 0.0456 | 0.713 |
| 8–12 | coarse | 0.500 | 0.0784 | 4.97 | 0.0466 | 0.444 |
| 8–12 | fine | 0.473 | 0.0769 | 4.71 | 0.0469 | 0.423 |
| 12–20 | coarse | 0.273 | 0.0564 | 4.26 | 0.0949 | 0.273 |
| 12–20 | fine | 0.272 | 0.0581 | 4.25 | 0.0970 | 0.264 |

The fitted exponential slopes fall substantially between early and late windows on both grids. There is no single corroborated rate. The power fit is also window-dependent. Refinement lowers the defect but does not stop it exceeding the operating threshold in the tested interval.

Largest fine-grid residual component by Fourier mode, with energy and momentum magnitudes in the kernel row normalization, before strength scaling. The fields are real, so mode signs can switch between conjugate partners of equal magnitude.

| t | Integer mode | Energy | Longitudinal momentum | Transverse momentum |
|---:|---|---:|---:|---:|
| 4 | (-6, -3, -2) | 9.48e-4 | 1.18e-4 | 1.75e-3 |
| 8 | (-6, -3, -2) | 5.75e-3 | 2.63e-3 | 3.58e-2 |
| 12 | (-6, -3, -2) | 7.17e-2 | 9.61e-3 | 1.63e-1 |
| 20 | (2, 2, 2) | 1.48e+0 | 1.92e-1 | 1.21e-2 |

## Synchronous linear vacuum comparison

For each nonzero retained mode let Q be the vacuum acceleration matrix in units of |K|², where K is the Y1 lattice wavenumber. On both grids the numerical identity Q²(Q+I)=0 held to within 5.29e-15 in matrix norm, with rank(Q)=3 and rank(Q²)=2. The normalized first-order generator therefore has oscillatory tensor branches and a defective zero sector in the unrestricted evolution. Numerical diagonalization of that sector is sensitive to roundoff. The apparent positive rates from direct diagonalization were at most 5.28e-4 in simulation time units, far below the finite-window fitted rates. The vacuum constraint kernel had dimension eight and was invariant to 2.37e-15 after scaling velocities by |K|. No matching exponentially growing physical re-slicing mode was identified.

In axes aligned to K, let S be the trace of the two transverse spatial metric components and L the longitudinal metric component. The unrestricted vacuum gives S¨=0 and L¨=−|K|²S/2. Thus violating the linear constraints can already produce polynomial metric drift, while the true tensor modes oscillate. The second-order test reuses that generator with an imperfect time-dependent source, unlike the original GC34E vacuum test. These results establish failure of constraint control for this source, grid range and interval. They leave all-order dynamics and other gauges open.

Recorded CPU use was 1.02 core-hours for the five evolutions plus a small vacuum audit. Raw data include [coarse](coarse.json), [fine](fine.json), [half timestep](half.json), [continuum](continuum.json), [compensator-off](off.json) and [vacuum](vacuum.json). Reproduce using [the driver](../../src/exgc35c.py) and [vacuum audit](../../src/exgc35c_vacuum.py). The moving-knot gate was not met.
