# EX-GC34U — an isolated continuum benchmark for reciprocal scalar gravity (EXPLORATORY)

2026-11-24 | ChatGPT | PI go-ahead: “the next reciprocal test you mentioned in 34T can be ran”. The [plan](../../spec/EX-GC34U_plan.md) and log were committed before execution (e51dc2e and 6956218). PI closure pending.

**U0, U1 and U2 PASS on the corrected numerical implementation.** With isolated boundaries a spherical massive scalar sources the metric and evolves on that same metric. The resulting change in the scalar also changes the enclosed gravitational mass. The independently evaluated mass-flux equation converges, including the first coupling derivative that probes the next gravitational order.

## Model and order

This is a radial discretization of the continuum Einstein–Klein–Gordon system, with the polar-areal metric \(-\alpha^2 dt^2+a^2dr^2+r^2d\Omega^2\). The radial constraints are solved from the evolving scalar at every RK4 stage. The initial field and canonical momentum are the same positive Gaussian profile. Exterior time is fixed by matching the lapse to Schwarzschild at the outer radius, where the field remains negligible. The positive total mass is retained. Spherical symmetry makes this an inexpensive continuum reference and excludes gravitational-wave degrees of freedom.

Write \(a^{-2}=1-2\lambda M/r\). The geometric mass is \(m_{\rm geom}=\lambda M\), so the derivative \(M_1=\partial_\lambda M|_0\) contributes to the metric at order \(\lambda^2\). The diagnostic differentiates the radial mass constraint in the scalar evolution direction, using a separate quadrature, and compares it with the independently specified Einstein flux
\[
\dot M=\frac{r^2\alpha\Pi\Phi}{2a^3}.
\]
The flux equation is unused in the evolution. Its coupling derivative is checked using the centered difference at \(\lambda=\pm0.005\). Negative couplings are used only to estimate the weak-coupling derivative. The fixed-scalar control evaluates the nonlinear mass constraint on the flat scalar history, so subtracting it isolates the contribution from changing the scalar rather than from the metric factor in the mass constraint alone.

## Registered results

- **U0:** At \(\lambda=0.01\), the initial fine-grid geometric mass is 0.00336510, the maximum compactness is 0.001689 and relative mass drift is \(8.88\times10^{-6}\). The boundary scalar amplitude is below \(3\times10^{-22}\) for the positive-coupling runs. Doubling the outer radius changes the scored interior profiles by at most \(1.1\times10^{-13}\) relative.
- **U1:** At \(h=0.03125\), the largest raw relative flux discrepancy is \(3.22\times10^{-4}\), and the largest coupling-derivative discrepancy is \(3.79\times10^{-4}\). Both improve by approximately 16× from \(h=0.125\) to 0.03125, consistent with second-order convergence. This tests the spherical constraint propagation condition rather than the full tensor identity on a 3D lattice.
- **U2:** At \(t=4\), the feedback contribution to the mass derivative inside \(r\le8\) has normalized norm 0.04681. The profiles estimated with the two coupling strengths differ by \(2.49\times10^{-5}\) relative, and the two finer grids by 0.00105. The differentiated Hamiltonian constraint predicts the full inner mass-derivative profile within \(9.59\times10^{-5}\) relative. The signed feedback contribution to \(M_1\) at \(r=2\) is \(+0.027293\). The corresponding leading contribution to geometric mass is \(\lambda^2\) times that coefficient.

## Implementation correction and limits

The first execution used a centered pointwise derivative of the radial flux with a separate origin stencil. It developed numerical growth near the origin and non-finite values on the fine grid. Its [failed output](first_attempt.json) is retained, with non-finite values written as JSON null. Replacing that divergence with a conservative spherical shell flux, including the center half-cell, removed the failure. The physical problem, grids, couplings, observables, and gates were unchanged. The registered verdict refers to the corrected implementation and the raw first attempt remains separate.

This benchmark establishes reciprocal weak scalar gravity for the isolated, spherically symmetric continuum equations. It supplies a next-order test for a future Y1 discretization. The nonlinear staggered action and its discrete conservation identities remain to be constructed and checked.

Run `python src/exgc34u.py`. Python 3.13.2 and NumPy 2.2.4 on Windows. Corrected run 4.06 CPU seconds, less than 10 CPU seconds including the failed attempt, within the 2 core-hour cap. The cross-check shares the scalar evolution and radial field data, but uses an independent radial quadrature and the unused flux equation.

[Numerical record](results.json) · [Runner](../../src/exgc34u.py) · [GC34T prerequisite](../EX-GC34T/EX-GC34T_note.md) · [Background on polar-areal Einstein–scalar evolution](https://pmc.ncbi.nlm.nih.gov/articles/PMC5253925/)
