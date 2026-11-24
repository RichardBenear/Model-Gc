# EX-GC34V — next-order Y1 source compatibility (EXPLORATORY)

2026-11-24 | ChatGPT | PI authorization in the GC34U closure request. [Plan](../../spec/EX-GC34V_plan.md) and log registered before execution (1a12627 and 05bb856). PI closure pending.

**V0 PASS, V1 FAIL.** Linear Y1 remains consistent. The tested quadratic completion fails the exact source gate in both harmonic directions. No second-order metric was evolved.

The test evaluates the \(1+3\) scalar and \(2+2\) metric interactions at harmonic 4 of a single complex Klein–Gordon mode. Using a positive-frequency interaction coefficient avoids a homogeneous channel in this diagnostic. A real packet also has conjugate modes and a positive mean density, so the diagnostic does not implement isolated lattice boundaries or solve the mean-density problem for a periodic packet.

The tested effective mixed source is \(F_2=T_{1,\mathrm{scalar}}+T_{1,\mathrm{metric}}-G_2\), with \(\kappa=1\). It uses GC34O's exactly balanced flat source for the scalar variation and GC34R's scalar force. It uses phase-centered Y1 derivatives in the continuum metric-dependent matter term and in the quadratic Einstein expression, taking derivatives after products at their output harmonics.

The first-order forced residual was at most \(8.5\times10^{-16}\). Both continuum control residuals were below \(5\times10^{-17}\). At finite spacing:

| Direction | h | Normalized divergence | Normalized antisymmetric part after lowering | Separate geometric Bianchi defect |
|---|---:|---:|---:|---:|
| Axis | 0.5 | 3.57e-3 | 2.20e-2 | 3.64e-3 |
| Axis | 0.25 | 8.60e-4 | 5.37e-3 | 8.46e-4 |
| Axis | 0.125 | 2.13e-4 | 1.33e-3 | 2.07e-4 |
| Diagonal | 0.5 | 2.19e-3 | 2.25e-2 | 3.14e-3 |
| Diagonal | 0.25 | 4.67e-4 | 5.32e-3 | 6.98e-4 |
| Diagonal | 0.125 | 1.12e-4 | 1.31e-3 | 1.69e-4 |

**V2 diagnostic.** The geometric and matter terms separately violate their expected connection-corrected balances. The relative product-rule defect \( \|K(2q)-2K(q)\|/\|2K(q)\|\) falls from 0.0192 to 0.00120. The staggered symbol satisfies \(K_a(2q)=2K_a(q)\cos(q_a/2)\). Merely substituting it into continuum product identities does not preserve them exactly. The residuals converge roughly quadratically, but the requested linear source must be symmetric and conserved at finite spacing.

## Consequence for Y1

This rules out the tested completion, including a direct evolution with that candidate's source. It does not rule out an action with appropriately derived lattice products, field placements, and nonlinear symmetry transformations. An exact identity to use in building that candidate is the link product rule \(D(fg)=\bar f\,Dg+\bar g\,Df\), where each bar is the average over the two endpoints. Whether the resulting gravity and matter terms satisfy the combined identity is open. GC34U remains the isolated reciprocal continuum benchmark. A real-packet lattice test also needs consistent boundary conditions and constraint initialization.

Run `python src/exgc34v.py`. Python 3.13.2, NumPy 2.2.4, Windows, 1.30 CPU seconds. The continuum control shares the tensor assembly but uses continuum derivatives, analytic scalar stress, and the continuum scalar variation. No gate changed. [Numerical record](results.json) · [Runner](../../src/exgc34v.py)
