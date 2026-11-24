# EX-GC34S — direct reciprocal sourcing fails the Y1 constraints (EXPLORATORY)

2026-11-24 | ChatGPT | PI go-ahead: “go ahead with the remaining step”. [Plan](../../spec/EX-GC34S_plan.md) and log were registered before production. PI CLOSED 2026-11-24 (Richard: “I sign off on all of the open items”).

## Registered result

**S0 PASS, S1 DIRECT RESOURCING RULED OUT, S2 PASS.** The scalar correction computed in GC34R changes the source. Feeding that change into the unchanged linear Y1 metric equations gives a constraint violation that survives timestep and grid refinement. This rules out this direct prescription for reciprocal evolution. GC34R's leading-feedback result retains its original scope.

| Grid | Timestep | Max original constraint residual | Max constraint residual with resourcing at strength \(10^{-4}\) |
|---|---:|---:|---:|
| \(h=0.5\) | 0.025 | \(1.24\times10^{-8}\) | \(7.45\times10^{-5}\) |
| \(h=0.5\) | 0.0125 | \(7.78\times10^{-10}\) | \(7.45\times10^{-5}\) |
| \(h=0.25\) | 0.025 | \(1.34\times10^{-8}\) | \(7.62\times10^{-5}\) |

The flat scalar's local balances hold within \(5.1\times10^{-15}\) in the audited norms. At \(t=2\) on the fine grid, the tangent energy and momentum source defects are 2.476 and 0.939 relative to their respective zeroth-order time-derivative norms before multiplication by the weak-field parameter. Halving the timestep changes the measured normalized tangent defects by at most \(2.16\times10^{-8}\) relative. At strength \(5\times10^{-5}\), the source defect divided by strength agrees with the tangent defect to within \(5.0\times10^{-5}\) relative.

## Local product-law cross-check

Let \(\phi_0\) be the flat packet and \(\psi\) its leading response in GC34R. Define \(D=-m^2+\Delta_h\), and write the scalar correction equation as \(\ddot\psi=D\psi+F\), with \(F\) the acceleration induced by the first-order metric. The GC34O flat-action source evaluated on \(\phi_0+\epsilon\psi\) then has a flat scalar equation residual \(\epsilon F\). The first-order local energy balance error is \(h^3\dot\phi_0 F\). The first-order momentum balance error on a link is \(h^3\overline F_a G_a\), where \(\overline F_a\) is the endpoint average of \(F\) and \(G_a\) is the original field gradient on the link. These products reproduce the measured tangent defects within \(1.5\times10^{-14}\) relative. The defect is the exchange with the metric missing from the flat source prescription.

## What a reciprocal test still needs

At the order where the corrected scalar sources gravity, the continuum expansion requires both the metric dependence of the matter source and the quadratic gravitational terms. Schematically, if a parameter \(\epsilon\) multiplies \(\kappa\), and \(h_1,h_2\) are the coefficients of \(\epsilon,\epsilon^2\) in the metric,

\[
G^{(1)}[h_2]=\kappa T^{(1)}[\phi_0,\psi,h_1]-G^{(2)}[h_1,h_1].
\]

The lattice counterpart needs a matched source and gravity completion whose combined right-hand side obeys the constraint identities at that order. A separate choice is required for the periodic homogeneous mode. GC34P–S omit it from the linear solve. In a periodic expansion about a static flat background, the spatial average of the linear Hamiltonian constraint is zero, so positive mean matter density cannot be included in that constraint without changing the background or the boundary or source prescription. A reciprocal experiment must account for this consistently at the next order.

The next milestone is to verify the next-order conservation identity and homogeneous-mode prescription before trusting an evolved reciprocal solution. The current run does not assign a physical interpretation to the inconsistent metric correction.

## Reproducibility

Run \`python src/exgc34s.py\`. Shared dependencies are the GC34R scalar tangent force and packet, GC34O matter sources, and GC34P Y1 solver. The product-law check is a separate algebraic path against the full discrete source derivative. No registered criterion changed. Python 3.13.2 and NumPy 2.2.4 on Windows. 151.7 CPU seconds (0.042 core-hours against a 2 core-hour cap).

[Full numerical record](results.json) · [Runner](../../src/exgc34s.py)
