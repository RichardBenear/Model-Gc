# EX-GC34T — next-order conservation and the periodic mean (EXPLORATORY)

2026-11-24 | ChatGPT | [Plan](../../spec/EX-GC34T_plan.md) and log registered before execution (cf15388 and 0309851). PI closure pending.

**T1 PASS and T2 PASS.** This verifies the continuum identities and a homogeneous expanding benchmark. It does not supply a second-order lattice action or validate an inhomogeneous reciprocal solution on Y1.

## T1. What must cancel

Use mixed-index tensors and multiply the gravitational coupling by \(\epsilon\). Expand
\[
G=\epsilon G_1+\epsilon^2(G^{(1)}[h_2]+G_2)+\cdots,\qquad
T=T_0+\epsilon T_1+\cdots .
\]
Here \(G_1=G^{(1)}[h_1]\), \(G_2=G^{(2)}[h_1,h_1]\), and \(T_1\) includes the metric dependence of the matter tensor and the scalar response. Write
\[
\mathcal C_{\Gamma_1}(X)_\nu=\Gamma_1^\mu{}_{\mu\lambda}X^\lambda{}_\nu
-\Gamma_1^\lambda{}_{\mu\nu}X^\mu{}_\lambda .
\]
Then the expanded Bianchi and matter identities give \(\partial_\mu G_2^\mu{}_\nu=-\mathcal C_{\Gamma_1}(G_1)_\nu\) and \(\partial_\mu T_1^\mu{}_\nu=-\mathcal C_{\Gamma_1}(T_0)_\nu\). The effective source \(\kappa T_1-G_2\) has zero flat divergence if the same first-order metric and source satisfy \(G_1=\kappa T_0\).

The symbolic runner builds \(G_1,G_2\) from the connection and Ricci tensor for the registered periodic manufactured metric. It finds exact cancellation in all four components. The manufactured source is defined from \(G_1\). This tests the geometric identity, while the matter identity follows separately from covariant conservation on the scalar equations of motion.

As a negative control, deleting a constant positive density from the leading matter tensor alone leaves an energy mismatch proportional to \((2t+1)\sin x+\tfrac12\cos x\) in the registered example. Its mean vanishes, but the mismatch is present in nonzero modes. Thus deleting a Fourier mean on the source side is insufficient as a next-order conservation prescription. The test leaves open a separately modeled, covariantly conserved compensator.

## T2. What to do with the mean

The flat periodic linear Hamiltonian constraint is
\[
\partial_i\partial_j h_{ij}-\Delta h_{ii}=2\epsilon\kappa\rho .
\]
The spatial integral of the left-hand side is zero. A positive mean density cannot satisfy this constraint on a static flat periodic background. In the homogeneous expanding three-torus benchmark, the exact constraint and evolution equations instead give
\[
3H^2=\epsilon\kappa\rho,\qquad
\dot H=-\tfrac12\epsilon\kappa\dot\phi^2,\qquad
\ddot\phi+3H\dot\phi+m^2\phi=0 .
\]

| Coupling parameter | Initial \(H\) | Max relative constraint residual |
|---:|---:|---:|
| \(10^{-4}\) | 0.00577350 | \(4.95\times10^{-12}\) |
| \(2.5\times10^{-5}\) | 0.00288675 | \(4.88\times10^{-12}\) |

The factor-of-two change in \(H\) confirms the square-root scaling. An expanding background is viable in this homogeneous example. It changes the hierarchy of the weak-field expansion: at fixed density the expansion rate is order \(\sqrt{\epsilon}\), so it cannot be restored as an order-\(\epsilon^2\) add-on to the previous flat calculation. Setting \(H=0\) gives a unit relative constraint error at both couplings.

## Implication for the next reciprocal test

For a localized packet, the shortest physically clean route is a small continuum benchmark with isolated boundary conditions and a consistent initial constraint solve, so the field of the positive total energy is retained. A periodic route needs an expanding background or a compensating sector whose stress and exchange are included in both orders. After choosing the background, verify the full second-order effective source on that packet before evolving the correction; an eventual Y1 implementation must verify its own discrete identity. This study does not establish that the continuum cancellation survives a chosen lattice discretization.

Run `python src/exgc34t.py`. Python 3.13.2, SymPy 1.14.0, 1.06 CPU seconds against the two-core-hour cap. No registered criterion changed.

[Results](results.json) · [Runner](../../src/exgc34t.py)
