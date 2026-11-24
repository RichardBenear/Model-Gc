# ChatGPT review of EX-GC34N

2026-11-23 | Independent review of the exploratory results

**Verdict.** Keep Y1 as the linear vacuum benchmark. Its staggered geometry addresses the vacuum problems found in J1, and its scalar momentum balance improves substantially. The scalar still fails the registered local momentum gate. Resolve that source mismatch before attempting a coupled scalar–gravity evolution.

## Vacuum and prescribed sources

Y1 puts the shift and spatial displacement on links, normal metric components at vertices, and shear on plaquettes. Six metric components replace the seven Kuhn edge lengths. The vacuum kernel implements the linearized Einstein equations using the consistent staggered derivative \(K_a=2\sin(q_a/2)/h\). Both coordinate symmetries hold to roundoff. The two gravitational-wave polarizations and the massless scalar share \(\omega=|\mathbf K|\), including its finite-spacing directional dependence.

The static benchmark gives \(\gamma=1\) for the tested modes. The prescribed moving-source curvature error is \(1.04\times10^{-4}\) at \(kh=0.05\) and drops by about four on halving the spacing. This addresses J1’s first-order error for the prescribed source. The higher-precision root calculation is a reasonable resolution of the apparent growth from numerical splitting of zero roots. The registered time-domain run also passed with constraint-compatible seeds.

## Scalar matter

The energy identity held to roundoff. The local momentum residual dropped from roughly 14–20% in GC34M to 0.2–0.4% at \(h=0.25\). The fine refinement step is consistent with second-order convergence, but three of six coarse-step comparisons missed the registered factor of 3.5. E-M2 and the overall study therefore fail as registered. A finer rerun would clarify convergence. Even clean second-order convergence would leave a finite-spacing mismatch against the exact local source identities of the vacuum equations.

## One-dimensional algebraic check

Let \(G_v=(\phi_{v+1}-\phi_v)/h\). In the one-dimensional limit, the gradient contribution to the normal stress, in the normalization \(2s_{aa}/h^3\), is \((G_v^2+G_{v-1}^2)/4\) in the code. Differentiating the existing link momentum using the discrete Klein–Gordon equation gives exact local balance if that contribution is \(G_v G_{v-1}/2\). The difference is

\[
\frac{G_v^2+G_{v-1}^2}{4}-\frac{G_v G_{v-1}}{2}
=\frac{(G_v-G_{v-1})^2}{4}.
\]

This is second order for a smooth field and offers a concrete explanation for the axial residual. This calculation covers the one-dimensional limit. The transverse terms and the full three-dimensional matter action remain to be derived.

## Recommended next study

Hold Y1’s vacuum equations and the flat scalar evolution fixed. Derive the stress needed by the existing link momentum in all three directions, including the shear terms. Check whether it comes from a single local matter action that preserves the exact energy identity. Use the one-dimensional identity as an algebraic control and the generic-direction packet as the first numerical test. Gate a coupled matter–gravity run on the full local source conservation check.

## Scope and records

I reviewed the registered plan, results, and implementation on Claude’s branch and checked the one-dimensional identity algebraically. I did not rerun GC34N.

- [Registered plan](https://github.com/RichardBenear/Model-Gc/blob/claude/project-thread-2twuh8/spec/EX-GC34N_plan.md)
- [Results and verdicts](https://github.com/RichardBenear/Model-Gc/blob/claude/project-thread-2twuh8/results/EX-GC34N/EX-GC34N_note.md)
- [Vacuum kernel and scalar source code](https://github.com/RichardBenear/Model-Gc/blob/claude/project-thread-2twuh8/src/exgc34n.py)
