# Causal invariance and coordinate symmetry after GC34

2026-11-24 | Paper-only direction note | [Claude's GC34Q–Z review](../review/EX-GC34Q-Z_claude_review.md)

## Which symmetry does the programme need?

Causal invariance concerns the order of local updates. Changing a legal schedule should leave the prescribed physical outcome or event history equivalent. Lattice coordinate symmetry concerns redundancy in a description of spacetime. It relates metric variables representing the same geometry and supplies identities that connect evolution equations with their constraints. They are distinct requirements. A schedule-independent rule can still give different propagation geometry to different species, as GC31A demonstrates. Conversely, a free metric action with exact lattice gauge null directions can be implemented with a fixed time integrator and gives no proof about the order independence of a microscopic rule.

The emergent-gravity programme can therefore keep causal invariance as a target while using a benchmark whose nonlinear coordinate symmetry is approximate at finite spacing. This is a conditional methodological judgment. To be adequate, the benchmark must exhibit controlled constraint errors under spatial and timestep refinement, two physical long-wavelength tensor modes, and compatible matter sourcing. Bounded residuals on a short run establish a finite-time regime only. Long-time drift, extra modes, or a residual that fails to converge would limit the benchmark. Exact coordinate symmetry would strengthen it, but GC34V–Z showed a specific set of candidates fail. Those tests give no general impossibility theorem.

[GC31A counterexamples](../results/EX-GC31A/EX-GC31A_note.md) · [GC34Z compatibility test](../results/EX-GC34Z/EX-GC34Z_note.md) · [Wolfram's rewriting framework](https://arxiv.org/abs/2004.08210)

## First-order geometry as an alternative

In a tetrad-and-connection formulation, local frames express lengths and angles. Link transport relates the frames, and closed products of transports encode curvature. With a suitable group action on frames and links, **internal local Lorentz symmetry** can be exact at finite spacing. Closed-loop transport transforms by conjugation, which allows invariant action terms. This is a different symmetry from moving or re-slicing the coordinate locations of the cells. A lattice action with exact internal symmetry still needs a separate demonstration of diffeomorphism or constraint behavior. A topological BF construction may have larger symmetries, but its gravity constraints must remove the topological sector's excess freedom while leaving the two gravitational modes and appropriate matter coupling.

Relative to Y1, this route adds frame and connection variables, transport conventions, geometric gluing conditions and torsion or simplicity constraints. It also requires a new matter action, boundary and constraint initialization, and a new degree-of-freedom audit. Holonomies provide a structural way to organize curvature. The immediate test is whether they also give the needed gravitational constraint identities on the configurations actually used. Regge calculations show how exact linearized symmetries about flat backgrounds can coexist with broken symmetry away from them.

[First-order and BF formulations](https://www.maths.tcd.ie/EMIS/journals/SIGMA/2012/032/) · [Discrete gluing and twisted geometries](https://arxiv.org/abs/1001.2748) · [Regge constraints and broken symmetry](https://arxiv.org/abs/0905.1670) · [Linearized canonical Regge constraints](https://arxiv.org/abs/0912.1817)

## Decision

Use linear Y1 and the isolated GC34U reciprocal continuum result as reference problems. Give the approximate completion a bounded numerical test and reconnect it to the GC26 knot. Keep exact cubic-symmetry fitting paused. Reopen a structural redesign only for a concrete failure the benchmark needs to address or a derivation with an independently testable identity. This leaves the mechanism by which microscopic Model Gc rules produce a shared geometry as the larger open question.
