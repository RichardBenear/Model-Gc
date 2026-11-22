# EX-GC30B2 — Weyl tetrad universality (EXPLORATORY)

**PI CLOSED 2026-11-22 (Richard: "I think those runs can be closed"; Claude review review/EX-GC30A-B3_claude_review.md).**

Topological Weyl nodes survive all 30 registered backgrounds (120 node instances), but their reconstructed effective metrics do not generally agree. A four-node single-band-pair model has exactly equal spatial metrics on a symmetry-restricted family, including shared shear and axial deformations. Breaking its extra root-reflection symmetry while preserving inversion and all topological windings splits the metrics. Two inequivalent orbital flavors fail even in the natural homogeneous vacuum and only agree after calibration; density-like common mass shifts break that agreement. Topology protects the existence and integer chirality of nodes, not the common geometry. This is a reference lattice result, not a no-go for other dynamical mechanisms of universal metric emergence.

Design/run ChatGPT, 2026-11-22; PI authorized GC30B2. Preregistered spec/EX-GC30B2_plan.md (56da0991) and LOG (11e58b75). One thread, default 2 core-hour cap, corrected production 19.79 CPU seconds. No GC30C or other later stage executed. Preliminary import failed before any calculation because missing j loop in direction list; corrected to the registered axes/face-diagonal list, no physics or criteria change. <1 CPU s additional startup.

## Models and independent measurements

Reference Bloch Hamiltonian in unit lattice/hbar H(k)=d0 I+d.sigma, d_xy=T(sin kx,sin ky), d_z=t_z f(cos kz)+3(2−cos kx−cos ky), d0=u sin kx+v cos kz. S: a single two-band Hamiltonian with f(x)=(x−m)^2−a^2 (baseline m=0,a=.5,t_z=1), four nodes at (0,0,+/− arccos(m+/−a)); F: two conserved orbital flavors, each f(x)=x−m_a with m_a=.2,.6, two nodes each; natural t_z=1, deliberately calibrated t_z=1/sqrt(1−m_a²). No individual post-calibration retuning. Distinct conserved flavors are deliberately block diagonal; no hybridization, interactions or quantum vacuum backreaction modeled. All changes are uniform model coefficients, not self-consistent gravitational strain fields. Model lacks microscopic boost symmetry (potential exit from Weinberg–Witten premise, not a demonstrated tensor graviton).

Linear expansion around a node: E_±−E_node=w.q ± sqrt(q^T B q), J=partial d/partial k, B=J^T J, w=grad d0. This defines the normalized inverse-metric cone, with coefficients [-1, -w^T; -w, B−ww^T]. Conformal scale, spin connection and massive clocks not determined. Numerical diagonalization of complex 2x2 Bloch matrix at ±h in ten directions; symmetric band splitting squared/(4h²) averaged over signs gives quadratic spatial matrix; odd band mean gives tilt. Nine directions fit; body diagonal holdout. Richardson extrapolation of windows .02/.01 and .01/.005; compare analytic J/w. Max B error 8.97e-10 relative, w error 1.09e-12; max window differences 1.35e-8 and 1.64e-11; holdouts 2.23e-10 and 4.89e-13. All far below .002 threshold. All 120 nodes isolated: transverse zero condition requires kx,ky=0 or pi; Wilson term at pi >=6 exceeds possible negative contribution, roots strictly inside (-1,1), nonsingular J. Node residual <1e-12. Cube-surface map d/|d| winding from oriented triangle solid angles, 6 and 12 subdivisions and radii .02/.04, 480 checks. Max winding integer residual 4.0e-15, minimum boundary norm .01439. Charges (degree of d-map, Berry band sign convention separate) sum zero for each case, sign matches det J; all cones type I (max w^T B^-1 w=.01). This proves local winding for specified models, not real-space soliton stability.

Metric mismatch = largest pairwise Frobenius norm of B difference / mean tensor norm across all four nodes; tilt mismatch = difference / mean characteristic speed (sqrt(mean trace/3)). Both must <=1% to pass kinematic gate; energies independently recorded. No species-specific coordinate or time rescaling permitted. All node pairs measured in same k coordinates, but absolute node energy differences can also matter to having a shared vacuum.

## Results

| Arm / background | Max spatial mismatch (%) | Max tilt mismatch (%) | Winding preserved |
|---|---:|---:|---|
| S baseline, transverse or axial scale, shear, combined | ~0 | ~0 | yes |
| S root separation a=.4/.6 | ~0 | ~0 | yes |
| S root center m=±.025 | 3.124 | 0 | yes |
| S root center m=±.05 | 6.249 | 0 | yes |
| S root center m=±.1 | 12.511 | 0 | yes |
| S m=.05 plus shear/axial | 6.903 | ~0 | yes |
| S common sine tilt u=.1 | ~0 | ~0 | yes |
| S common cosine tilt v=.1 | ~0 | 18.091 | yes |
| F natural untuned | 19.623 | 0 | yes |
| F tuned, common shear or axial | ~0 | ~0 | yes |
| F common offset −.025 / +.025 | 2.067 / 2.144 | 0 | yes |
| F common offset −.05 / +.05 | 4.060 / 4.371 | 0 | yes |
| F common offset −.1 / +.1 | 7.841 / 9.083 | 0 | yes |
| F offset .05 plus shear/axial | 4.770 | ~0 | yes |
| F flavor 2 transverse scale 1.1 | 15.993 | 0 | yes |
| F sine tilt | ~0 | 0 | yes |
| F cosine tilt | ~0 | 19.596 | yes |

S: 8/16 total cases pass both metric tests; all 7 root-reflection symmetric untilted cases plus sine-tilt case pass. S common offset cases break extra root-reflection but preserve inversion, no topology changes. F: 4/14 total pass (calibration and common shear/axial/sine tilt); natural, relative and calibration excluded from robust generalization, remaining 3/11 common tests pass. Both fail full registered C3 robustness gate irrespective of counting controls. Cosine tilts yield nonzero energy spread .1 (S), .04 (F); they remain type I and preserve inversion. Thus same spatial tensor is not sufficient for a common spacetime cone. S symmetry restricted result exactly reproducible, not generic proof of equal metrics. F natural spatial mismatch from B_z=.96 vs .64, no independent inter-flavor protecting symmetry.

## Why topology is insufficient

For S, at cos kz=m±a, B_z=4t_z²a²[1−(m±a)²]. Difference in B_z has magnitude 16t_z²a³|m|. At m=0 the k_z→pi−k_z root reflection maps the two families, and reflection k_z→−k_z relates partners; in our restricted separable Hamiltonian B has no xz/yz entries. Equal squared slopes follow exactly; common shear changes T^T T for every node equally. At m≠0, root-reflection symmetry is lost but inversion H(−k)=sigma_z H(k)sigma_z remains when u=0, and each integer winding remains unchanged. A nonzero m is an allowed microscopic perturbation, not an instability. This is an explicit counterexample to topology alone fixing the tensor. Inversion relates coordinates at opposite momenta, not equality of all coordinate components in arbitrary Hamiltonians. Cosine tilt also gives w_z opposite signs at inversion partners; sine tilt u gives same w_x, breaking inversion but retaining root-family equivalence in this ansatz. Even symmetry-related excitations need spacetime transformation behavior checked, not just energy degeneracy.

For F after calibration, B_z=[1−(m_a+delta)^2]/(1−m_a²). Uniform shift delta changes the two slopes differently. Distinct microscopic flavors can have equal velocities at one point without sharing their response derivative. Compare GC30B lambda=n*g12 tuning; analogous structural failure. No inference that all interacting Weyl or tetrad theories must fail. The positive result says what to protect: a common tensor in every distinct sector, not merely an integer topological invariant. Real-space stable particles, fermionic matter, gauge waves, massive probes and composites must eventually be independently tested. No tensor collective mode, emergent diffeomorphism symmetry, energy sourcing, or 1/r potential established.

## Proposed direction (not authorized)

Do not proceed to GC30C with this model. Do not spend more on this exact free-band counterexample. Consider GC30B3: dynamical metric locking rather than static common microscopic parameters. Primary reference Anber & Donoghue, https://arxiv.org/abs/1102.0789, finds interaction-induced IR velocity convergence in examples (often slow logarithmic); Roy, Juricic, Herbut, https://arxiv.org/abs/1510.07650, finds common IR speed in specific interacting critical fermion/boson/gauge models. This suggests an alternative mechanism in which initially inequivalent sectors may flow toward a common cone, not one imposed or protected solely by topology. First analytical screen: derive/verify published beta functions and attraction basin; distinguish irreducible anisotropy/tilt from isotropic speed; assess finite scale range, masses/thresholds, weakly coupled neutral sectors, backreaction and radiative stability. A common uniform IR cone is still not a dynamical metric or gravity. If this fails, reserve genuinely pregeometric causal/entanglement approaches rather than invent further tuned local media. Sakharov remains downstream once geometry and universal coupling are justified.

## Sources and reproducibility

Volovik & Zubkov, emergent Weyl spinors, https://arxiv.org/abs/1402.5700; Cortijo & Zubkov, strain-dependent tetrad in cubic tight binding, https://arxiv.org/abs/1508.04462. These are motivation, not replication of their specific Hamiltonian. Explicit models above constructed for GC30B2. Run `python src/exgc30b2.py`; compact case files cases_0.json through cases_5.json include all fit-point observables, node locations/energies, winding, tensors, tilt, errors, pairwise comparisons, package versions and CPU time. Same-author analytic check is not independent replication. PI closure pending.
