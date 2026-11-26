# EX-GC26B2 — reflection-symmetric numerical inertia (EXPLORATORY)

**Summary.** The eight-orientation reflection-symmetric stencil almost eliminates unwanted transverse momentum, but does not fix the longitudinal energy-flow/momentum mismatch. A new exact semidiscrete energy-current calculation explains this mismatch in both the old and new stencil: adjacent-site velocity enters energy flux, whereas the canonical momentum estimator uses same-site velocity. Fine-grid errors still exceed the .5% ratio target. A longer coast passes the speed-loss criterion but fails the global momentum criterion. A known weak external force produces direction-symmetric but biased acceleration; reducing coupling to the tail improves the inferred mass from .783 to 1.056 E0 without passing the gate. Thus precision inertia/gravity calibration is not ready, even though all ten cases conserve energy and retain geometric degree. No new physical inertia claim. Recommend targeted current-consistent higher-order discretization and a localized weak-force response convergence test rather than progressing to precision 26E.

PI Richard Benear; design, implementation and runs ChatGPT; review Claude pending. PI CLOSED 2026-11-26. Plan prospectively published to main 897db1c5c307cc9d2f75397b045daed324b70c59, base f0b4e8b. Plain massless Skyrme theory, numerical construction ours, no G or chemistry. See `spec/EX-GC26B2_plan.md`, `review/EX-GC26B_response.md`.

## Review and method

Claude reproduced 26B with its original NumPy solver. Axial velocity relative to P/H, not H or P by themselves, exposes numerical error. The two-orientation all-forward/all-backward stencil preserves permutation and full inversion but not individual reflections. Here average all eight combinations of nearest-neighbor difference signs, with independent orientation masks and weight 1/8. No higher-order derivative or pion mass. The same orientation set is used in the static energy and exact kinetic action. C++ derives the constrained acceleration from the discrete Lagrangian, including external potential when used. RK4 projection, dt=.05h. All rest states freshly relaxed from hedgehog, no asymmetric rest carried forward.

Independent NumPy reference sums all eight orientations for energy, gradient, forced acceleration. Agreement 7.77e-16, 4.16e-17, 2.84e-14 respectively. Tangent finite difference relative gradient error 1.34e-8. Single-coordinate mirror equivariance: acceleration 2.49e-14, H 0, P 2.34e-19. Vacuum energy and acceleration zero even with external force, checkerboard energy .311024 >0, numerical directional H+V rate zero to printed precision. This validates implemented symmetry and discrete equations, not the continuum approximation. Geometric degree evaluated by three independent regular targets as in A2. No rescore of B.

## Rest calibration

| h | nominal L | E0 | E2/E4 | max projected gradient | degree |
|---:|---:|---:|---:|---:|---:|
| .2 | 7 | 1.202402 | .879589 | 7.95e-6 | −1 |
| .15 | 7 | 1.217104 | .930178 | 7.17e-6 | −1 |

Both reach 300-iteration cap, not formal convergence. Energies nearly identical to 26B, but independent static error remains, including substantial Derrick imbalance. Fine actual L=6.975. The primary diagnostic is fitted energy-centroid velocity divided by mean projected P/H, regardless of launch velocity or rest-energy error. Compare H/(γ(u)E0), P/(γ(u)uE0), inferred M/E0, transverse P, global momentum drift and geometric endpoint degree; use first/last 1 time unit for coast and matched T=.8 cases. Regional continuum-stress flux as in B remains approximate, not a discrete conservation identity.

## Matched short boosts

| h | u | direction | v/(P/H) | M/E0 | H/(γ(u)E0) | transverse P/E0 max | global ΔP/E0 max |
|---:|---:|---|---:|---:|---:|---:|---:|
| .2 | .1 | axis | .985832 | 1.014643 | 1.000271 | 2.85e-9 | 5.47e-5 |
| .2 | .3 | axis | .982277 | 1.021430 | 1.002491 | 1.88e-9 | 6.24e-4 |
| .2 | .3 | diagonal | .984900 | 1.018209 | 1.002530 | 5.04e-9 | 2.08e-4 |
| .15 | .1 | axis | .992871 | 1.007315 | 1.000148 | 6.91e-10 | 4.81e-5 |
| .15 | .3 | axis | .991040 | 1.010730 | 1.001377 | 1.86e-9 | 2.19e-4 |

Compare B h=.2 axis .3: transverse P/E0=.003673, now 1.88e-9 (reduction ~2 million-fold). At same h and u=.3 directional mass split falls from 1.72 to .322 percentage points, but the old diagonal estimate had smaller absolute error (1.004325). All primary .5% ratio gates fail; only fine axial .1 has mass error <1%. Refinement roughly halves error but reflection symmetry alone is not sufficient. H drift ≤1.42e-11 across short runs, all degrees −1. No third resolution as preregistered; no retrospective gate change.

## Post-hoc exact energy-current diagnostic (not a gate change)

For the discrete action L = Σ_i L_i(φ_i,φ_neighbors,v_i), local H_i=v_i·π_i−L_i gives dH_i/dt=Σ_j(v_i·∂L_j/∂φ_i−v_j·∂L_i/∂φ_j), because the constraint force is orthogonal to v_i. With fixed vacuum walls and no external force, dΣ_i x_i H_i/dt = Σ_i,j (x_j−x_i)v_j·∂L_i/∂φ_j. For our nearest-neighbor derivative, the displacement s h cancels the derivative coefficient s/h. Thus exact energy current J_a=Σ_i,s w (∂L_i/∂d_a)·v_neighbor. The continuum-motivated canonical momentum estimator is the same expression with v_i rather than v_neighbor. J/H equals instantaneous energy-centroid velocity exactly for the semidiscrete unforced system; no exact translation Noether symmetry exists here. This is not a redefinition of physical P or mass to pass a gate.

Numerically verified on a random tangent field: finite-differenced energy moment derivative agrees with analytic J to 4.26e-11. Initial data for each recorded boost yield:

| h | u | direction | J/P at launch |
|---:|---:|---|---:|
| .2 | .1 | axis | .985832 |
| .2 | .3 | axis | .982153 |
| .2 | .3 | diagonal | .984904 |
| .15 | .1 | axis | .992860 |
| .15 | .3 | axis | .990981 |

Coarse observed fitted ratios .985832, .982277, .984900 match the initial J/P to at most 1.25e-4. Thus the main short-coast mismatch can be accounted for by distinct node-vs-neighbor velocities in the two observables, independently of rest mass or centroid fitting. This survives restoration of reflection symmetry; spatial truncation, not merely missing mirrors, is implicated. Exact current is useful for diagnosing and designing a better stencil, but matching H centroid by construction cannot serve as independent inertia evidence. Source `src/exgc26b2_current.py`, diagnostic `discrete_current.json`, compute included.

The same expression applied to all twelve inherited B launches (using their saved rest fields and original two-orientation mask) also predicts fitted ratios to within 1.69e-4. In particular the old diagonal h=.2,u=.3 gives J/P=.998796 vs fitted .998785, whereas restoring mirrors gives .984904. Thus the old diagonal agreement is not a general continuum property; it is specific to the two-orientation directional stencil. See `inherited_current.json`. This is a post-hoc mechanism diagnostic, not a replacement for the original prospective criteria.

## Five-time-unit free coast

At h=.2 L=7 axis u=.3, early and late fitted speeds .29733996 and .29546640; reduction .6301%, passes 1% speed target. Displacement 1.48121, less than many knot widths. H conservation 4.40e-11; maximum global ΔP/E0=.001979, fails .001. Transverse momentum/E0 ≤1.17e-8. Regional R=3 drift .004844 before flux correction, .000900 after; R=5 .001833 before, .001002 after. This is clear outward transport as well as residual lattice/wall effects, not simply global drag. Finite-domain tail remains; one box and duration cannot identify an asymptotic drag coefficient. Full fit speed/P/H=.980932; M/E0=1.023253. Degree unchanged. Long precision gate partially met, not passed. Fine long coast omitted because short primary ratio fails at .15; no additional qualifying claim.

## Weak external force: linear deficit pair

The applied potential V=-f∫x(1−φ0) in normalized units gives F=fQ for a continuum translation, Q=∫(1−φ0), with no boundary term if vacuum at the wall. Discrete translation is not an exact symmetry, so the momentum impulse comparison tests accuracy of this force definition. Force amplitude f chosen for F(0)/E0=.025; T=1.6, initially static, both signs. Fit quadratic full energy-density centroid including free intercept/initial velocity. Total H+V conserved; bare H increases under forcing. Regional flux residuals without source term are not force ledger residuals.

| sign | fitted a | M_force/E0 | impulse mismatch | H+V relative drift |
|---:|---:|---:|---:|---:|
| +1 | .0321051 | .782922 | .030056 | 3.26e-14 |
| −1 | −.0321052 | .782921 | .030056 | 3.25e-14 |

Both fail mass and impulse targets; sign reversal agrees to ~2 ppm, but both are biased. At T=1.6 bare H rises .829%, while P/E0=.03899 corresponds to small-velocity kinetic energy ~.076% E0; Q grows 2.38%, signaling deformation/energy redistribution. Full centroid displacement .03724; fitting regional acceleration gives .012517 (R3), .024734 (R5), .032105 (full). Applying the same *total* force normalization to those regional accelerations gives M/E0≈2.008, 1.016, .783. These are not alternative legitimate mass measurements because the regional force and flux differ; they demonstrate sensitivity to extended tail and distributed source. No claim that the physical mass is .783 E0. A core-localized squared-deficit control is an exploratory execution amendment, not a replacement for the failed pair.

Squared-deficit control, V=-f∫x(1−φ0)^2: modified force independently validated against NumPy (acceleration max error 2.84e-14, H+V derivative zero to printed precision). Same h=.2,L=7,T=1.6, nominal F/E0=.025, positive sign. M/E0=1.056264, a=.0237051, impulse discrepancy .039728, Q growth .588%, H rise .163%, H+V drift 5.76e-14. Regional accelerations .022972 (R3), .023727 (R5), .023705 (full), far closer than linear source. Fit max residual falls from .000435 to .0000293. Thus suppressing the extended tail sharply reduces centroid bias and energy redistribution, but still fails both the mass and impulse criteria; improvement is evidence of coupling sensitivity, not a qualified measurement. No-force T=1.6 control: velocity 1.16e-9, displacement 1.65e-9, P drift/E0 5.05e-9, energy drift 2.50e-14; rest drift negligible. All forced/control endpoints degree −1.

## Disposition and next steps

- PASS: exact mirror-symmetry validation, checkerboard and vacuum controls, transverse momentum all boosts, topology all 10 dynamics cases (three regular values, no inadmissible tetrahedra), total energy including external potential drift ≤4.40e-11. Short boost energy within 1%, directional split greatly reduced.
- FAIL: primary speed/(P/H) tolerance .5% at both resolutions; dynamic mass within 1% only one of five short cases. Long global momentum drift fails while speed loss passes. Fine long coast not run after primary failure; optional third resolution deferred.
- FAIL: independent force mass and input momentum ledger at both signs and both couplings. Core-localized coupling materially improves mass and spatial centroid robustness, but is unqualified; no selective pass.
- Proposed next (PI decision): develop reflection-symmetric higher-order static *and kinetic* operator with analytic bond current, validate as before and target the .5% canonical P vs exact J gate. Alternatively evaluate one finer h≈.125 axial case to confirm the trend; do not claim extrapolation from two h. For force calibration use a core-localized coupling, test force amplitude and duration convergence, distinguish core motion from energy redistribution with topological-density centroid and signed impulse ledger, include complete external-source term in regional flux. Do not redefine P=J simply to satisfy the velocity identity. No precision 26E until force and coast calibrated. Qualitative exploration is separate from a precision claim.

## Reproducibility

All recorded CPU seconds 3080.18=.8556 core-hours; conservative total <.90 including compilation/analysis, well within 2 cap. Linux Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0. One numerical thread, at most four workers, per-process CPU time limits. `src/exgc26b2.py validate`, `static --h .2 --L 7 --maxiter 300` (also .15); `run --h H --L 7 --u U --axis 1 0 0` or `1 1 1`, `--T 5` coast, `--T 1.6 --u 0 --sign ±1` forced pair, sign 0 rest, `validate --power 2` and `--power 2 --sign 1` for core-localized control; default power 1 (absent in original compact JSON). `src/exgc26b2_current.py` and optional argument path to prior B rest states for inherited-current diagnostic; equation and compact data above. `src/exgc26b2_figure.py` independently recomputes recorded ratios, force mass and impulses from compact data and produces figure. Raw rest/end arrays 183 MB retained in local cloud worktree `local_raw/EX-GC26B2`, hashes and sizes in `raw_manifest.json`, not committed or transferred to PI PC; regenerable. Full trajectories not retained. No remote overwrites or closure sign-off; Claude review pending.
