# EX-GC26A2 — independent numerical qualification (EXPLORATORY)

**Summary.** A2 reproduces the original degree-one Skyrmion and resolves the apparent fractional topological degree and lopsided-start discrepancies. The geometric degree is exactly −1 for all resolved states tested. Most apparent static translation energy variation is from the outer boundary and interpolation rather than a large periodic lattice barrier. All six short directional boosts meet the A2 momentum threshold, but leakage remains direction-, resolution- and box-dependent; a larger box halves it. A new higher-order static discretization meets the Derrick ratio at h=.15 but is not fully converged, and its dynamics have not been implemented. A2 therefore supports exploratory 26B qualification with the inherited solver, not a precision inertia claim or a retrospective pass of A. The default 2 core-hour cap was respected. PI decides closure and next stage; Claude reviews independently.

PI: Richard Benear; design/implementation/runs ChatGPT; review pending Claude. Protocol committed prospectively a0b5c58 from affe4ba; `spec/EX-GC26A2_plan.md`. No chemistry or G: standard massless SU(2) Skyrme model, our numerical construction, not LaViolette's. Source `src/exgc26a2_qualify.py`, `src/exgc26a2_stencil.cpp`; figure `src/exgc26a2_figure.py`. Reference `src/exgc26a_skyrmion.py` remains unchanged. No new physical mass measurement.

## Validation and independent implementation

For the inherited all-forward/all-backward energy and exact derivatives, independent C++ vs Claude NumPy on a random small field agrees to 3.89e-16 in normalized energy and 2.78e-17 in gradient. Tangent finite difference relative errors 5.71e-8 (original), 1.05e-8 (higher order). Accelerated independent C++ evaluation of the same constrained acceleration agrees to 2.49e-14 against NumPy, then used with inherited RK4 and projection. The higher-order stencil assigns positive energy to a checkerboard test. This is separate-author cross-validation of the discretization, not independent physical replication. Validation JSON, all settings and outcomes retained. L-BFGS uses normalized per-site vectors and radial penalty (zero for unit norm), exact chain rule; terminal projected gradients and iteration cap recorded, rather than assert formal optimizer convergence.

## Static reproduction, topology and lopsided start

| h, half-box | original E/12π² | FD degree magnitude | geometric degree | E2/E4 |
|---|---:|---:|---:|---:|
| .25, 5 | 1.183846 | .875134 | −1 | .775520 |
| .2, 5 | 1.206077 | .928981 | −1 | .860966 |
| .15, 5 | 1.221179 | .962493 | −1 | .911933 |
| .2, 7 | 1.202362 | .930736 | −1 | .879643 |

Reference .2 E=1.20607713873529 vs 1.20607714173089 here (2.48e-9 relative), original FD degree agrees 2.51e-6 absolute. The L-BFGS optimizer hit the 300/350 iteration cap, but target energy/degree independently reproduced and projected gradient small; no claim of exact stationary solution. The noisy original seed, run 700 L-BFGS iterations, E=1.20608767 (8.73e-6 relative above symmetric), FD B=−.928990 and geometric B=−1. Meets A2 lopsided agreement target even though optimizer reached cap. The original A C3 remains historical not met. h=.15 grid is even-sized (68 points, actual half box 5.025), so node alignment differs from .2; geometric tests and whole-cell shifts explicitly recorded.

Geometric method: Freudenthal triangulation into six oriented tetrahedra per cell. Map each by affine interpolation of its four unit S³ vertices followed by normalization. A sufficient hemisphere certificate (all vertices have positive dot product with normalized vertex sum) holds for every tetrahedron in the tested static and final dynamic fields. For three fixed regular values near the antipode, solve the 4×4 vertex matrix against target direction, count strictly positive barycentric solutions weighted by orientation and determinant sign. Degenerate tetrahedra (det below 1e-13) skipped; regular values away from vacuum and all three counts agree; vacuum control returns 0 and parity reflection +1. This measures degree of an admissible interpolated map, not an exact invariant of every possible lattice update. All reported degree-one fields give −1 (orientation convention); final dynamics also −1. The original FD estimator remains noninteger; do not retrospectively change A's thresholds. For background on geometric S³ tetrahedra see B. Svetitsky, A. G. Schramm, `Topology and metastability in the lattice Skyrme model`, Phys. Rev. D 62, 114020 (2000), https://arxiv.org/abs/hep-lat/0008003. Our regular-value counting differs from their volume calculation.

## Pinning first: apparent versus periodic effect

Frozen relaxed fields translated by cubic interpolation, normalized and fixed at box faces, yield relative energy ranges 2.27e-4, 2.21e-4, 1.83e-4 at h=.25,.2,.15 L=5; .2 L=7 gives 6.04e-5. Integer-cell translations (which should have no periodic pinning difference away from boundaries) themselves change energy by ~1.6–2e-4 at L=5. Independent radial profile sampling (from Claude's radial ODE, sparse profile interpolation; not re-relaxed) confirms substantial whole-cell variation caused by clamping nonzero tail at the faces. Larger boxes reduce it. Thus the raw L=5 pinning test fails its 1e-4 target; do not treat its range as pure lattice barrier.

Post-hoc discriminator: smoothly taper the relaxed field to exact vacuum between radii L−2 and L−1, then periodically translate via Fourier phase, project to S³ and repeat with cubic interpolation. Taper alters the tail so this is not a rigorous pinning bound for the original fields. However, whole-cell shifts now preserve energy at roundoff and the residual is very small:

| h, L | original Fourier energy range/E | original cubic energy range/E |
|---|---:|---:|
| .25,5 | 6.70e-8 | 1.11e-4 |
| .2,5 | 4.12e-8 | 5.72e-5 |
| .15,5 | 1.37e-7 | 2.04e-5 |
| .2,7 | 2.35e-9 | 1.82e-5 |

The Fourier phase-shift result is not monotonic in h (taper and grid alignment differ). It establishes that the particular smooth compact profiles have negligible periodic barriers compared with raw figures; the original discrepancy is not evidence of order-1e-4 lattice pinning. Cubic interpolation substantially overestimates the range. Grid pinning alone does not explain the observed momentum drift. The stronger evidence for a boundary contribution is the dynamic comparison below. Source actual per-offset energies retained, so any alternative curvature analysis can be performed; wall clamping and interpolator contamination prevent interpreting a raw curvature as a physical force.

## Higher-order static candidate

New operator averages eight independent forward/back orientations of third-order biased 4-point derivative; odd-order bulk truncation terms cancel (formally fourth order on smooth interior fields). Fixed vacuum ghost extension includes a 3-cell halo. Unlike the inherited two-orientation choice, independent sign averaging gives axis-reflection symmetry and avoids central checkerboard nulls. It changes the discrete action, not the continuum model.

| h, L | E/12π² | E2/E4 | FD degree magnitude | geometric degree |
|---|---:|---:|---:|---:|
| .25,5 | 1.291293 | 1.188246 | .927908 | −1 |
| .2,5 | 1.265936 | 1.085025 | .947769 | −1 |
| .15,5 | 1.249420 | 1.001139 | .967227 | −1 |
| .2,7 | 1.259815 | 1.092926 | .948480 | −1 |

At h=.15 the Derrick ratio is within .114% (original 8.81%), and E differs from the radial reference 1.23136 by 1.47%, meeting numerical point targets narrowly. But all runs reached 250-iteration limit and normalized projected force max is 4.25e-5 at h=.15, 4.83e-4 at h=.2; do not assume full convergence. Hard vacuum boundary and nonzero radial tail produce particularly large artificial energy for this high-order stencil. Compact Fourier pinning ranges ~0.9–1.7e-6 (larger than original but still below target); raw pinning at L=5 fails and even L=7 slightly exceeds target. No kinetic operator/boost run for new stencil. A2 higher-order result is promising but **unqualified for inertia**, not a reason to replace inherited solver wholesale.

## Momentum qualification: original conservative dynamics only

Initialize v=−0.1 (unit direction · ∇)φ tangent to S³, not exact Lorentz boost; monitor energy centroid, canonical momentum (same observable as A), FD degree and geometric degree at T=1.5. No energy/inertia fit. All cases completed, all geometric final degrees −1.

| h, L | direction | dt | max ΔP norm/E0 | max ΔP norm/P0 norm |
|---|---|---:|---:|---:|
| .2,5 | axis | .010 | 4.185e-4 | .004211 |
| .2,5 | diagonal | .010 | 3.661e-4 | .003736 |
| .15,5 | axis | .0075 | 3.707e-4 | .003710 |
| .15,5 | diagonal | .0075 | 3.434e-4 | .003461 |
| .2,5 | axis | .005 | 4.185e-4 | .004211 |
| .2,7 | axis | .010 | 2.075e-4 | .002093 |

All meet A2 1e-3 E0 short-window target, worst energy drift 6.18e-13 relative. Halving dt changes momentum drift by 3.1e-13 E0, so temporal integration not the driver. Refining h from .2 to .15 modestly reduces leak (axis 11.4%, diagonal 6.2%), enlarging box at .2 halves it (50.4%). This suggests finite-domain/tail effects are significant even before a wavefront could cross the box from its center; massless Skyrmion tail extends everywhere. No boundary-flux ledger here, so not a proof of origin. At .2 initial axial P has transverse components ~−.000842 (each ~.7% of axial); diagonal direction has equal components. A finite-grid anisotropy is explicit, consistent with two-orientation stencil not individual-reflection invariant. Centroid velocities ~.0973 at .2, not exactly imposed tangent speed .1. This is also why a direct initial P/.1 ratio is not an independent inertial mass measurement. A2 is shorter than A C6's 3-unit gate, no retroactive C6 pass.

## Disposition

- PASS: independent original numerical reproduction, geometric degree, lopsided energy/degree, short clean directional A2 momentum gate and timestep check.
- CONDITIONAL: compact-profile pinning and boundary control; raw frozen-field translations at L5 outside tolerance, not resolved into a pure physical pinning force.
- PARTIAL: higher-order static candidate passes h=.15 Derrick point but not broad convergence or dynamics qualification; original C2 remains not met.
- Recommend 26B as **exploratory gated**: start with inherited h=.2,L=7 and h=.15,L=7 (new settled rest), static translational kinetic tensor and stress/virial ledger (axis, diagonal), then small ± impulses and centroid response at independently measured velocities; require convergence of mass ratio and full momentum including flux. Longer clean coasts and box refinement before physical conclusions. Do not conflate boosted solution of a known relativistic action with independently invented inertia. Keep higher-order dynamic extension separate until static convergence and boundary treatment are established. No gravity or orbits.

## Compute / reproducibility

All recorded case CPU totals 1243.8 s = .3455 core-hours; compilation/validation and report overhead <.02 core-hours, conservative total <.37 core-hours, well within 2 core-hour cap. Linux, Python 3.12.14/NumPy 2.3.5/SciPy 1.17.0. Numerical library threads=1, parallel independent processes max 4 except brief overlap of short diagnostics. State files retained locally in `local_raw/EX-GC26A2` (not committed), SHA-256 in each static JSON; regenerate via source, deterministic seeds. Run `python src/exgc26a2_qualify.py validate`, then `static --h H --L L --mode 0 --maxiter 300` (mode 0 L7 maxiter 350; mode1 maxiter 250 after corresponding mode0 state), `lopsided --maxiter 700`, `radial`, `compact`, `boost --axis 1 0 0` / `1 1 1` with recorded settings. Figure script consumes compact JSON. Local branch `ex-gc26a2-chatgpt`; no remote write attempted. PI CLOSED 2026-11-26.
