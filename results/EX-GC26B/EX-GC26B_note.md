# EX-GC26B — numerical inertial mass calibration (EXPLORATORY)

**Summary.** The inherited conservative Skyrme solver reproduces continuum relativistic inertia to within 1% at h=.15, L=7, u=.1 in both tested directions, but not across the full velocity range: axial errors rise to 1.08% at .3c and 2.43% at .5c. All 15 cases preserve topological degree and conserve discrete energy extremely well. Refining space roughly halves mass errors, while enlarging the box or halving dt has almost no effect on the axial mass calibration. Regional stress-flux residuals improve under refinement but the fastest axial run misses the target. This places the principal numerical limitation in spatial discretization and velocity initialization/measurement rather than time integration. A low-speed operating point is promising; the full precision gate is not met. No new physical inertia or gravity claim. Review and PI closure pending.

PI: Richard Benear. Design, implementation and runs: ChatGPT; Claude review pending. Base 7f7d086, pre-run plan published to main as d4c98a108ea1746f3105774d6cf5e009622f3cc8. Plain massless Skyrme model, no chemistry or gravity. Lorentz invariant continuum action predicts M=E0/c²; numerical calibration, not a new physical derivation. See `spec/EX-GC26B_plan.md` and response to Claude in `review/EX-GC26A2_response.md`.

## Method and validation

Retain the validated inherited all-forward/all-backward stencil and A2 C++ acceleration, RK4 projection. In units c=1, normalized energies and momenta divided by 12π². On rest data, integrate kinetic tensor for v = −u_i ∂_i φ with centered tangent gradients using the original kinetic metric. Also compute integrated continuum spatial stress with matched one-sided derivatives. For the matched prescription, M^match = E I + ∫T, a useful algebraic check, not a discrete conserved Noether momentum theorem. The actual rest fields have residual stress due to the grid and box.

For launch speed u and unit direction n, sample rest φ at x′=x+(γ−1)n(n·x), with φ_t=−γu n·∇′φ; cubic interpolation, normalized tangent projection, clamp at walls. This is the continuum exact Lorentz boost prescription applied to an imperfect finite rest lattice; no claim of exact discrete boost symmetry. The derivative is computed from the same interpolant via symmetric offset ε=1e-4. H and canonical P match the inherited Hamiltonian to 0 and 5.20e-18 respectively on a random tangent field. Boost derivative vs finite-time difference at interior sites agrees to 8.59e-9; vacuum zero. Matched static stress identity error 1.39e-17. A2 source-state SHA256 exactly matches the committed A2 record.

For dynamics, calculate actual velocity from linear regression of full energy-density centroid over T≈.8. M/P estimate is mean P_parallel/(γ(v)v), with measured v, not imposed u. Report H/(γ(u)E0) and P/(γ(u)uE0) as independent checks against launch data. Also report H/(γ(v)E0), centroid fits for R≈3,5, and geometric degree pre/post. Regional numerical momentum ledger: integrate outward T_ij on cube faces by trapezoidal quadrature and compare P_R(t)−P_R(0)+∫F_R dt. Spatial derivative and surface quadrature are continuum-stress approximations on the lattice; not an exact discrete ledger. No flux through the exterior clamp is asserted. Short time only.

## Static and controls

| h | nominal L | E0 | E2/E4 | centered tensor eigenvalue min/E0 | max/E0 | matched tensor max eigenvalue/E0 |
|---:|---:|---:|---:|---:|---:|---:|
| .2 | 7 | 1.202362 | .879643 | 1.001329 | 1.007624 | 1.042868 |
| .15 | 7 | 1.217099 | .931332 | 1.001472 | 1.004803 | 1.023182 |
| .2 | 9 | 1.201201 | .885559 | 1.000408 | 1.006619 | 1.041550 |

All degree −1 and hemisphere certificate held. The two newly optimized states reached the iteration cap, not formal optimizer convergence (projected gradient max 1.80e-5 and 5.00e-6). Centered tensor agrees within 1%, but matched stress-based diagnostic does not: apparent agreement partly depends on choice of translation derivative. Do not interpret favorable centered values alone as proof of calibration.

The rest control h=.2 L=7 has fitted speed −1.14e-9, global P drift 3.76e-9 E0 and regional residual ≤3.33e-10 E0. At .3c axial, halving dt from .01 to .005 changes M/E0 by 1.31e-10, while changing L=7→9 shifts M/E0 from 1.021540 to 1.021339; spatial lattice error dominates this mass estimate. Contrast A2's lower-speed tangent kick, where L=5→7 halved momentum drift. Neither finding proves the precise origin of all leakage.

## Dynamic cohort (L=7, dt=.05h)

| h | u/c | direction | measured v/c | M/E0 | H/(γ(u)E0) | transverse P/E0 max | global ΔP/E0 max | max regional residual/E0 |
|---:|---:|---|---:|---:|---:|---:|---:|---:|
| .2 | .1 | axis | .100029 | 1.014734 | 1.000271 | .001019 | 5.52e-5 | 1.69e-5 |
| .2 | .1 | diagonal | .100075 | 1.000394 | 1.000291 | 4.95e-10 | 4.39e-5 | 5.34e-6 |
| .2 | .3 | axis | .297437 | 1.021540 | 1.002493 | .003673 | 6.39e-4 | 5.60e-4 |
| .2 | .3 | diagonal | .298811 | 1.004325 | 1.002712 | 7.84e-10 | 2.34e-4 | 1.28e-4 |
| .2 | .5 | axis | .483584 | 1.048755 | 1.007275 | .009136 | .004828 | .004828 |
| .2 | .5 | diagonal | .492110 | 1.018001 | 1.007777 | 1.12e-9 | .001010 | .000890 |
| .15 | .1 | axis | .100114 | 1.007337 | 1.000145 | .000507 | 5.06e-5 | 3.70e-6 |
| .15 | .1 | diagonal | .100135 | 1.000195 | 1.000154 | 5.69e-10 | 4.80e-5 | 8.32e-7 |
| .15 | .3 | axis | .299045 | 1.010754 | 1.001358 | .001811 | .000226 | .000132 |
| .15 | .3 | diagonal | .299679 | 1.002193 | 1.001473 | 9.21e-10 | .000154 | 2.96e-5 |
| .15 | .5 | axis | .492447 | 1.024256 | 1.003981 | .004421 | .001124 | .001124 |
| .15 | .5 | diagonal | .496447 | 1.009044 | 1.004209 | 1.68e-9 | .000363 | .000216 |

All 12 main cases plus rest, box and dt controls completed, geometric degree exactly −1 at both endpoints (all three targets, no inadmissible tetrahedra), worst relative energy drift 3.66e-11. Actual finest L=6.975, T=.8025; cubes R=3.075 and 5.025. H/(γ(v)E0) also within 1% for six finest cases. Full centroid includes radiation; restricting to R≈3 reduces measured speed by about 1–1.6%, illustrating sensitivity of mass estimates to definition and tail exclusion. No continuum extrapolation from only two resolutions.

The flux correction is small compared with residual in the fast axial case; the corrected residual still exceeds 1e-3 E0 at .5c. At h=.15 u=.3 axial, regional residual falls from 5.60e-4 to 1.32e-4 E0 relative to h=.2, substantially better. At the finest h and .1c both regions meet the tolerance even without correcting flux. At .1c some flux-corrected values exceed raw values (both tiny); no claim that flux correction universally improves the ledger. Diagonal momentum symmetry reflects the lattice's permutation symmetry, not independent proof of continuum isotropy.

## Disposition against prospectively specified targets

- PASS: centered static kinetic tensor within 1%, with stress-matched caveat; H/(γE0) within 1% on finest h; energy conservation and endpoint topology all 15 cases; spatial refinement improves mass, high-speed regional residuals and transverse P.
- PARTIAL: dynamic mass within 1% at h=.15 in all diagonal cases and the axial .1 case (4/6); axial .3 and .5 fail. Axial vs diagonal mass difference <1% at .1 and .3, fails .5 (1.52 percentage points). Transverse P target passes diagonal and low-speed axial; fails .3 and .5 axial. Measured speed differs from launch by up to 1.51% at finest axial .5, passing 1% in other fine cases.
- PARTIAL: all fine regional residuals improve or are below 1e-4, but axial .5 misses 1e-3 threshold (1.124e-3). This is not an exact lattice momentum accounting identity.
- Overall full precision gate NOT MET. Recommend *conditional exploratory* low-speed continuation only at h=.15,L≥7, speeds ≤.1c, with directional checks, independent force/acceleration response and longer coast before gravitational mass comparison. A targeted 26B2 could use a third resolution or an all-eight-orientation reflection-symmetric stencil (with consistent kinetic operator and exact dynamics, validated anew), not simply substitute the A2 higher-order static operator. Do not proceed to precision 26E on the current evidence alone. PI assigns next stage.

## Reproducibility and cost

Recorded CPU seconds 2002.99 = .5564 core-hours; validation, compilation, reporting and recomputation overhead conservatively total <.02 additional core-hours, well below 2. Linux, Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0; one BLAS thread per process, max four workers. `src/exgc26b_mass.py validate`, `static --h .2 --L 7 --source <A2 state>` (when absent regenerate via A2), `static --h .15 --L 7 --maxiter 350`, `static --h .2 --L 9 --maxiter 350`, `boost --h H --L L --u U --axis 1 0 0` or `1 1 1`, `--dtfac .025` for timestep control, `--u 0` rest test. `src/exgc26b_figure.py` independently recomputes fitted mass and launch energy from compact time series and checks against original summary; figure `paper/figures/exgc26b_mass.png`. All raw rest/end fields (337 MB) retained in local cloud worktree `local_raw/EX-GC26B` with manifest and SHA256; not committed, not yet copied to PI PC, regenerable. End fields only, not full dynamics/velocity checkpoints. Review source, compact records and code in repo. PI closure pending.
