# EX-GC35E — constraint projection of the reciprocal Y1 scalar run

2026-11-24 | EXPLORATORY | ChatGPT | [Registered plan](../../spec/EX-GC35E_plan.md) · [Claude review](../../review/EX-GC35CD_claude_review.md). PI closure pending.

**Post-projection PASS. Pre-projection FAIL. Tidal accuracy FAIL.** Projection enforces the sourced constraints at the accepted time samples. This is a separate result from the size of the corrections and the accuracy of the projected trajectory.

## Projection and controls

The homogeneous treatment, physical grids, mode support, first-order evolution, scalar response and next-order source are from GC35C. The correction is applied only to the second-order metric and velocity at accepted RK4 steps. Per mode, write the affine residual as Cx−b. The change is the weighted minimum-norm solution of CΔx=−(Cx−b), using the tensor Frobenius norm and measuring velocity in units of |K|. This weights off-diagonal tensor entries twice and leaves the TT subspace untouched. The unprojected second metric is evolved in parallel from the same initial fields and RK stage sources.

A pilot that weighted the six stored components equally introduced spurious TT changes. The production projection uses the tensor norm above. Its full complex projector has a relative TT leakage bound of 7.68e-16 over the tested modes. A production check compared reconstructed projected and free TT curvature on the same fine grid, with a maximum relative difference of 6.47e-11 at the registered readouts. A band-limited 18³ quadrature evaluates the products while translations and spatial derivatives use the physical coarse and fine spacings. The largest relative deviation from the original free GC35C residuals at t=4 and 20 across the two grids and continuum was 2.55e-7.

| Run | Peak before projection | Peak after projection | Peak free residual | Sum of metric correction norms | Sum of velocity correction norms |
|---|---:|---:|---:|---:|---:|
| coarse | 6.42e-3 | 2.31e-16 | 1.56e+0 | 3.70e-1 | 1.26e+1 |
| fine | 1.86e-3 | 2.57e-16 | 3.68e-1 | 9.30e-2 | 3.30e+0 |
| continuum | 6.89e-10 | 2.62e-16 | 1.57e-7 | 2.47e-6 | 1.33e-6 |
| half | 1.53e-5 | 5.96e-19 | 5.70e-4 | 5.09e-3 | 2.33e-2 |

The residuals multiply by the perturbative strength 10⁻³ and divide by the original constraint-source norm. The correction columns sum unscaled coefficient norms and do not measure a net shift. Spatial refinement reduced the peak pre-projection defect by 3.45×. Halving the coarse timestep changed the pre-projection t=4 value from 3.01e-5 to 1.53e-5. A smaller pre-projection residual from a shorter reset interval is expected even when the underlying source has the same conservation defect.

## Tidal accuracy

Tidal comparison uses all six spatial components with the tensor norm. The registered observable is the GC35A metric contribution −∂²ₜ h⁽²⁾/2. The reconstructed estimate takes second differences of the accepted metric positions and incorporates position and velocity resets. It uses a centered stencil at interior readouts and a second-order backward stencil at the final time. The nominal acceleration uses the unchanged spatial RHS evaluated on the projected state. The impulse column gives the norm of the velocity correction divided by the step and by twice the accepted continuum tidal norm. It is a scale diagnostic. The full coordinate tidal check adds the term +∂ₜh⁽¹⁾·∂ₜh⁽¹⁾/4 in synchronous coordinates and uses the same nonzero Fourier modes.

| t | Reconstructed vs projected continuum | Nominal vs nominal continuum | Full coordinate tidal error | Reconstructed TT error | Free vs free continuum | Velocity impulse scale |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 2.31% | 4.42% | 2.40% | 2.34% | 5.27% | 3.48% |
| 8 | 4.82% | 23.5% | 9.75% | 4.87% | 32.4% | 22.8% |
| 12 | 8.11% | 37.1% | 11.1% | 6.85% | 120% | 36.0% |
| 16 | 9.12% | 41.6% | 11.1% | 6.26% | 260% | 40.3% |
| 20 | 8.94% | 39.9% | 9.97% | 8.81% | 398% | 38.7% |

The table uses the finer lattice. The continuum projection changed its own reconstructed tidal tensor by at most 0.000298%. At t=4 the coarse timestep-halving cross-check changed the nominal result by 0.000181% and the reconstructed result by 0.515% relative. The half-step endpoint uses a backward estimate, while the main run uses a centered estimate at t=4. This is a finite-timestep sensitivity check, with different stencils at that readout. At t=20 fine reconstructed and full coordinate errors were 8.94% and 9.97%.

The weighted projection is a numerical modification of the accepted evolution. The post-projection roundoff residual does not verify an exact conservation identity of the supplied source. The pre-projection and tidal gates test whether the altered evolution is accurate enough for the benchmark.

Production and the quadratic replay used 0.302 core-hours. The initial weighting pilot used 0.540 core-hours. [Driver](../../src/exgc35e.py) · [Curvature replay](../../src/exgc35e_curvature.py) · [Analysis](../../src/exgc35e_analysis.py) · [Tidal results and gates](tidal.json). Run the driver with `coarse`, `fine --n 48`, `continuum --continuum`, and `half --dt .025 --end 4`, then run the curvature replay and analysis. The raw result snapshots for each driver configuration contain the per-step residuals and correction norms.
