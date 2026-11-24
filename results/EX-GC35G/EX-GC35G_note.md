# EX-GC35G — projected Y1 scalar benchmark at h=.125

2026-11-24 | EXPLORATORY | Design and runs ChatGPT. [Plan](../../spec/EX-GC35G_plan.md) · [Claude review](../../review/EX-GC35EF_claude_review.md). PI closure pending.

**Pre-reset PASS. Post-reset PASS. Tidal threshold FAIL.** Halving the physical spacing from .25 to .125 reproduced the predicted nearly fourfold reduction. The reconstructed tidal error was **2.29% at t=16 and 2.25% at t=20**. The late readouts level off over this window. The largest sampled error still exceeds the original 2% threshold.

| t | Reconstructed tidal error | TT error | Full synchronous coordinate tidal error | Improvement from h=.25 | Velocity impulse scale |
|---:|---:|---:|---:|---:|---:|
| 4 | 0.577% | 0.589% | 0.600% | 4.00× | 0.885% |
| 8 | 1.21% | 1.22% | 2.44% | 3.99× | 5.77% |
| 12 | 2.04% | 1.72% | 2.79% | 3.98× | 9.10% |
| 16 | 2.29% | 1.54% | 2.78% | 3.98× | 10.2% |
| 20 | 2.25% | 2.25% | 2.51% | 3.98× | 9.81% |

The run uses the existing weighted projector, fixed physical Fourier support, 18³ product quadrature, RK4 at Δt=.05 and perturbation strength 10⁻³. The 96³ label is the physical derivative and translation spacing for L=12. The product quadrature is unchanged from GC35E. The comparison uses the existing projected continuum snapshots. The full coordinate column includes +∂ₜh⁽¹⁾·∂ₜh⁽¹⁾/4. The impulse uses the per-step velocity correction divided by twice the step and by the continuum tidal norm. The result does not isolate dispersion phase error from other spatial discretization terms.

The peak pre-reset constraint residual was 4.83e−4, 3.86× smaller than at h=.25. The post-reset peak was 2.76e−16. The free control reached 9.10e−2. Cumulative unscaled metric and velocity correction norms were 0.02329 and 0.83593. The projected and free reconstructed TT fields on the same grid agreed to 4.27e−11 relative at the readouts. The t=20 free tidal error was 95.5% against the free continuum control, and the fine projected nominal acceleration error was 10.1%. The nominal and reconstructed quantities differ because the reset changes the accepted trajectory.

The production run and comparison used 180 CPU seconds. [Raw evolution](finer.json) · [Comparison](comparison.json) · [Recipe](../../docs/EX-GC35_projected_Y1_recipe.md)
