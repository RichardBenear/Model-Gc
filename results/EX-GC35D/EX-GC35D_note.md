# EX-GC35D — Skyrmion relaxation and static Y1 stress compatibility

2026-11-24 | EXPLORATORY | ChatGPT | [Registered plan](../../spec/EX-GC35D_plan.md) · [Claude review](../../review/EX-GC35AB_claude_review.md). PI closure pending.

**Full-source compatibility failed.** All eight cases reached the tighter projected force tolerance. Refinement substantially reduced the mismatch, while tighter optimization and larger walls changed little at h=.2. The finest measured residuals were 0.867% and 0.698%, above the 0.1% gate for moving-knot runs.

| B | h | Box half width | Sites | Degree integral | Max. projected force density | Full residual | Stress divergence / energy-source norm | Tr ∫T / E |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.25 | 6.00 | 49³ | -0.87738 | 8.19e-6 | 2.96% | 0.07908 | 11.8% |
| 2 | 0.25 | 6.00 | 49³ | -1.7957 | 8.78e-6 | 2.24% | 0.05658 | 9.02% |
| 1 | 0.2 | 6.00 | 61³ | -0.93012 | 9.60e-6 | 1.65% | 0.04420 | 6.76% |
| 2 | 0.2 | 6.00 | 61³ | -1.8786 | 9.31e-6 | 1.31% | 0.03300 | 5.50% |
| 1 | 0.15 | 6.00 | 81³ | -0.96306 | 8.60e-6 | 0.867% | 0.02324 | 3.95% |
| 2 | 0.15 | 6.00 | 81³ | -1.9348 | 7.69e-6 | 0.698% | 0.01758 | 3.31% |
| 1 | 0.2 | 7.60 | 77³ | -0.93094 | 8.52e-6 | 1.64% | 0.04403 | 6.25% |
| 2 | 0.2 | 7.60 | 77³ | -1.8800 | 8.93e-6 | 1.30% | 0.03269 | 5.04% |

The normalized force divides the gradient of the energy in Skyrme units by h³ and removes the component normal to the local unit sphere. The stopping tolerance was 10⁻⁵, with a per-case compute limit of 650 CPU seconds. Models use fixed-vacuum walls. The Fourier solve uses the nonzero modes with |k|≤2.05, matching the integer mode set for the width-6 runs and the physical cutoff across widths. The periodic FFT length includes the final grid point, so the actual periods across spacings differ slightly. The energy-only residual was below 4.9e-16 in every case. The direct stress-divergence check uses the phase-aligned Y1 lattice derivative. The full integrated stress tensors and separate quadratic and quartic energies are archived with each case.

At h=.2 the raw GC35B stopping gradients divided by h³ were 2.65e-4 and 3.10e-5. Remeasuring those 350-iteration states on the same mode set here gave 1.6539% and 1.3096%. Tightening the force tolerance changed the mismatch by only 0.00330% and 0.00249% relative. Increasing the half width from 6 to 7.6 changed it by 1.01% and 1.08% relative. The endpoint orders under spatial refinement were 2.41 and 2.28. With three spacings and nonzero boundary errors, these are observed trends rather than verified continuum orders.

The stress trace equals E₄−E₂ for this stencil at each sampling point before integration. Its residual decreases under refinement, as does the divergence of the placed stress. An independent metric differentiation audit reproduced the quadratic and quartic stresses to roundoff and reproduced the plaquette averaging and phase for a Fourier wave. Those controls establish the intended stress convention and placement. They leave open whether a different action-level matter/geometry coupling would improve the finite-spacing Ward identity. The mismatch at the available resolutions limits use of the full Skyrme source. It cannot be explained by incomplete h=.2 relaxation or the tested wall shift alone.

Recorded CPU use for relaxations and matched baselines was 0.803 core-hours, plus final source measurement and short audits, within the two-core-hour cap. The moving-knot gate was not met. See the [driver](../../src/exgc35d.py), [baseline driver](../../src/exgc35d_baseline.py), and [stress audit](stress_audit.json). Regenerate the main cases with `python src/exgc35d.py TAG H B [--half 7.5]`.
