# EX-GC35K — bending by a matter-sourced Y1 lapse and spatial metric

2026-11-24 | EXPLORATORY | Design and runs ChatGPT | [Plan](../../spec/EX-GC35K_plan.md) · [Claude review](../../review/EX-GC35GHI_claude_review.md). PI closure pending.

**The full-metric bending was 2.009 times the lapse-only bend on the finest packet grid.** The finest group-angle result was **3.135°** toward the mass. The matched continuum finite-path prediction was 2.834° and the leading asymptotic 4GM/b prediction was 2.903°. The factor near two is reproducible in the test. A precision match to the absolute continuum angle is still open.

## Matter source and event rule

Solve the static linear Y1 equations for a Gaussian three-dimensional density with GM=.095, width 1.4, and zero mean mode on a 160³ periodic grid of spacing 1. The Fourier solution uses the Y1 derivative symbol, the lapse perturbation n=−κρ/(2|K|²) on nonzero modes, and the isotropic spatial solution γᵢⱼ=(1−2n)δᵢⱼ. Five axis and off-axis mode probes gave a relative full-kernel residual of 1.07e−16. Sample the equatorial metric N=1+n and spatial scale A=sqrt(1−2n). The lapse-only control keeps the same N and sets A=1.

Propagate a localized complex scalar packet in the 2D equatorial metric with impact parameter b=7.5 and carrier 4 in units with c=1. In these coordinates the spatial link lengths are scaled by A. In two spatial dimensions their conformal factors cancel in the scalar face flux and leave A⁻² in the proper-time equation, φ̈=(N A²)⁻¹ div(N grad φ), where the second time derivative is in local proper time. The site events use equal proper increments and the lapse sets their coordinate spacing. Local periods use 2,048 clock quanta with a maximum lapse rounding error of 2.44e−4. At h=.12 the sites generated 2,364–2,498 events through t=50. The flat control generated 2,500 each.

A clock interface reads fixed neighboring versions. The first run used the GC35I linear past predictor and produced large amplitude and angle errors in the spatially varying clock field. For example, at h=.16 the lapse-only field norm grew to 2.93× its initial value and the measured group angle was 6.90°. The follow-up uses a quadratic approximation. When the upper bracketing version precedes the output event, it can be used with two earlier versions. Otherwise the update extrapolates from three earlier versions, with a linear startup fallback. The chronological event queue executes from the fixed dependencies. Reversing the ordering of equal-time events gave bitwise-identical output for the h=.16 full case. This is a limited scheduling check, distinct from a competing-rewrite test.

## Deflection and numerical controls

The table lists magnitudes toward the lens. The group-angle estimate uses the forward packet spectrum and the lattice group direction. The phase angle uses the spectral wavevectors. The continuum finite-path estimate integrates the transverse gradient of log(A/N) along straight rays, averages over the packet waist, and uses the measured longitudinal span in the same periodic metric.

| Packet spacing | Step parameter | Clock and geometry | Group bend | Phase bend | Finite-path continuum estimate |
|---:|---:|---|---:|---:|---:|
| .16 | .04 | Uniform clock, lapse-only | 1.634° | 1.527° | 1.435° |
| .16 | .04 | Uniform clock, full | 3.283° | 3.069° | 2.828° |
| .16 | .04 | Varied clock, lapse-only | 1.617° | 1.511° | 1.435° |
| .16 | .04 | Varied clock, full | 3.263° | 3.050° | 2.828° |
| .16 | .02 | Varied clock, lapse-only | 1.632° | 1.525° | 1.434° |
| .16 | .02 | Varied clock, full | 3.279° | 3.066° | 2.827° |
| .12 | .02 | Varied clock, lapse-only | 1.560° | 1.502° | 1.437° |
| .12 | .02 | Varied clock, full | 3.135° | 3.019° | 2.834° |

The uniform-clock control advances in common coordinate steps and uses the same metric in the wave operator. At h=.16, halving the varied-clock step made its angle agree with this control to about 0.004° and brought the full field norm ratio from 0.974 to 0.997. The refined flat direction was within 1e−6 degrees of zero. The full group-angle excess over the matched finite-path estimate fell from 16.0% to 10.6% between h=.16 and .12. The leading carrier dispersion of the centered spatial stencil predicts a group-angle multiplier sec²(kh/2) for slight bends. At h=.12 this gives 3.00° for the full finite-path estimate, compared with 3.135° measured. Finite-width propagation, ray-path corrections, finite spacing, and the periodic source are limitations of the comparison.

The metric is solved from the static density externally. The packet neither generates nor updates the density. The event rule therefore tests propagation with matter-sourced static clocks and spatial geometry in the equatorial section. Reciprocal metric dynamics and causal equivalence under contested inputs remain open.

Combined recorded CPU time was 825 seconds. [Final driver](../../src/exgc35k.py) · [First predictor](linear_predictor.json) · [Timestep and common-clock controls](controls.json) · [Finest packet results](refinement.json). Reproduce the last three runs with `python src/exgc35k.py TAG --h .12 --step .02 --mode flat|lapse|full`.
