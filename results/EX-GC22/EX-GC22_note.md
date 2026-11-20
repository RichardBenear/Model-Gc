# EX-GC22 — intrinsic drag and the weak-field limit (EXPLORATORY)

PI: Richard Benear. Design, implementation, computation, and analysis: ChatGPT. Our construction, not LaViolette's. Pre-run plan `spec/EX-GC22_plan.md`, local commit 890c21e before computation. Branch `ex-gc22-chatgpt` based on df7e73a. No registered claim or PI closure.

## Answer so far: steady self-wake drag can be suppressed

The prescribed-source arm supports the weak-source hypothesis. At fixed G flux a weaker source can respond with almost the same acceleration per **assigned** inertial mass while its steady self-recoil drag approaches zero. This separates the wake from the externally imposed response in the current equations. It is not yet Newtonian inertia: the explicit head-on absorption drag is disabled, just as in EX-GC19/20. The full particle arm is in progress.

## A. Prescribed moving sink

The Gaussian sink of width 0.8 has total consumption M = 5.32 K/0.2. A constant compensating mean source is used only for the periodic Fourier diagnostic. The field is solved directly in the translating frame, eliminating initialization transients. No chemical deformation, flow cohesion or lattice pinning is involved. For nonzero wavevector k and ω = k_z v,

```
den = c² |k|² − ω² − iω/τJ
Ghat = −(1/τJ − iω) shat / den
Jzhat = i c² k_z shat / den
```

The uniform G mode is 1; the uniform J mode is zero. The force observable is acceleration `2∫s Jz/G / (20∫s)`, the unweighted prescribed velocity is not solved by the equations. Linearizing at small v with G=1 inside the force gives `γ = (2/20M) ∫ |shat|² k_z²/(D |k|⁴) d³k/(2π)³` with `D=c²τJ`. For an infinite Gaussian this is `M / (6·20·π^(3/2)·σ·D)`. Thus a linear-in-source, inverse-D decay rate is expected from the steady equations themselves, without assuming an additional absorption drag. Finite velocity and the denominator G give additional nonlinear terms.

At v = 0.05, side 32, h = 0.5:

| K | D=1 γ | D=8 γ | D=80 γ | D=80 peak depletion | D=80 external acceleration |
|---:|---:|---:|---:|---:|---:|
| 0.2 | 0.0113806 | 0.00125349 | 0.000234769 | 0.00481135 | 0.00200656 |
| 0.05 | 0.00229202 | 0.000282016 | 0.0000350504 | 0.00120284 | 0.00200163 |
| 0.01 | 0.000436367 | 0.0000548248 | 0.00000575718 | 0.000240568 | 0.00200033 |
| 0.001 | 0.0000431721 | 0.00000544748 | 0.000000547574 | 0.0000240568 | 0.00200003 |

The external acceleration uses uniform additional flux 0.02 with G held fixed (linear superposition diagnostic). The K=0.2, D=80 steady gamma is only 1.4% of the EX-GC20 reported transient 0.017; the models are not identical (fixed Gaussian vs deforming Gc core, periodic vs Dirichlet reservoir). At D=80, K=0.2, replacing G with 1 in the force gives γ=0.000108890, whereas the full J/G gives 0.000234769; that nonlinear contribution is not captured by a naive inverse-D extrapolation. Both terms fall with K. At K=0.001 their difference is only 3.13e-9.

Controls: v=0 gives |a|<6e-20, reversal changes the sign at roundoff, Fourier residuals <4e-15 normalized. h=0.25 at side 32 agrees within 2e-7 relative. Increasing side to 64 raises γ by 8.5% (D=1) and 3.4% (D=80), so box-size error is non-negligible but does not change the qualitative result. Compare absolute numbers with that caveat. `fixed.json` contains all 36 main cases and 10 control cases. Independent analytic small-v limit agrees with the G=1 force evaluation on the same finite box at small v; its continuum infinite-box value differs due to excluded low modes.

## B. Full chemical particle

In progress. The initial G field is recomputed for each distinct transport coefficient with a Dirichlet Poisson solve, then X/Y/W and the symmetric recoil flow are relaxed together. The subsequent time evolution uses the original explicit chemistry and face flux evolution. No knot is included (separate issue).
