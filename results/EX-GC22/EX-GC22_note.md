# EX-GC22 — intrinsic drag and the weak-field limit (EXPLORATORY)

PI: Richard Benear. Design, implementation, computation, and analysis: ChatGPT. Our construction, not LaViolette's. Pre-run plan `spec/EX-GC22_plan.md`, local commit 890c21e before computation. Branch `ex-gc22-chatgpt` based on df7e73a. No registered claim or PI closure.

## Answer: self-wake drag can be suppressed, but the full flow model has a translation defect

The prescribed-source arm supports the weak-source hypothesis. At fixed G flux a weaker source can respond with almost the same acceleration per **assigned** inertial mass while its steady self-recoil drag approaches zero. In full Gc, weakening K also reduces decay of the internal flow by 117–187×, but increasing D by 80× at fixed K does not reliably reduce it. The strongest wave particle has almost the same flow decay as the diffusive particle even though its peak depletion is 61× smaller. The new instantaneous balance identifies a specific confound: the flow V is not transported with the moving chemical pattern, so averaging over the shifting consumption field produces apparent velocity decay even with zero integrated recoil. Applying the missing uniform-transport term to that balance cancels 99.9% of the initial decay in the strong wave case. This is an analytic diagnostic of a candidate solution, not a dynamic test of the correction. The actual chemical centroid still slows on the coarse grid. It is not yet Newtonian inertia: the explicit head-on absorption drag is disabled, just as in EX-GC19/20.

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

## B. Full chemical particle (8/8 main cases complete)

The initial G field is recomputed for each distinct transport coefficient with a Dirichlet Poisson solve, then X/Y/W and the symmetric recoil flow are relaxed together for 30 time units. The final five-unit change in M is at most 0.116%; the discrete G equilibrium residual is ≤6e-11. Wave J starts at −D∇G, unlike the original mismatched EX-GC20 initialization. The subsequent time evolution uses the original explicit chemistry and face flux evolution. No knot is included (separate issue). Bare drag off, κ=1, flow response τ=20, h=0.5, box ±6, initial V boost +0.05, 20 time units. Decay is −slope(log Vcore) over t=10–20; early fit t=2–10 also recorded.

| K | D | Initial G depletion | Late γ | V(20)/V(0) | M(20)/M(0) | Late centroid speed / mean flow |
|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | 1 | 0.351568 | 0.00997815 | 0.747950 | 0.940404 | 0.381671 |
| 0.05 | 1 | 0.123003 | 0.00355558 | 0.918541 | 0.960860 | 0.497719 |
| 0.01 | 1 | 0.0275250 | 0.000826086 | 0.982249 | 0.967399 | 0.563317 |
| 0.001 | 1 | 0.00282906 | 0.0000855199 | 0.998189 | 0.968750 | 0.579274 |
| 0.2 | 80 | 0.00572520 | 0.0104327 | 0.791262 | 0.948024 | 0.587184 |
| 0.05 | 80 | 0.00169473 | 0.00252819 | 0.944106 | 0.957352 | 0.528510 |
| 0.01 | 80 | 0.000351502 | 0.000547163 | 0.988437 | 0.966851 | 0.568520 |
| 0.001 | 80 | 0.0000354407 | 0.0000558554 | 0.998833 | 0.968699 | 0.579778 |

1. Weak K robustly suppresses the measured flow decay. In the weak regime γ is approximately linear in K, not simply proportional to depletion across D values. At K=0.2, D=80 and 1 give nearly equal late γ despite the much flatter wave G field. At K=0.001, the full decay is ~102× the steady prescribed-sink γ at D=80, compared with ~2× at D=1; the profiles and domains differ, so this is a diagnostic contrast, not an exact subtraction.
2. All cores survive; M remains within 6% and W peak within 2% of their respective settled values. At K=0.001 the D=80 internal velocity loses only 0.117% in 20 units, but the W centroid moves 0.663 rather than 1.0, and its speed falls from 0.0392 (t=0–5) to 0.0290 (t=10–20). Thus stable flow is NOT established particle coasting. Coarse-grid pinning/slip is plausible from EX-GC14/15, but the present experiment does not identify it conclusively. Fine-grid full dynamics remain unresolved.
3. Flux response survives. Four additional trials have positive imposed Jz=0.02, E=−0.02/D, correct initial flux; four matched no-kick nulls. Early linear Vcore accelerations (t≤5):

| K | D=1 | D=80 |
|---:|---:|---:|
| 0.2 | 0.00244117 | 0.00191899 |
| 0.001 | 0.00199051 | 0.00198681 |

The nominal value 2 Jz/τ=0.002 is recovered to 0.7% for both weak variants, under the assigned mass-per-consumption rule. Actual weak centroid acceleration is only ~0.00150 (quadratic fit). Null velocities <3e-19, centroid drift <3e-17. We have not established physical equivalence or phase-space momentum conservation. Reducing K reduces the active source M as well; the external flux is independently supplied, not generated by another equally weak particle.

## C. Unexpected mechanism: an Eulerian velocity field that is not carried with the particle

This is a post-hoc analysis of the eight settled checkpoints, not another production run. EX-GC19/20 advect the chemical deviations with V, but evolve V locally without the convective derivative V·∇V. The reported Vcore = ∫s V / ∫s is a weighted average, not directly a material velocity or a conserved momentum. Its exact derivative includes two terms:

```
dVcore/dt = ∫s ∂t V / ∫s + ∫(∂t s)(V − Vcore) / ∫s
```

Immediately after a uniform boost, the local inward self-flux and radial flow remain symmetric. Net recoil and cohesion integrals are zero to numerical error, but the chemical pattern starts moving through the unadvected radial velocity field. This changes the weight s and generates a negative acceleration in the measured Vcore. `initial_balance.json` evaluates all terms directly from the saved rest states with the same finite-difference derivative; no time fit or free parameter. At K=0.2, D=80 the predicted initial decay is 0.0124056, versus observed early fitted 0.0124617. Almost all of this initial decay arises from reweighting, not a nonzero net opposing force. That is why the EX-GC20 push-off control (which removes the radial flow itself over time) is not a unique diagnosis of accretion drag. At weak K the reweighting term shrinks essentially in proportion to K, independently of D, consistent with the measured trend.

To diagnose a candidate repair, compute just the transport term −v_boost ∂z V in the same initial balance. It cancels reweighting to 99.895% at K=0.2, D=80: initial γ would be 0.0000130716 instead of 0.0124056 (factor ~949), without switching off recoil. At D=1 the cancellation is 91.9% (strong G profile remains fixed at the initial instant). At K=0.001, D=80 the corresponding rates are 0.0000707409 and 3.15e-9. These are instantaneous balance evaluations only; no dynamics with the corrected rule have been run. They isolate a specific transport defect rather than prove a fully physical repair.

## D. Recommended next: EX-GC23 — conservative collective momentum/transport

Before a Lorentz-invariant redesign, test a material rather than Eulerian velocity derivative, ∂t V + V·∇V, and ultimately a conservative momentum density with a consistent mass/energy and exchange ledger. Cohesion must be checked for net momentum conservation rather than assumed. Re-run the matched strong wave kick with and without transport at h=0.25, measure the W centroid as well as weighted V and integrated momentum, evaluate the recoil/cohesion/reweighting/transport terms throughout, and test a larger box to separate waves and walls. Use identical external flux trials to verify that restoring transport did not remove attraction. Then confront the distinct EX-GC17 bare absorption drag (which EX-GC22 has deliberately left disabled), perhaps by a reversible field-mediated exchange instead of material accretion. Neither the weak-source limit nor the initial cancellation proves Lorentz invariance or emergent inertia.

## Reproducibility, deviations and limits

- `src/exgc22_drag.py fixed`, `full --workers 2`, `smoke`; C++ chemistry stencil independently compared with NumPy to <1e-10 on a test array, both transport integrators smoke-tested. Linux C++ compiler required for accelerated path. No inherited solver modified. `python tests/test_regression.py` passes all three radial checks.
- Fixed and full records, `summary.json` derived independently from saved time series by `src/exgc22_summary.py`, `initial_balance.json` from `src/exgc22_initial_balance.py`, `paper/figures/exgc22_drag.png` via `src/exgc22_figure.py`. Initial-balance inputs are the eight `local_raw/EX-GC22/settle_*.npz` checkpoints (included in transfer package, excluded from git, hashes retained). Recreate from code if necessary; no claim of external independent replication.
- First full batch aborted after settling D=1 cases on missing h argument; both branches corrected after smoke test. No scientific results were excluded. The pre-run plan includes the amended initialization cadence and worker-overlap deviation. Measured successful production CPU time 456.9 s (0.127 core-h); including failed attempt bounded conservatively below 0.35 core-h, well within 2 core-hour cap. Platform and versions in JSON. Two worker processes (brief initial overlap with fixed arm as logged).
- Fixed arm includes 36+10 cases; full arm includes all 8 kick cases plus 4 flux and 4 no-kick nulls. The full runs use positive signs only (fixed sign reversal explicitly tested), single deterministic trajectory per case. No fine-grid full run, because of execution time; fine-grid evidence is limited to the fixed arm. Full cases are no-knot and only 20 units, with one τJ=20 relaxation time in wave cases; late fit not necessarily steady. No claim of indefinitely maintained motion. The actual motion/flow mismatch at h=0.5 makes physical coasting unresolved.
- GitHub read access succeeded but `git push` failed due to missing authentication and the authorized connector's separate branch creation returned 403. No remote mutation succeeded. Results and complete local commit chain are provided in a Git bundle for PI/Claude import; no main branch changes. PI decides future ownership and closure.
