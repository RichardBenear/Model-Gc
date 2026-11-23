# EX-GC34A-r — corrected rerun of EX-GC34A (EXPLORATORY)

Claude, 2026-11-23. PI decision (Richard, 2026-11-23, decision card "How should we proceed after the 34A readout
flaw?"): "Rerun as 34A-r". The card stated the total would exceed the 2 core-hour cap (about 2.6). Revised total
estimate: about 2.8 core-hours (34A 1.45, plus this rerun about 1.3).

## Why

EX-GC34A's registered readout (linear energy centroid in a periodic box) failed because packet tails wrapped around
the box (`results/EX-GC34A/EX-GC34A_note.md`). **34A stays as recorded: NOT MET.** This rerun is a new registration
with the readout fixed before any data.

**Correction to the 34A note's M1 recomputation.** On a uniform background a linear field's spectrum is conserved, so
the spectral velocity there just restates the dispersion formula. It is not an independent test of propagation. The
M2 spectral numbers (curved backgrounds, where the wave numbers change) are unaffected.

## What changes

Only the position readout. Everything else is identical to `spec/EX-GC34A_plan.md`: fields, backgrounds, boxes,
packets, integrator, windows, tolerances and pass criteria. That includes the same M1 box and the same M2 box.

- **Readout: periodic (circular) energy centroid**, per axis:
  x̄ = lo + (L/2π) arg Σ e(x) exp(2πi(x − lo)/L), unwrapped in time. For a distribution symmetric about its mean this
  is exactly the mean even when the distribution wraps around the box. Implemented as `circ_centroid` in
  `src/exgc34a.py` (series `cc…` in the records).
- **Rays get the identical definition**: the same circular mean over the energy-weighted ray positions with the same
  box (`box=` in `src/exgc34a_rays.py`). So any residual asymmetry bias is the same in the field and in the rays.
- **Scoring:** `src/exgc34ar_eval.py`, the 34A evaluator with this readout. Output `results/EX-GC34Ar/summary.json`.
- **Method check before registration:** flat background, not a registered case. The packet started at x = 3 and
  crossed the box seam. Circular readout 0.94104 vs rays 0.94112 (−0.008%); the old linear readout gave −1.5.

## Separate recomputation (`src/exgc34ar_check.py`)

- **M1:** from each saved final field, the packet position is the peak of the energy density after FFT Gaussian
  smoothing (width 1), with a parabolic sub-cell fit and periodic distance to the known start, x0 = 0. Then
  v ≈ (x_peak − x0)/T. Compared with the fitted circular-centroid velocity. This is a peak estimator, not a mean: it
  checks the readout, not the solver equations.
- **M2:** the late spectral direction, as in 34A's `exgc34a_check.py`, compared with the late-window circular-centroid
  direction.

## Compute

Same 57 runs as 34A (about 1.0 core-h) plus rays and scoring (about 0.3). Results in `results/EX-GC34Ar/`. Final
fields are in `local_raw/EX-GC34Ar/`, with SHA-256 values committed.
