# EX-GC26E2 — EX-GC26E1 repeated at h = 0.15 (EXPLORATORY)

Claude, 2026-11-21. PI go-ahead (Richard, 2026-11-21): "I sign off on EX-GC26E1. Go ahead and use my PC for 26E2."
and "Use up to 8 cores on my PC." Runs on the PI's PC under the 40 core-hour PC allowance (about 33–35 left). No new
cloud cap. Massless Skyrme reference model, not Model Gc; the construction is ours.

## What changes from 26E1

Only the grid spacing changes: h = 0.15 instead of 0.2. Everything else is identical and frozen before the runs:

| Setting | Value |
|---|---|
| Code | `src/exgc26e1.py` with `--h 0.15` (Windows-safe CPU limit), driver `src/exgc26e2_run.py`, analysis `src/exgc26e1_analysis.py E2` |
| Box | L = 7, same box for both species |
| Field | g = 0.0125, half-cosine ramp over t ∈ [0, 1], T = 5, dt = 0.05h |
| Coupling | e^{2Φ} on the quadratic term only |
| Control constant | f = 2g·E₂(B=1)/Q(B=1), from the h = 0.15 B = 1 state |
| Relaxation | B = 1 hedgehog, 300 iterations; B = 2 rational map z², up to 600 iterations |
| Runs | the same six: energy, control and null for each species |
| Workers | up to 6 in parallel (at most 8 allowed), one thread each |

## Question

How does the B = 2 ÷ B = 1 fall-rate ratio R move from h = 0.2 to h = 0.15?

- **Binding energy gravitates** (expected): R stays near 1 or moves toward 1, along with the lattice prediction
  (2E₂/E)_{B2} ÷ (2E₂/E)_{B1}.
- **Binding energy does not gravitate:** R moves toward 2E(B=1)/E(B=2), which was 1.041 at h = 0.2. That alternative
  value is recomputed at h = 0.15.
- **Control:** it must again follow its predicted (Q/E) ratio. If it does not, the result is not interpretable.

Report all four readouts as in 26E1: baryon and energy centroids, each with t ≥ 1 and t ≥ 2 windows. Two resolutions
still allow only a *candidate* conclusion, with no continuum extrapolation claim.

## Compute

About 3–4 core-hours, estimated by h⁻⁴ scaling from 26E1's 1.04. Wall time on 6 workers is roughly 1–1.5 hours. Raw
states stay on the PC in `local_raw/EX-GC26E2`; compact records and the figure are committed from the PC through a
branch. Do not touch other processes on the PC (`live_trader.py`).
