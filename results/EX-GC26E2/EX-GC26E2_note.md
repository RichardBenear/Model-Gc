# EX-GC26E2 — EX-GC26E1 repeated at h = 0.15 (EXPLORATORY)

Claude, 2026-11-21. PI go-ahead (Richard, 2026-11-21): "Go ahead and use my PC for 26E2", "Use up to 8 cores". **PI CLOSED 2026-11-21.**
Plan `spec/EX-GC26E2_plan.md`, pushed before the run (3ec3738). The run was on the PI's PC: 6 workers, records in
commit 9e5eb3c, merged to main. Code unchanged from 26E1: `src/exgc26e1.py --h 0.15`, driver `src/exgc26e2_run.py`,
analysis `src/exgc26e1_analysis.py E2`. Figure `paper/figures/exgc26e2_fall.png`.

This uses the massless Skyrme reference model, not Model Gc. The construction is ours.

## Conclusion

On the finer grid, the bound pair and the single knot still fall together. The B = 2 ÷ B = 1 ratio is 0.988–1.009
across the four readouts, against a lattice prediction of 1.006 and a continuum value of 1.

The alternative "binding energy does not gravitate" predicts 2E(B=1)/E(B=2) = 1.043. All four readouts are now 3.4%
to 5.5% below it, more than the 2% spread between readouts. At h = 0.2 the energy-centroid readout (1.020) could not
separate the two cases; at h = 0.15 it has moved to 1.009, toward 1 and toward the lattice prediction.

So, as a candidate result from two resolutions: **binding energy gravitates and falls with the body**, in this
energy-coupled scalar theory. The non-energy control still shows a large, correctly predicted violation (1.32–1.41,
against a predicted 1.39).

## h = 0.2 → 0.15

| | h = 0.2 | h = 0.15 |
|---|---:|---:|
| Derrick ratio E₂/E₄, B = 1 / B = 2 | 0.880 / 0.902 | 0.932 / 0.942 |
| Binding energy fraction | 3.9% | 4.1% |
| Lattice-predicted a/g, B = 1 / B = 2 (energy coupling) | 0.936 / 0.948 | 0.965 / 0.970 |
| Measured a/g, B = 1 (baryon, t ≥ 1) | 0.920 | 0.976 |
| Measured ÷ prediction, B = 1 / B = 2 (baryon, t ≥ 1) | 0.983 / 0.966 | 1.011 / 0.994 |
| **Ratio B2/B1, energy coupling: baryon t ≥ 1 / t ≥ 2** | 0.995 / 1.003 | **0.988 / 0.998** |
| **Ratio B2/B1, energy coupling: energy centroid t ≥ 1 / t ≥ 2** | 1.020 / 1.020 | **1.009 / 1.009** |
| Lattice-predicted ratio | 1.013 | 1.006 |
| No-binding alternative 2E₁/E₂ | 1.041 | 1.043 |
| Control ratio: baryon / energy centroid (predicted) | 1.345 / 1.445 (1.425) | 1.316–1.325 / 1.407–1.409 (1.390) |
| Momentum gain ÷ impulse (all runs) | 0.964–0.969 | 0.981–0.983 |

Two further points:

- The momentum transfer matches 26B3 at both resolutions (0.965 and 0.982), so it is a common lattice effect.
- The predicted fall rates move toward the continuum value of 1 as the Derrick imbalance shrinks. The measured single
  rates follow them.

## Checks

- **Separate recomputation.** I refitted the accelerations from the raw time series with my own least-squares code,
  not the analysis script. All eight a/g values and all eight ratios match the summary to 4 digits.
- **Validation on the PC.** The PC kernel agrees exactly with the 26B2 kernel. It agrees with the independent NumPy
  sum to 8e-16 (energy), 3e-17 (gradient) and 4e-14 (acceleration). The ramp ledger mismatch is 5e-5 and the energy
  rate is 6e-11.
- **Topology and conservation.** Every run keeps its degree (−1 or −2) at both ends. Energy after the ramp is
  conserved to ≤ 2e-12.
- **Null runs.** Drift is ≤ 6e-9 in position and ≤ 1.5e-9 E in momentum.

## Deviations

- **Compiler.** The PC has no `c++` on PATH, so the kernels were built with MSVC cl 19.44.35211 (x64), flags
  `/O2 /std:c++17 /LD` with explicit exports. The cloud build used g++ -O3 -march=native. No source was edited. PI chose
  this option ("A").
- **Compute.** 4.39 core-hours, above the plan's 3–4 estimate. Platform: Python 3.13.2, NumPy 2.2.4, SciPy 1.16.1,
  Windows 11. The PC allowance has about 29–31 of 40 core-hours left.

## Limits

- **Two resolutions.** No continuum extrapolation is claimed.
- **Test body only.** The field is imposed: the knots do not source it (26C), and there are no two-body forces (26F).
- **One fall direction**, along the axis.
- **Fit-window sensitivity.** The baryon centroid still rings after the ramp, so its a/g shifts by 1–2% between fit
  windows. The B2/B1 ratio is far less sensitive.

## Next

26E3: an isospinning B = 1 and a breathing B = 1. Both carry internal kinetic energy; the theorem applies to them only
on time average. ChatGPT flagged that radiation from the spinning state must be checked first. Alternatively, move to
26C (the knot sourcing its own field). PI decision.
