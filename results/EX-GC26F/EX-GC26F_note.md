# EX-GC26F — two knots attracting through their self-generated field (EXPLORATORY)

Claude, 2026-11-21. PI go-ahead (Richard, 2026-11-21): "go ahead with 26F". **PI CLOSED 2026-11-22** ("I sign off: 26F and close the series").

- **Plan and code:** plan `spec/EX-GC26F_plan.md`, committed with the LOG row before any PC run (d51adae). Code
  `src/exgc26f.py`, driver `src/exgc26f_run.py`.
- **Records:** from the PI's PC (1b8ad8a), with `summary.json` from the driver.
- **Separate recomputation and figure:** `src/exgc26f_figure.py` → `recompute.json` and
  `paper/figures/exgc26f_attraction.png`.

Massless Skyrme reference model with a dynamical Nordström scalar, not Model Gc. The construction is ours.

## Conclusion

**Two knots attract through the fields they generate themselves, close to Newton's law, and the attraction falls off
as 1/d². By the registered criterion the result is SUPPORTED, but only just, and the pull is systematically 4–10% weak.**

- **Strength.** By the registered measure (baryon centroid, quadratic fit over the whole run), the pull is
  **0.901 of the prediction at d = 6 and 0.960 at d = 9.2**. The window was 0.90–1.10, so d = 6 passes by 0.001.
- **The shortfall is real.** Every readout falls below 1, and none is above 0.96. The separate recomputation (no fit)
  agrees: baryon 0.90 and 0.95. The energy centroid reads lower, at 0.84–0.91 at d = 6 and 0.89–0.92 at d = 9.2.
- **Inverse square.** The ratio of pulls, d = 6 ÷ d = 9.2, is 2.20 (baryon) and 2.38 (energy). The prediction is
  2.35; a 1/d force would give 1.53. 1/d² is clearly favoured, and 1/d is excluded.
- **Bookkeeping is clean.** Total energy is conserved to 2e-10. Total momentum stays ≤ 1.2e-8 of E, and the two
  knots move exactly symmetrically. The degree stays −2.

## Results (prediction: relative acceleration −2 (2E₂/E) A / d², 2E₂/E = 0.936, A = 0.011098)

| d | readout | measured a_rel | predicted | ratio, fit [0, 10] | ratio, fit [1, 10] | ratio, no-fit recomputation (t ≥ 5) |
|---:|---|---:|---:|---:|---:|---:|
| 6 | **baryon centroid** | −5.20e-4 | −5.77e-4 | **0.901** | 0.898 | 0.898 |
| 6 | energy centroid | −5.21e-4 | −5.72e-4 | 0.912 | 0.960 | 0.842 |
| 9.2 | **baryon centroid** | −2.36e-4 | −2.46e-4 | **0.960** | 0.958 | 0.949 |
| 9.2 | energy centroid | −2.19e-4 | −2.45e-4 | 0.894 | 0.903 | 0.909 |

**Gravitational displacement by t = 10** (separation change, coupled − control):

| d | measured | predicted |
|---:|---:|---:|
| 6 | −0.0259 | −0.0289 |
| 9.2 | −0.0117 | −0.0123 |

**How to read the controls.**

- Without gravity the separation wanders by ±0.01. This comes from the direct Skyrme push and the product-ansatz
  knots relaxing and breathing (right panel of the figure).
- The coupled-minus-control difference is nonetheless a smooth parabola (left panel).
- The driver's "a_rel_control" values of ±2e-4 to 5e-4 are quadratic fits to that wandering. They are not a
  constant force.

## Why the pull is weak: candidates (not tested)

1. **Gravity changes the direct Skyrme force.** Each knot sits at Φ ≈ −0.014 (W ≈ 0.97), which alters the
   repulsion between the knots slightly. The subtraction does not remove that change.
   - At d = 6 the Skyrme force is about twice gravity, against about 1.5 times at d = 9.2. The deficit is also larger at
     d = 6 (10% against 4–5%), which fits.
   - The plan estimated 3–6%.
2. **Centroids on the lattice move about 2% slower than momentum ÷ energy** (26B/26D: speed ÷ (P/H) = 0.98).
3. **The knots breathe in the product ansatz.** The passive charge 2E₂/E was measured for a settled knot; it was not
   re-measured for a breathing one.

A test of candidate 1 would repeat d = 6 at a second G_N, or at a larger d. Both are PI decisions.

## What this does and does not show

- **The chain is now complete for a weak scalar field.** Each knot sources a 1/r field (26C). The other knot falls in
  that field at its energy-set rate (26E). The two act together as mutual attraction with the 1/d² law, and without the
  drag Model Gc had (26D). Active and passive mass agree to within 4–10%, and the deficit is likely an artefact of the
  method, not a breach of the equivalence principle.
- **Not a precision test.** One resolution (h = 0.2), weak coupling (G_N = 0.01), and two separations. The direct
  Skyrme force is larger than gravity and is removed only by the control.
- **Scalar (Nordström) gravity limits remain**: no light bending.

## Compute

- **PC:** 23,108 cpu-s = **6.42 core-hours**. Four runs of 3601–8240 cpu-s, 2 h 19 min of wall time on 4 cores.
  - This is over the planned ~5 core-hours: the cloud timings underestimated the PC.
  - Python 3.13.2, NumPy 2.2.4, SciPy 1.16.1, Windows 11.
  - PC allowance left: about 21–23 of 40 core-hours.
- **Cloud:** about 0.3 core-hours (static pre-check and smoke tests).

## Next (PI decision)

- Write up the 26A–F chain: the white-paper revision deferred until EX-GC26 had results.
- Or test the deficit (candidate 1), or move to 26E3 (internal motion).
