# EX-GC34A-r — one shared link structure: one light cone, and textbook light bending (EXPLORATORY)

Claude, 2026-11-23. PI decision (Richard, decision card): "Rerun as 34A-r". Plan `spec/EX-GC34Ar_plan.md` (7a4f44a,
pushed before the rerun); setup and criteria identical to `spec/EX-GC34A_plan.md`, only the position readout changed
(periodic/circular energy centroid, applied identically to fields and rays). Scoring `src/exgc34ar_eval.py` →
`summary.json`; recomputation `src/exgc34ar_check.py` → `check.json`; figure `src/exgc34a_figure.py` →
`paper/figures/exgc34ar_cone.png`. Our construction, not LaViolette's Model G. EX-GC34A itself stays NOT MET.
**PI CLOSED 2026-11-23** (Richard: "go ahead and close 34A-r, 34B, 34C and proceed with 34D").

## Conclusion

**When every field has to move through the same lattice links, every field sees the same light cone. On prescribed
backgrounds the solver reproduces the textbook light-bending ratios.** 68 of 69 registered cases pass. The one miss is
the nonlinear sigma model in the GR-like background, which overshoots the ray prediction by 5.3% of the tolerance
reference, against an allowed 5%.

- **One cone (M1, 54/54).** Three backgrounds (flat, anisotropic with off-diagonal metric, the same metric with a
  different lattice split) and three directions. Massless, massive, two-component (both modes) and sigma-model fields
  all move at the shared-structure prediction: within 0.19% (linear) and 0.32% (sigma), with direction within
  8e-4 rad.
- **The control is caught.** The field with its own link weights is +4.7% to +4.9% off the shared prediction and
  within 0.2% of its own. A per-field speed is exactly what this structure removes, and the readout sees it.
- **Light bending (M2, 14/15).** Deflection of massless waves past a prescribed weak field Φ = −0.05/√(r² + 4):

| background | field (all massless types) | rays | ratio to lapse-only |
|---|---:|---:|---:|
| Nordström (conformal) | −0.00007 | 0 | 0.005 |
| lapse only (PPN γ = 0) | −0.01497 to −0.01502 | −0.01457 | 1 |
| GR-like (γ = 1) | −0.03001 to −0.03010 | −0.02933 | 2.005 |

  Massive waves bend more, as massive rays predict. In the conformal background they bend while light does not, which
  is the classic Nordström signature.
- **The fields agree with each other more closely than with the rays.** Massless types differ by ≤ 3e-4 relative
  among themselves. All of them bend 2–3% more than the ray prediction, the same in every case. That common offset is
  the eikonal approximation behind the rays, not a field disagreement. It is what tipped F4 in the GR-like background
  just over the line (|Δα| = 7.7e-4 against a tolerance of 7.3e-4). Registered result: **F4/B3 FAILS**. Since F4
  agrees with F1 to 3e-4, I don't read it as a nonlinear-field effect.

## What this does and doesn't show

- **The common cone here is structural, by construction.** That was the point: no field carries its own speed
  coefficient, so nothing can be tuned. The test verifies the code, the metric-to-link mapping and the light-bending
  readout. It is not evidence that nature works this way.
- **The backgrounds were prescribed.** Whether matter *produces* a lapse plus spatial metric with γ = 1 (GR-like
  bending) or only a scalar (Nordström, zero bending) is the question for 34B and 34C.
- **The lattice split is a lattice artefact, and it is predicted.** B1 and B1′ (same metric, different split of the
  nine weights) differ by up to 1% in packet speed at kh = 0.6. Each matches its own ray prediction to within 0.3%.

## Separate recomputation (`check.json`) — weaker than registered precision

- **M1 peak estimator.** It agrees with the centroid velocities to 1.7% (median), worst 2.3%, every field in the same
  direction (slow). The peak estimator assumes motion from t = 0 at one speed, so it carries the start-up transient
  and the packet's dispersive skew. It confirms the readout to about 2%, not to the registered 0.3%. The estimator is
  the weak point, and I'd use a better one in 34B.
- **M2 late spectral direction vs late circular-centroid direction differ by about 20%.** This comparison was
  registered but is apples to oranges. For a lensed, spread packet, the circular centroid of an asymmetric
  distribution is not its mean. That is why the plan applies the same definition to the rays. The meaningful
  comparison is spectral vs ray (linear mean) late direction, from 34A's post hoc diagnostic: −0.01262 vs −0.01238
  (B2), −0.02535 vs −0.02483 (B3) and +0.00003 vs 0 (B4). These agree within 2%.

## Compute and platform

Runs 1.21 core-h plus rays, scoring and recomputation about 0.35. 34A-r about 1.55; with 34A (about 1.45), about
3.0 core-hours in total, against the 2.8 stated at the decision. Linux cloud container (4 cores), Python 3.13.18,
numpy 2.5.3, scipy 1.18.1. Final fields are in `local_raw/EX-GC34Ar/` (not committed), with SHA-256 values in
`manifest_sha256.json`.

## Next

34B: put the knot onto the shared links (including the Skyrme term, which needs pairs of links), give the lapse and
link weights their simplest local dynamics, and test GC33A's source identity. A plan needs PI go-ahead. Lesson for 34B
and 34C: use wrap-immune readouts from the start, and compare like with like (ray readout = field readout).
