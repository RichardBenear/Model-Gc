# EX-GC34H — the cube-mean kinetic term B does not hold up under constraint reduction (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "go ahead with the B run first."
- The plan `spec/EX-GC34H_plan.md` and the `explore/LOG.md` entry were pushed before any run (a07101c). A pre-run
  edit to the zero-cluster definition is bee2875.
- Code: `src/exgc34h.py`.
- Registered outputs: `reduced` → `reduced.json` (H1–H4, H6, C1, S1) and `evolve` → `evolve.json` (H5).
- Post hoc output: `posthoc` → `reduced_posthoc.json`.

This is our construction, not LaViolette's Model G. It follows priority 1 of ChatGPT's review
(`review/EX-GC34F-34G_chatgpt_review.md`).

## Conclusion

**B is not a viable kinetic term as registered. It fails H2, H4, H5 and H6.** H1, H3 (where scored), C1 and S1
(121 of 122) pass. The one edge combination with no kinetic energy is confirmed. Eliminating it does not leave a clean
two-wave system, for three reasons:
- **Drift sector.** The zero-frequency sector keeps a long Jordan chain. Post hoc, it contains constraint-satisfying
  solutions that grow polynomially in time. They are made mostly of the non-metric part of the 7th edge, the
  body-diagonal mismatch.
- **Curvature.** The metric part of the zero sector carries curvature of order (kh)². That is the same size as the
  lattice error of an exact displacement gauge mode where that error is non-zero. But it is not zero on the
  (1,−1,0)-type diagonals, where the displacement control is exact.
- **Wave speed.** In a subset of directions one polarisation is several percent fast at kh = 0.05: 5.8% along
  rand14, ∝ (kh)². This confirms 34G.

## Registered results

| criterion | result |
|---|---|
| H1 one non-kinetic edge combination | **PASS** (249 of 249 cases: one singular value below 1e-12, the next about 0.25) |
| H2 reduction regular (2×2 block σ ratio ≥ 1e-10) | **FAIL** (127 of 249; see below) |
| H3 no growth, two low physical waves (on the 122 reduced cases) | **PASS** (122 of 122); internal residual ≤ 1e-10 PASS |
| H4 zero sector pure gauge (r ≤ max(10 r_disp, 1e-10)) | **FAIL** (12 of 122) |
| H5 reduced time evolution | **FAIL** (both directions) |
| H6 wave speed within 1e-3 at kh ≤ 0.05 | **FAIL** (8 of the 34 reduced cases) |
| C1 control: T growth reproduced | **PASS** (Im ω/k = 1.42211, relative difference 4e-11) |
| S1 roots match 34G | 121 of 122 (rand08 at kh = 0.2: 2.3e-6, threshold 1e-6) |

**H1.** The non-kinetic combination exists in every case.
- Its lapse coupling |C_n,e u|/|C_n,e| is 0.2–0.35 in generic directions.
- On the axes it is about 1e-16.

**H2 failed for two different reasons:**
- **On the axes, a real degeneracy (19 cases).** The lapse does not couple to u at all, and the full pencil is
  singular there (34G found normal-rank deficit 1). That is the exact re-slicing symmetry on the axes.
- **Off-axis at small kh, a threshold I chose badly (108 cases).** The block's singular-value ratio scales as
  (kh)⁴ purely from the units of n and a, even though the coupling is O(0.3). The cases pass from kh = 0.05 up. This
  is a method flaw, recorded as a FAIL as registered.

**H4.**
- The zero-sector curvature measure is r ≈ 1e-4 at kh = 0.05, growing ∝ (kh)². That compares with r ≈ 1 for a
  wave.
- r is about 2–3 × r_disp where r_disp is non-zero, so those cases pass.
- It fails on the (1,−1,0)-type diagonals and three other cases. There r_disp is exactly zero, but r stays at about
  1e-4.

**H5.** In both directions the seeded runs differ from the wave run by 4e-5 to 5e-4 of the tidal amplitude (limit
1e-6). The momentum residual is 2.5e-8 (limit 1e-8). Halving the step changed the result in opposite directions in
the two runs. The reduced F has a ± pair of eigenvalues of about 1e-4 k² (generic) or 7e-4 i k² (110) inside the
zero cluster, which is the signature of a split Jordan block. The evolution therefore mixes numerical splitting with
secular drift. It does not show clean wave propagation.

**H6.**
- The (1,1,0)-type diagonals are 1.7e-3 fast in one polarisation, and rand21 is 1.1e-3 fast.
- Most random directions are 1e-4 to 1e-3 off.
- The worst directions (rand14, rand15, rand29) were among the H2 failures, so they appear only post hoc.

**My stated expectations:**
- H1 passes: held.
- H2 passes except possibly on the axes: wrong, because of the threshold scaling.
- H3 passes: held where scored.
- H4: no firm prediction. My lean towards "gauge up to lattice error" holds only partly.
- H5 fails unless H2–H4 pass: as expected, it failed.
- H6 fails: held.

## Post hoc (labelled; not scored)

For this pass:
- the block threshold is removed;
- on the axes, n comes from the twice-differentiated lapse constraint;
- the zero cluster is |λ| ≤ 1e-2 k²;
- directions of the singular pencil are excluded;
- the zero sector is filtered on all 11 equations.

1. **The fast polarisation is real.** It matches 34G's confirmed roots, with residual about 1e-14. Along rand14 the
   errors are:

   | kh | speed error |
   |---:|---:|
   | 0.0125 | 0.0046 |
   | 0.025 | 0.0155 |
   | 0.05 | 0.058 |
   | 0.1 | 0.22 |
   | 0.2 | 0.73 |

   That is about 23 (kh)². With all directions included, 23 cases at kh ≤ 0.05 exceed 1e-3. The error does not track
   the lapse–u coupling (correlation −0.2).
2. **The zero sector has a long Jordan chain, with secular growth.** In generic directions, three of the five
   constraint-satisfying zero-frequency solutions:
   - are dominated by the body-diagonal mismatch (the non-metric part of the 7th edge), at 0.15–2 of the metric norm;
   - grow by a factor of 400–30 000 between t = 0 and t = 2/k. They grow polynomially, not exponentially.

   Their own curvature measure is only about 1e-4, so the growing part is mostly lattice gauge plus non-metric edge
   content. Even so, an undamped polynomial drift of a lattice variable is not acceptable in a dynamical model. It
   also explains why every time-domain look at B, in 34G and here, was inconclusive.
3. **At kh ≤ 0.025 the reduction runs out of double precision.** Off-axis, the reduced F has a norm of about 3e4 to
   7e4 k². Its defective zero block then splits into ± pairs larger than 1e-2 k². That gives 66 spurious "growing"
   roots post hoc: symmetric ±λ, growing relative to k as kh decreases. They are numerical, but they mean small-kh
   statements about B rest on 34G's argument-principle spectra, not on this reduction.
4. **On the axes, the extra null direction is lapse-dominated** and present at every ω. That is the exact re-slicing.
   The two waves there match 34G to 1e-12.

## What this means

- Neither 34G's fitted terms nor the cube-mean B gives a clean linear lattice graviton. B avoids exponential growth,
  but it pays with:
  - a non-dynamical 7th-edge variable whose constraint chain allows polynomially growing non-metric solutions;
  - direction-dependent wave speeds of several percent at kh = 0.05.
- This matches ChatGPT's priority 2: stop tuning fits. The two routes left are:
  - (a) remove the 7th (body-diagonal) edge from the dynamics, i.e. a variable set with exactly 6 metric components
    per vertex;
  - (b) derive a kinetic term from an exact lattice re-slicing identity, checked for mode count before any longer
    run.
- Neither is approved.

## Compute

About 12 CPU-seconds in total (cap 2 core-hours):
- registered reduced run: 3.4 s;
- registered evolution: 5.2 s;
- post hoc: 3.3 s.

Platform: Linux cloud, Python 3.13, numpy 2.5, scipy 1.18.

## Deviations

- **Precision bug.** The first `reduced` run left mpmath at default precision, which corrupted the 34D
  finite-difference Hessian: entries came out around 1e8. I fixed it by setting `mp.dps = 50`, as 34G does, before
  the recorded run. The faulty run was not kept.
- **H2 threshold.** It was scale-dependent, as described above. It is recorded as a FAIL as registered.
- **Post hoc changes.** Everything in the post hoc section was done after seeing the registered results.
