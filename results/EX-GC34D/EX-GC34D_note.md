# EX-GC34D — a lattice-native covariant gravity energy (Regge) gives GR light bending (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "go ahead and close 34A-r, 34B, 34C and proceed with 34D".
Plan `spec/EX-GC34D_plan.md` (f0f9ebf, pushed before the analysis). Code `src/exgc34d.py` → `analysis.json`;
`src/exgc34d_bend.py` → `bend.json` and `D4_*.json`; `src/exgc34d_shift.py` → `shift.json`. Post hoc:
`src/exgc34d_posthoc_diag.py` → `posthoc_diag.json`, `src/exgc34d_posthoc_iso.py` → `posthoc_iso.json`. Our
construction, not LaViolette's Model G.
**PI CLOSED 2026-11-23** (Richard: "you can close 34D and go on to 34E").

## Conclusion

**A Regge gravity energy (lapse × deficit angles) on the cubic lattice selects GR's static weak field exactly, at
every wavelength the lattice holds, not only in the continuum limit.** Each edge shrinks by the mean lapse of its two
ends, up to a lattice relabelling (vertex displacement):

  l_e = l⁰_e (1 − n̄_e) + gauge,  n = lattice Newtonian potential with G = κ/8π.

That is linearised GR's isotropic solution, so γ = 1 exactly. The fit holds to 1e-13 at 24 wavevectors spread over
the whole Brillouin zone (post hoc). 34C's open question, whether a lattice can supply the covariance rule, has a
positive answer: Regge calculus does it. The only lattice effect is in the Newtonian potential itself, a 0.02%
dispersion at kh = 0.05.

**Two registered light-bending criteria failed in the packet test (D4).** The packet bent 2.056 × the lapse-only
amount and the rays 2.062 ×, against a window of 2 ± 3%. In a post hoc rerun on the same geometry in isotropic
coordinates, the ratios were 1.989 (field) and 2.000 (rays). So the excess is a coordinate effect of the finite
early/late readout windows. The solved metric carries a relabelling part, which an asymptotic or straight-ray
readout does not see but a direction measured 10–18 units from the source does. Separately, rays bent 3–3.5% more
than the field in both backgrounds, which is the opposite sign to 34A-r. That gap put the Regge case outside the
registered field-versus-ray tolerance. I haven't traced its cause.

**A moving source needs a shift-like variable (D5, analytic).** Without one, the equations still give GR's
gauge-invariant field for every Fourier mode with k·v ≠ 0. But the gravitomagnetic part must then be carried by the
spatial metric, with amplitude ∝ 1/(k·v). For modes with k ⊥ v it is missing altogether: the momentum constraint is
violated by κρv. So on the shared lattice, the field of a moving body cannot stay bounded and co-moving unless the
structure gains a shift-like link variable. That is the next structural question.

## Registered results

| part | criterion | result | verdict |
|---|---|---|---|
| D1 | Schläfli per tetrahedron ≤ 1e-10 | 6e-30 | PASS |
| D1 | lapse–lapse block of K exactly 0; K Hermitian | 0; 4e-16 | PASS |
| D1 | ‖KG‖/‖K‖ ≤ 1e-9 (incl. zone boundary) | ≤ 3e-16 | PASS |
| D1 | direct Regge sum vs Fourier kernel ≤ 1e-4 | 1.0e-5 (8³ lattice, ε = 1e-3) | PASS |
| D2 P1 | \|γ − 1\| ≤ 0.01 at kh = 0.05, three directions, shrinking with k | 0, 7.2e-5, 1.2e-4; ∝ k² | PASS |
| D2 P2 | k²h³n̂/(κM̂) = −½ within 1% at kh = 0.05 | −0.50010, −0.50005, −0.50004 | PASS |
| D2 P3 | exactly one O(1) gauge-fixed eigenvalue | one at 120; the other seven ∝ k² | as expected |
| D3 P4 | straight-ray ratio within 2% of 2 (Gaussian, knot B1, knot B2; b = 6, 8, 12) | 2.000000000000 in all 9 | PASS (see note below) |
| D4 P5 | packet and ray ratio Regge/lapse-only within 3% of 2 | 2.056 (field), 2.062 (rays) | **FAIL** |
| D4 P6 | \|α_field − α_ray\| ≤ 5% of \|α_ray(lapse)\| = 8.1e-4 | lapse 4.7e-4; Regge 1.06e-3 | lapse PASS, **Regge FAIL** |
| D5 P7a | ω ≠ 0: reduced (no-shift) system solvable, momentum constraint automatic, Riemann = GR | all three hold (symbolic) | as expected |
| D5 P7b | spatial vector part ∝ 1/ω; k ⊥ v breaks it | h_xy = −2κρv_y/(k²v_x); at k ⊥ v the momentum constraint misses by κρv_y | as expected |

**D3 is stronger than registered, and less informative.** For a ray along a lattice edge direction, h_tt = −2n holds
*exactly* for every mode perpendicular to the ray (post hoc: also along the face diagonal (1,1,0)). That's because
the lattice relabelling modes do not stretch an edge whose ends are at the same phase, so the edge's length change is
exactly the conformal one. The 1e-12 agreement is therefore an identity, not a numerical coincidence. The source
shape (Gaussian, B = 1 hedgehog, B = 2 torus) cannot matter.

**The registered D2 readout understated the result.** The O(k²) departure of γ(k) from 1 comes from my continuum
projector P = δ − k̂k̂, which doesn't remove the lattice relabelling modes exactly. Fitting the lattice modes directly
(post hoc) shows the solution is exactly conformal at every k.

## D4 detail (34A-r machinery, F1 massless packet, impact parameter 6)

| background | α field | α rays | ratio field | ratio rays |
|---|---:|---:|---:|---:|
| lapse only (solved n, flat space) | −0.01565 | −0.01612 | 1 | 1 |
| Regge solution as solved (registered) | −0.03217 | −0.03323 | 2.056 | 2.062 |
| same geometry, isotropic coordinates (post hoc) | −0.03111 | −0.03223 | 1.989 | 2.000 |
| 34A-r prescribed GR-like bump, for reference | −0.03001 | −0.02933 | 2.005 | 2.013 |

Energy drift 9e-6 in each run. The post hoc row is not a re-scoring: P5 and P6 stay FAIL as registered. What it
shows is that the excess lives in the coordinate part of the metric, since the geometry is the same. A readout of
direction at finite distance is coordinate-dependent at order n(r) ≈ 0.003 rad, which is the size of the excess. The
plan should have used an asymptotic or straight-ray readout for a solved (not prescribed) metric. That was my
design slip.

## D5 detail (linearised GR, k along x, dust velocity (v_x, v_y, 0), ω = k v_x)

- Full theory (with shift) and reduced theory (h_0i ≡ 0, lapse and h_ij equations only) are both solvable for
  ω ≠ 0. All linearised Riemann components agree, so the reduced theory has the same physics there. The momentum
  constraint holds automatically, by the Bianchi identity and matter conservation. One residual gauge parameter
  remains.
- The cost: the reduced solution has h_xy = −2κρ v_y/(k² v_x). This diverges as k·v → 0, and in real space it is a
  spatial-metric disturbance that builds up along the path. In GR this is the zero-shift gauge, which is legitimate
  but secular. On the shared lattice, where the link weights are what fields move through, accumulating strain is
  not a harmless choice of coordinates.
- Exactly at k ⊥ v the reduced system solves the static equations, but the momentum constraint is off by κρv_y. The
  gravitomagnetic field of that mode is lost.
- Check: the static Newtonian limit gives R₀₁₀₁ = −κρ/2, that is n = −κρ/(2k²).

## Limits

Linear only. Static except for the D5 algebra. The knot sources carry energy only. The geometry lives on the Kuhn
edges (axes, e_i + e_j, body diagonal), and the nine field links are built from the metric those edges imply. The
periodic box removes the mean source, equally in every background. The field-versus-ray offset (3–3.5%, opposite in
sign to 34A-r) is not explained.

## Compute

About 0.7 core-hour of the 2 cap: kernel, checks and Fourier solves 0.05; registered packets and rays 0.4; post hoc
0.25. Linux cloud, Python 3.13.18, numpy 2.5.3, scipy 1.18.1, mpmath 1.3.0, sympy 1.14.0. The 1.9 GB solved field and
the packet end states are in `local_raw/EX-GC34D/` (not committed; regenerated by `src/exgc34d.py` in about 3 min).

## Suggested next (not approved)

**34E: add a shift-like link variable**, a per-link "flow" that fields are dragged along. Ask whether the pure-clock
plus Regge structure then gives the full ADM form, by checking the gravitomagnetic field of a moving source and the
speed of linear gravitational waves on the lattice. A separate option is a lighter consolidation step: write
34A–34D into the white paper as "light bending and the gravity action".
