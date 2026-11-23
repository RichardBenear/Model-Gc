# EX-GC34F — the 34E growing mode is real, survives the continuum limit, and comes from the per-tetrahedron kinetic term (EXPLORATORY)

**PI CLOSED 2026-11-23** (Richard: "I sign off on 34F and 34G"). Reviewed by ChatGPT 2026-11-23: `review/EX-GC34F-34G_chatgpt_review.md`.

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "run 34F". The plan `spec/EX-GC34F_plan.md` and the
`explore/LOG.md` entry were pushed before any run (920041e). Code: `src/exgc34f.py`.
- Registered outputs: `f1` → `spectrum.json`, `f2` → `evolve.json`, `f3` → `mechanism.json`.
- Post hoc output: `f1p` → `spectrum_posthoc.json`, plus a labelled extension block inside `evolve.json`.

This is our construction, not LaViolette's Model G. Diagnostics only: any fix is a separate study (34G).

## Conclusion

**The instability is genuine, and it does not go away on a finer lattice.**
- Growing mode in the time domain:
  - Started on the growing eigenvector, the constrained time evolution grows at exactly the predicted rate
    (1.42211 k against 1.42211 k).
  - The lapse and momentum constraints stay satisfied to 1e-12.
  - The growth shows in the tidal (curvature) field as well as in the mode amplitude.
- Directions:
  - Every wavevector direction tested has one growing pair, except the lattice axes and the three (1,−1,0)-type face
    diagonals. In those exceptional directions the re-slicing symmetry is exact, so the mode is pure gauge.
  - The growth rate depends strongly on direction: Im ω/k ranges from 0.35 to 15.6 over 30 random directions.
  - At fixed physical wavenumber, the rate tends to a constant as the spacing shrinks (kh from 0.2 down to 0.0125).
- Size of the effect:
  - The mode's curvature-carrying part is small, 0.047–0.049 (kh)² of its size.
  - But it grows fast. A seed 1e-8 the size of a gravitational wave overtakes the wave's tidal field after about
    3 wave periods (generic direction, kh = 0.05).

**Mechanism.**
- What the mode is:
  - It is the would-be re-slicing direction.
  - Off the special directions, the lattice gives it a small kinetic energy and a small curvature (Regge) energy, both
    O((kh)²) and of opposite sign. Their ratio, not their size, sets ω², which is why the growth rate is O(k) at any
    spacing.
  - The lapse–Regge coupling does not set the rate (sensitivity 0).
  - In the zero-shift gauge, ω² is set by the edge kinetic term against the Regge edge term, with sensitivities −1
    and +1.
- Diagnostic variants (for 34G, not fixes):
  - **Fitting the strain rate per cube** (19 edges) instead of per tetrahedron removes the growing mode in every
    direction tested, and keeps the 34D static field and the 34E gravitomagnetism. But it adds a third, slower
    propagating mode (0.73–0.85 c), so it is not GR either.
  - **Changing the DeWitt trace coefficient** (λ = 1/2, 1/3) also removes it. That changes the continuum theory.
  - **Removing the lapse–Regge coupling** removes it but destroys the Newtonian limit.

So the 34D–34E tensor benchmark cannot be kept as it stands. The per-tetrahedron kinetic term is the thing to
replace, and a replacement must keep exactly two modes.

## Registered results

| part | criterion | result | verdict |
|---|---|---|---|
| F1 symmetry | which of 48 cube symmetries map the spectrum onto itself (≤ 1e-10) | group order 1 as run: the root sets included split zero roots and differed in length | **FAIL as registered** (method; see post hoc) |
| F1(a) | growing roots reproduce `review_unstable.json` to 1e-6 in ω/k | 111, generic ≤ 6e-8; 110 and near-100 differ by up to 5e-4 (kh = 0.2) | **FAIL** for 110 and near-100: the review values there were less accurate (see below) |
| F1(b) | axis directions decided | decided: no growing root on the axes (exact re-slicing direction, pure gauge) | PASS |
| F1(c) | every root passes residual ≤ 1e-10 and singular-value contrast ≤ 1e-6 | 1789 failures, mostly genuine roots rejected by the contrast check, plus split zero roots counted as roots | **FAIL as registered** (method) |
| F1(d) | exactly two propagating modes, kh ≤ 0.05, \|c − c_light\| ≤ 1e-3 | two in every direction; max 2.5e-4 | PASS |
| F1(e) | kh = 3 axis root, growing-branch multiplicity (reported) | the kh = 3 root is a lattice cutoff branch (ωh ≈ 2.3–2.45), present on the axes at every kh; one growing pair per direction | reported |
| F2(a) | lapse and momentum residuals ≤ 1e-10 throughout | ≤ 3.3e-12 | PASS |
| F2(b) | run (i) grows at Im ω within 1% in amplitude and E_ij, at both steps; steps agree to O(dt²) | 1.42216/1.42212 k vs 1.42211 (generic); 2.12391/2.12385 vs 2.12384 (110); error falls ×4 per halving | PASS |
| F2(c) | wave amplitude change ≤ 1e-8 | ≤ 5e-13 | PASS |
| F2(d) | seeded component: tidal growth at Im ω within 1%; crossing time t× within 1% of prediction | growth 1.42210 and 2.12384 k: pass. t× (74 and 46 time units) lies beyond the registered run length 15/Im ω | **growth PASS; t× not reached** as registered |
| F3 | a named term set reproduces the growth rate to 1% in 111 and generic | axis-edge gauge, kept blocks kinetic shift–shift, kinetic shift–edge and Regge edge–edge: 111 within 0.13% (kh = 0.05) to 2% (kh = 0.2); generic off by 35% | **FAIL as registered** (generic) |

**Why the F1 method failed, and the post hoc check.** In the zero-shift gauge there is an exact six-fold zero root:
three static displacements, each with a chain e = t·G0 that violates only the momentum constraint. Numerically it
splits into roots of size 1e-6 to 4e-2 k. The registered "zero" threshold (1e-8 k) let these through, and they were
classified as growing or decaying modes. The singular-value contrast also proved meaningless here, because the
kernel's singular values span many orders of magnitude.

Post hoc (`spectrum_posthoc.json`; not a re-scoring), every root was confirmed with the argument principle. The
winding number of a rank-completed determinant around the root must equal the number of roots found inside. That is
a separate computation from the eigensolver, done with two random completions. The zero cluster's multiplicity was
taken from the winding around ω = 0.
- 1470 of 1508 non-zero roots are confirmed, all with residual ≤ 7e-13 in all 11 equations and momentum residual
  ≤ 5e-14.
- The 38 unconfirmed are split zero roots (36), plus two near-axis cases at 1e-3 to 1e-4 rad where the numerical rank
  changes.
- On the confirmed roots, the symmetry group is the 6 axis permutations times inversion (12 elements), matching to
  ≤ 9e-10. The other 36 signed permutations differ at O(1). The Kuhn triangulation's body diagonal makes (1,1,0) and
  (1,−1,0) inequivalent.

**F1(a).** For 110 and near-100, the 34E review values came from minimising a singular value in a gauge with an exact
null direction (residuals there were 1e-9 to 5e-7). The F1 roots have residual ≤ 1e-13, and the F2 time-domain growth
for 110 (2.123839 k) agrees with F1 (2.123836 k) to 1.6e-6, not with the review value (2.12387).

**F2(d), post hoc extension** (runs extended to 40/Im ω; labelled in `evolve.json`). Crossing times:

| direction | measured t× | predicted t× | difference |
|---|---:|---:|---:|
| generic | 73.95 | 73.72 | 0.3% |
| 110 | 46.11 | 46.18 | 0.17% |

At equal mode amplitude, the growing mode's tidal field is 4.1e-12 (generic) and 2.2e-11 (110) of the wave's, at
kh = 0.05. Even so, it dominates after 26 e-folds: k·t× ≈ 18 (generic) and 12 (110).

## Spectrum (post hoc confirmed roots)

| direction family | growing Im ω/k at kh = 0.05 | as kh → 0 | propagating (ω/k) |
|---|---|---|---|
| axes (100, 010, 001) | none (re-slicing exact, pure gauge) | none | 0.99976, 0.99995, plus a cutoff branch ωh ≈ 2.45 |
| (1,1,0)-type | 2.1238 | 2.1215 (kh 0.0125) | 0.99984, 1.00014 |
| (1,−1,0)-type | none (re-slicing exact) | none | 0.99979, 0.99987 |
| (1,1,1) | 1.4155 | 1.4143 | 1.00002 (pair) |
| (1,±1,∓1) and (1,−1,−1) | 2.0018 | — | 0.99983, 0.99989 |
| generic (0.8, 0.5, 0.33) | 1.4221 | 1.4210 | 0.99991, 1.00009 |
| 30 random directions | 0.35–16.4 | 0.35–15.6 | two per direction |

**Approaching the axis** (kh = 0.05):

| path | 1e-1 rad | 1e-2 rad | 1e-3 rad |
|---|---:|---:|---:|
| in the (x, y) plane | 1.64 | 1.52 | numerical rank change |
| off-plane | 1.30 | 1.21 | 1.20 |
| toward (1,1,0) | 1.76 | 2.08 | 2.12 |

So the rate stays O(1)·k right up to the exact directions and is discontinuous there. It is smallest near the
(1,−1,0)-type diagonals: 0.35 at 8.7° from one.

**Larger kh.** Growing modes persist to kh = 1 (110: 13.4; 111: 2.76) and are absent at kh = 2 and 3.

**Static zero modes.** One static zero mode with one Jordan chain appears only on the axes and the (1,−1,0)-type
diagonals. That is the exact re-slicing gauge chain. No unexplained defective root was found.

**Curvature content.** The non-gauge fraction of the growing mode's metric perturbation is 0.049 (kh)² for 111 and
0.047 (kh)² for generic, at every kh from 0.0125 to 0.2.

## F3 mechanism details

**Sensitivities, d ln ω² / d ln(block scale), at kh = 0.05:**

| gauge | kinetic shift–shift | kinetic shift–edge | kinetic edge–edge | Regge lapse–edge | Regge edge–edge |
|---|---:|---:|---:|---:|---:|
| axis-edge | +1.00 | −2.00 | 0.00 | 0.00 | +1.00 |
| zero-shift | (absent) | (absent) | −1.00 | 0.00 | +1.00 |

The sensitivities sum to zero, as they must.

**Energy balance at the root.** In the zero-shift gauge, with |n| = 1, the edge kinetic and Regge edge contributions
are ±1.1e-6 k² at kh = 0.05. Both scale as (kh)² (fitted exponent 1.99), and they cancel. Their ratio fixes ω², and
its sign makes ω imaginary. The lapse–Regge contribution is 1e-18, so the lapse acts only through the constraint.

**Minimal model.** Kept blocks (axis-edge gauge): kinetic shift–shift, kinetic shift–edge, Regge edge–edge.

| direction | error at kh = 0.0125 | kh = 0.05 | kh = 0.2 |
|---|---:|---:|---:|
| 111 | 1e-4 | 1.3e-3 | 2.1e-2 |
| generic | 35% | 35% | 36% |

The minimal model fails for generic because dropping the edge kinetic block entirely is not a small change there,
even though its local sensitivity is zero. In the zero-shift gauge the minimal model loses the root altogether,
because the lapse constraint disappears.

**Diagnostic variants** (none is a fix; registered for 34G only after this):

| variant | growing modes (111, generic, 110, 100; kh 0.05, 0.1) | propagating | 34D static kernel | E4 O(v) error, kh = 0.05 |
|---|---|---|---|---|
| base (34E) | 1.42, 1.42, 2.12; none on 100 | 2 | unchanged | 4.2e-4 |
| (a) per-cube strain fit | none | 2 plus an extra mode at 0.85 c (111) and 0.73 c (generic) | unchanged | 4.9e-4 |
| (b) λ = 1/2 | none | 2 | unchanged | 2.1e-4 |
| (b) λ = 1/3 | none | 2 | unchanged | 2.6e-4 |
| (c) no lapse–Regge coupling | none | 2 | changed (5.4e-2) | 1.0 (fails) |

Re-slicing breaking (34E E2 measure): base and (a) ~1e-9 at kh = 0.05; (b) ~2e-6. In (b) the symmetry is broken at
the continuum level, as expected for λ ≠ 1. Variant (a) also breaks the exact re-slicing on the axes (nk = 0 there).

## What this means

- The 34E kinetic action is linearly unstable in the continuum limit: the growth rate is O(k) and independent of the
  spacing. It is not a lattice artefact that goes away on refinement. A time-stepped or non-linear version would blow
  up from round-off within a few wave periods at the faster directions.
- The cause is the per-tetrahedron strain-rate construction. It breaks the re-slicing symmetry differently in
  different directions. Where the breaking is non-zero, the near-gauge direction acquires small kinetic and curvature
  energies of opposite sign.
- **34G candidates.** Variant (a) shows the instability is not forced by Regge × lapse itself. A different strain-rate
  construction, or one that keeps re-slicing exact or derives the projection from a preserved constraint (per the
  34G criteria in the plan), is the natural next step. Variant (a) as built adds a third mode, so it fails those
  criteria.

## Compute

About 1 CPU-minute in total, far below the 2 core-hour cap:
- registered F1: 10 s;
- post hoc F1: 36 s;
- F2: 5 s;
- F3: 6 s;
- plus development runs.

Platform: Linux cloud, Python 3.13.18, numpy 2.5.3, scipy 1.18.1, mpmath 1.3.0.

## Deviations

- The F1 root classification and validation thresholds failed as registered (see above). They were corrected post hoc
  with an argument-principle check; nothing was re-scored.
- One development bug was fixed before the recorded F2 run: run (i) at "dt/2" first reused the dt step.
- Registered F2(d) could not evaluate t× within the registered run length; the extension is post hoc.
- The plan said the 30 random directions would be mapped to the fundamental domain. Because the registered symmetry
  test returned the identity only, they were sampled on the whole sphere. That is a superset, so no coverage was lost.
