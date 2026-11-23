# EX-GC34G — no tested kinetic term passes; one candidate is free of growth but has new problems (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "Go ahead with 34G". The plan `spec/EX-GC34G_plan.md` and the
`explore/LOG.md` entry were pushed before any run (705dc8f). Code: `src/exgc34g.py`.
- Registered outputs: `spec` → `candidates.json`, `evolve` → `evolve.json`, `check` → `check.json`.
- Post hoc outputs: `evolve_posthoc_B*.json`.

This is our construction, not LaViolette's Model G.

## Conclusion

**None of the seven fixed candidates passes as registered. Every one fails G1 (growth), G2 (mode count) or G3
(Jordan chains).** In more detail:
- **T, the 34E per-tetrahedron control,** fails as expected: growing modes up to 40 k at large kh, and 2.1 k at small
  kh.
- **A, the per-cube least-squares fit,** fails worse than 34F suggested. 34F only tested four directions. Over the
  full set A also grows, for example 1.47 k along (1,1,−1). It also has the third slow mode, at 0.35–0.85 c.
- **C, the per-vertex 14-edge fit,** has no low-frequency growth. But it has a superluminal third mode at 2.0–2.9 c in
  every off-axis direction, and growth at kh = 3.
- **The blends D(s)** have growing modes in every case (0.38 k to 12 k) and keep a slow third mode for s ≤ 0.25. So
  mixing does not cancel the two defects.
- **B, the per-cube mean of the six tetrahedra's K,** is the only candidate with neither problem. Post hoc, it shows
  no genuine growing root anywhere and no third low-frequency mode. Its registered G1 flags were all numerical: Im ω
  at 1e-8 to 1e-6 of |ω|, or numerically split zero roots. But it fails for other reasons:
  - **G2:** one polarisation is strongly anisotropic. The speed error grows as k², e.g. 1.015 → 1.058 → 1.22 at
    kh = 0.025, 0.05, 0.1 along one random direction; 34E had 1e-4 at kh = 0.05.
  - **G3:** a static zero mode with a Jordan chain appears in every direction. It looks like an exact zero-frequency
    re-slicing (likely gauge-like linear drift), but it is unproven.
  - **Time domain (post hoc):** the step-dependent wave drift, 3.4% → 0.8% when the step is halved, and seed growth of
    10⁵–10⁶ point to an integrator problem, not physics. The look is inconclusive.

**Structural finding.** A, B and C each leave one combination of the 7 edge rates per vertex with no kinetic energy:
the edge kinetic matrix has one zero singular value, at 5e-17 relative. That is because 7 edge rates per vertex feed
only 6 K components. That edge combination then becomes an extra constraint, which raises the index of the time
evolution and breaks the midpoint integrator. The per-tetrahedron construction is the only one tested that gives all
7 edges kinetic energy, and it is the unstable one. So the problem is not a matter of choosing a better fit. A stable
version probably needs either the 7th edge removed from the dynamics in a controlled way, or a kinetic term built to
keep the lattice re-slicing symmetry exact.

## Registered results

| candidate | G1 no growth | G2 two modes at c (kh ≤ 0.05) | G3 chains = exact symmetries | G5 benchmarks | G7 real-space | G4 time domain |
|---|---|---|---|---|---|---|
| T (34E) | FAIL (207 roots; max 40 k) | FAIL (125 of 198) | FAIL (41) | PASS (E4 4.2e-4) | PASS (3e-15) | not run (scored only if G1–G3 and G5 pass) |
| A per-cube fit | FAIL (137; max 2.3 k) | FAIL (118) | FAIL (106) | PASS (4.9e-4) | PASS | not run |
| B per-cube mean | FAIL (7 flags, all numerical post hoc) | FAIL (34) | FAIL (226) | PASS (6.3e-4) | PASS | not run |
| C per-vertex fit | FAIL (7) | FAIL (125, third mode 2–2.9 c) | FAIL (44) | PASS (8.9e-4) | PASS | not run |
| D(0.1) | FAIL (176) | FAIL | FAIL | PASS | PASS | not run |
| D(0.25) | FAIL (187) | FAIL | FAIL | PASS | PASS | not run |
| D(0.5) | FAIL (209; max 12.6 k) | FAIL | FAIL | PASS | PASS | not run |

Notes on the table:
- **G5 for every candidate:** the static kernel equals 34D exactly, and the displacement-gauge modes are null to
  1.5e-16.
- **G7:** the real-space direct sums match the Fourier forms to ≤ 3e-15 for all seven, including the ones that failed
  (more than the plan required).

**My stated expectations:**
- T fails G1: held.
- A fails G2: held, and A fails G1 too.
- B behaves like A: wrong. B is the only candidate without growth or a third mode.
- A blend might work: wrong. None did.

**G1's flags are partly numerical.** The registered growth threshold (Im ω > 1e-8 |ω|) catches round-off on
near-degenerate real pairs and numerically split zero roots. The post hoc robust count uses Im ω > 1e-4 |ω| and
|ω| > 0.1 k:

| candidate | robust growing roots |
|---|---:|
| T | 207 |
| A | 108 |
| B | 0 |
| C | 1 (kh = 3) |
| D(0.1) | 159 |
| D(0.25) | 169 |
| D(0.5) | 207 |

Every verdict except B's is unchanged.

## What this means

- The 34D–34E benchmark (exact static field, gravitomagnetism, two light-speed waves) cannot be kept with any of
  these simple strain-rate constructions. The per-tetrahedron one grows. The fitted ones either grow, add a third mode,
  or leave a non-dynamical edge combination.
- The most promising lead is B, the cube-mean. It needs three things:
  - (i) a time integrator suited to the extra constraint, such as index reduction for the non-dynamical edge
    combination;
  - (ii) a proof or test that its static chain is pure gauge;
  - (iii) a check of whether its k² anisotropy is acceptable.

  That would be a new registered study. 34G does not support it as a fix.
- The broader option is to construct the kinetic term so that re-slicing is exact on the lattice. That would be a
  derivation, not a fit. Regge-type literature (Bahr and Dittrich) suggests this is hard at the non-linear level, but
  it may be achievable at linear order.

## Compute

About 1.5 CPU-minutes, well under the 2 core-hour cap:
- spectra: 51 s;
- real-space check: about 10 s;
- time-domain runs: 8 s.

Platform: Linux cloud, Python 3.13.18, numpy 2.5.3, scipy 1.18.1.

## Deviations

- None in the registered runs.
- The post hoc additions are labelled in the files:
  - the robust growth count;
  - the G2 reading of B;
  - B's time-domain look, plus a second look with zero-cluster seeds excluded (identical result: no such seeds were
    present).
