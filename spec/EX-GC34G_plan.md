# EX-GC34G — a lattice kinetic term for the shift-link geometry without the growing mode (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "Go ahead with 34G". This plan and the `explore/LOG.md` entry are
pushed before any 34G computation. Cap: 2 core-hours on the cloud. I expect it to use under 0.5. Our construction,
not LaViolette's Model G.

## Question

34F found that the per-tetrahedron DeWitt kinetic term of 34E has a growing re-slicing mode at every spacing. A
per-cube strain fit removed the growth but added a third propagating mode (0.73–0.85 c). Is there a lattice kinetic
term, built from the same edge and shift variables with the clock rule and the DeWitt form unchanged (λ = 1), that
meets all the criteria fixed in `spec/EX-GC34F_plan.md` for a fix?

## Candidates (fixed now; no others are scored)

All use K_ij fitted from edge strain rates s_e = l⁰(dδl_e/dt − D̂·Δβ) as in 34E, with kinetic energy
(1/2κ) V (K:K − (trK)²) and the unchanged 34D lapse × Regge potential.

| id | construction |
|---|---|
| T | per tetrahedron, exact 6-edge solve, V = h³/6 (34E; control, expected to fail G1) |
| A | per cube, least-squares K from the cube's 19 edges, V = h³ (34F variant a) |
| B | per cube, K = mean of the 6 tetrahedra's K, V = h³ |
| C | per vertex, least-squares K from the 14 edges meeting the vertex (7 outgoing, 7 incoming), V = h³ |
| D(s) | blend s·T + (1 − s)·A, for s = 0.1, 0.25, 0.5 |

No projection is applied to any candidate, so criterion G6 below does not arise.

## Criteria (from the 34F plan's list for a fix; thresholds fixed here)

Methods are those of 34F as finally used, with the zero-shift-gauge roots confirmed by the argument principle
(`src/exgc34f.py`, `analyse_posthoc`). A root counts as found only if confirmed.

- **G1: no growing root.**
  - Directions: the 34F direction set, i.e. axes, all signed face and body diagonals, generic, near-axis, the
    30 random directions, and the three approach sequences.
  - Spacings: kh = 0.0125, 0.025, 0.05, 0.1, 0.2. Along 100, 110 and 111 also kh = 0.5, 1, 2, 3 and the zone
    boundary.
  - "Growing" means Im ω > 1e-8 |ω|.
- **G2: exactly two low-frequency modes.** At kh ≤ 0.05:
  - in every direction, exactly two confirmed real roots with 0 < ωh < 0.5;
  - each with |ω/k − c_light| ≤ 1e-3;
  - every other non-zero root with ωh ≥ 0.5, classified and reported with its kh scaling.
- **G3: no unexplained defective zero roots.** The number of static Jordan chains equals the number of exact
  re-slicing directions in every case.
- **G4: constraints preserved in time.** Generic and (1,1,0) directions at kh = 0.05:
  - Run a wave eigenvector, and the same wave plus every other non-zero eigenvector seeded at 1e-8. Use the implicit
    midpoint rule on the zero-shift-gauge descriptor system, at steps dt = 0.02/ω_max and dt/2, for 40 wave periods.
  - Lapse and momentum residuals ≤ 1e-10 throughout.
  - The wave amplitude changes by ≤ 1e-8.
  - The seeded components grow by no more than a factor 1 + 1e-6.
- **G5: benchmarks kept.**
  - The static kernel restricted to (n, edges) equals 34D to ≤ 1e-12.
  - The displacement-gauge modes are null to ≤ 1e-10.
  - The 34E E4 O(v) Riemann error is ≤ 1e-2 at kh = 0.05 in all 9 cases.
- **G7: separate implementation path, for every candidate that passes G1–G5.** The kinetic energy is summed directly
  in real space on a periodic 8³ lattice for random single-frequency fields. It must match the Fourier form to
  ≤ 1e-10, as in 34E E1.

**A candidate passes if it meets G1–G5 and G7.** Expectations, stated before running:
- T fails G1 (34F).
- A fails G2 (34F: third mode).
- For B, C and the blends D(s) I have no firm prediction. My guess is that B behaves like A. A blend might give the
  extra variable enough kinetic energy to push the third mode to the cutoff while keeping the growth away, but I'm
  not confident.
- If no candidate passes, that is a valid outcome. It would point to building the kinetic term from a principle
  (exact re-slicing on the lattice) rather than by fitting.

## Reported, not scored

For each candidate: the re-slicing breaking (34E E2 measure), the constraint counts, the third-mode speed and its kh
scaling, and the 34E wave-speed table (100, 110, 111 at kh = 0.05, 0.1, 0.2).

## Scope

Linear, flat background, continuous time. Choosing a candidate here does not make it physically preferred. Any
passing candidate would still need non-linear and coupled-matter tests later.

## Files

`src/exgc34g.py` → `results/EX-GC34G/candidates.json` (G1–G3, G5), `evolve.json` (G4), `check.json` (G7). It reuses
`src/exgc34e.py` and `src/exgc34f.py` unchanged.
