# EX-GC09 (EXPLORATORY): route 3, stage 2 – the Gc particle as tokens on a 3D lattice graph

Claude, 2026-11-19, at the PI's request ("go ahead", Stage 2 of docs/ROUTE3_RIBBON_DESIGN.md). Exploratory track.
Gc1 is our construction, not LaViolette's.

## Answer
**Yes: the Gc particle survives as whole tokens on a coarse 3D lattice, and it still draws a 1/r G field. But it
needs a lot of tokens.**
- On a cubic lattice of spacing 0.5 the deterministic particle exists and is steady. Its consumption M is within
  about 10% of the radial value (5.76 at h = 0.5, 5.20 at h = 0.4, 5.32 radial).
- As tokens, it lives for the whole run (T = 30) at Ω = 1000 and Ω = 300, lingers in a shrunken state at Ω = 100,
  and dies within about 10 time units at Ω = 30.
- Its far G deficit is that of a point sink of strength M on the same lattice: within 0.4–1.5% for the
  deterministic lattice, and 4% (Ω = 1000) to 10% (Ω = 300) for the time-averaged token field.

## Results
| Run (h = 0.5) | Box | M (time-average after t = 10) | Far deficit / (M × point-sink Green's function) | Core W tokens |
|---|---|---|---|---|
| deterministic | ±10 | 5.76 (t = 59) | 1.015 (3.5 < r < 9) | — |
| deterministic | ±6 | 5.88 (t = 29) | 1.004 (3.5 < r < 5) | — |
| tokens Ω = 10⁵ (validation) | ±6 | 5.61 at t = 5 (deterministic: 5.61) | — | 3.9 M |
| tokens Ω = 1000 | ±6 | 5.51 ± 0.17 (−6%) | 1.04 | ≈ 34 000 |
| tokens Ω = 300 | ±6 | 4.72 ± 0.16 (−20%) | 1.10 | ≈ 9 000 |
| tokens Ω = 100 | ±6 | 3.23 ± 0.29, still falling (−45%) | 1.50 (not steady) | ≈ 1 700 |
| tokens Ω = 30 | ±6 | dies at t ≈ 10 | — | 0 |

"Core W tokens" is the total number of W tokens at the end (Ω × total W). Fluctuations weaken the particle well
before they kill it: at Ω = 300 it already consumes 20% less than the deterministic one.

## What it means for route 3
- The chemistry runs on the same cubic graph that carries EX-GC08's ribbon rule, so Stages 1 and 2 now share a
  substrate. Nothing about 1/r gravity is lost by going discrete: G is conserved token by token, and the far field is
  the lattice point-sink field.
- A particle needs roughly 10⁴ core tokens to be robust here. That matters for Stage 3: the "tether" will be a large
  bundle of supply paths, not a few threads. EX-GC08 showed the belt-trick topology works for a dense shell, so this
  is not a problem in principle, but it sets the scale.
- Stage 3 still has to define which token paths form the tether and test whether a hand-rotated particle keeps the
  one-turn twist.

## Method
- `src/exgc09_lattice.py`: cubic grid of spacing h, 7-point Laplacian; X, Y, W reflecting at the faces; G held at
  G∞ = 1 in a reservoir layer just outside the faces. Initial state interpolated from the radial reference
  (results/EX-GC01/radial_ref.npz). Deterministic mode: explicit Euler. Token mode: tau-leap with n = Ω h³ c tokens
  per cell; binomial hops (6 directions) and first-order losses, Poisson production and higher-order steps, removals
  clipped to the counts available. dt = 0.0015 throughout.
- `src/exgc09_summary.py` (separate analysis path): reads the records only, solves the lattice Green's function of a
  point sink in the same box with scipy, and writes `summary.json` and `paper/figures/EX-GC09_token_lattice.png`.
- Records: `det_h*_L*.json`, `tok_h0.5_O*_s0.json` (time series of M, total W, peak W, total G, clip count; radial
  G profile, time-averaged from t = 10 for tokens).
- Platform: Linux cloud container (1 core), Python 3.13.18, numpy 2.5.3, scipy 1.18.1. Compute: about 1.5
  core-hours, including the discarded runs below; within the 2 core-hour cap.

## Deviations and limits
1. **A bias in the first token code was found and fixed; those runs were discarded.** The first version applied
   hops to the counts after the reactions. Because reaction and diffusion terms are each about 60 and nearly cancel in
   the core, that ordering shifts the steady state by about 20% (even at Ω = 10⁵ the token M fell to 4.9 while the
   deterministic M rose to 5.6). The fixed version takes hops and reactions from the same start-of-step state. A
   one-step check at Ω = 10⁹ agrees with the deterministic update to noise level, and Ω = 10⁵ tracks the
   deterministic M (5.61 vs 5.61 at t = 5).
2. **The first deterministic runs blew up** (dt too large for the reaction terms); rerun at dt = 0.0015 (h = 0.5) and
   0.0009 (h = 0.4). h = 0.25 was not run (too costly); the radial solution stands in for the fine-lattice limit.
3. **Shorter and smaller than logged:** T = 30 (not about 100) and a box of ±6 for the token runs, one seed each, to
   fit the cap on a single core. Survival at Ω = 300 and 1000 is over T = 30 only. The box of ±6 shifts M up by about
   2% relative to ±10, and limits the far-field test to 3.5 < r < 5 for tokens.
4. **Ω = 100 is contaminated by the tau-leap:** many removals had to be clipped (about 54 000 events), so its numbers
   mix fluctuation effects with step-size error. Ω = 300 and 1000 had 23 and 32 clipped events in total.
5. The token far-field excess (4–10%) may partly be lag: M drifted down during the averaging window, and G responds
   with a delay of a few time units at these distances.
6. **h = 0.5 vs 0.4 differ by 10% in M.** The lattice particle is not converged in h. The survival thresholds in Ω
   are for h = 0.5 only.
