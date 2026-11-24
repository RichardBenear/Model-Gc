# EX-GC35A/B: Claude review of direction and suggested next steps

Claude | 2026-11-24 | Requested by the PI (Richard): "ChatGPT ran your suggestions and the results are in 35A/B. What do you think of those and where do we go next?"

**Scope.** This is a direction review. I read `results/EX-GC35A/EX-GC35A_note.md`, `results/EX-GC35B/EX-GC35B_note.md`, `docs/EX-GC34_symmetry_and_direction.md`, `HANDOFF.md` and `explore/LOG.md`. I did **not** recompute any result. Both studies are pending PI closure. The suggestions below need the PI's go-ahead and an `explore/LOG.md` entry before running. Follows [the 34Q–Z review](EX-GC34Q-Z_claude_review.md).

## Summary

Both runs were done carefully and reported honestly. Neither yet answers the question it was meant to answer.
- **35A** passes its gates, but its own numbers show the next-order constraint defect growing fast within t ≤ 4. "Bounded or growing?" is therefore still open, and the data lean towards growing.
- **35B** shows that a knot's energy sources Y1 cleanly, binding included. The full source (energy plus stress) misses the compatibility gate by 1.5–2%. The most likely cause is that the knots are not fully relaxed or resolved, not a defect of Y1. That needs a direct check before anything is built on it.

Step 3 of the previous review (the paper note `docs/EX-GC34_symmetry_and_direction.md` and the white paper v0.5 update) is done, and I agree with its conclusion. It is acceptable to keep causal invariance as the target while using approximate lattice coordinate symmetry, *provided* the constraint errors stay controlled over long times. That proviso is what 35A has not yet shown.

## EX-GC35A: what the numbers say

Next-order constraint residual (from the note's table):

| Run | t=1 | t=2 | t=4 | growth t=2→4 |
|---|---:|---:|---:|---:|
| Coarse 24³ | 7.99e-6 | 2.74e-5 | 5.70e-4 | ×21 |
| Fine 48³ | 2.08e-6 | 7.64e-6 | 1.72e-4 | ×23 |
| Continuum control | 5.9e-10 | 1.4e-9 | 4.8e-9 | ×3.4 |

- The growth **accelerates**: about ×3.5 from t=1 to 2, then about ×22 from t=2 to 4. That is faster than any low power of t.
- The growth factor is **the same on both grids**. Refinement lowers the amplitude (×3.3, observed order 1.73) but not the rate. If that continues, the defect is not "small and h²-controlled". It is an O(1)-rate mode with an h-dependent seed, the same kind of problem as the 34E growing mode.
- The continuum control grows far more slowly (×3.4), so most of the lattice growth is a lattice effect.
- Fine-grid curvature error rises from 1.35% to 5.27% between t=2 and t=4.
- Caveat: a mean next-order energy of 5.4e-3 is left out of the metric solve (the periodic homogeneous mode), so part of the growth may come from the periodic compensator rather than from Y1 itself. I have not checked this; it is an inference.

The note itself says t ≤ 4 cannot decide long-time behaviour. I agree, and I think the short-window trend is a warning rather than reassurance.

## EX-GC35B: what the numbers say

- **Registered run failed for a resolution reason.** The B1 knot unwinds at h = 0.3. GC26 worked at h ≤ 0.2, so h = 0.3 was too coarse (a plan slip, not physics).
- **Post hoc h = 0.2:**
  - Energy-only sourcing is clean: tidal error 0.05% vs continuum.
  - Binding is included: B2 is 3.9% below 2×B1, and the monopole ratio equals the energy ratio. The note correctly says this is a normalisation check, not an independent active-mass measurement.
- **The full source misses by 1.55–2.09% (gate 0.1%).** A true static solution has divergence-free stress, which is exactly what the static Y1 equations need. A residual of this size points to a knot that is not quite a static solution:
  - both optimisations hit the 350-iteration cap;
  - the degree integrals are −0.93 and −1.88 instead of −1 and −2 (under-resolved cores);
  - the box has fixed-vacuum walls.

  The note says the run cannot separate these from a real incompatibility. Until they are separated, the full-stress failure should not count against Y1.
- Feedback was checked only as the initial acceleration at rest, not as an evolution.

## Suggested next steps

### Step 1: long-time 35A growth test (~0.5–1.5 core-h; ChatGPT, cloud is fine)

- Extend 35A to t ≈ 20–40 on 24³ and 48³ (96³ only if affordable).
- Fit the growth law of the next-order residual (exponential rate or power) on each grid.
- Decisive question: does the **rate** fall with h (acceptable: an approximate-symmetry error that vanishes in the continuum), or only the **amplitude** (unacceptable: a lattice instability)?
- Control: repeat with the homogeneous mode handled differently, or with the compensator off, to see whether the periodic compensation drives the growth.
- If a mode grows at a fixed rate, find its wavenumber and polarisation. Then compare it with the 34E re-slicing mode and with Y1's linear spectrum.

### Step 2: 35B stress-residual convergence check (~1–2 core-h; PC for the finest grid)

- Relax the B1 and B2 knots to tighter convergence (no iteration cap; gradient tolerance stated in advance).
- Measure the full-source residual against h = 0.25, 0.2, 0.15 and against box size.
- Report the integrated stress (virial/Derrick balance) alongside the residual.

If the residual falls like h² and with convergence, the full source is compatible and 35B's failure was preparation. If it stays at about 1–2%, the Skyrme stress mapping onto Y1 is wrong somewhere: check the metric variation of the Skyrme energy against the continuum stress tensor term by term.

### Step 3 (only after 1 and 2 pass): a moving knot on Y1

Repeat 26E/26F in Y1 geometry: a knot falling in a prescribed Y1 field, then two knots attracting through their own Y1 field. This is the first test that puts Model Gc particles and GR-like gravity together dynamically (~2–4 core-h, PC).

If step 1 shows a fixed-rate growing mode, I would pause steps 2 and 3. Fixing the Y1 evolution, for example with constraint damping or a constrained or projected evolution, would then be the priority, since nothing long-lived can be built on a growing defect.

Recommended order: steps 1 and 2 can run in parallel; step 3 only after both pass.
