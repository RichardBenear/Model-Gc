# EX-GC35E/F: Claude review and next steps

Claude | 2026-11-24 | Requested by the PI (Richard): "review GC35E and F and give me your next steps".

**Scope.** This is a direction review. I read both notes and did no recomputation, apart from one ratio table taken directly from the committed `results/EX-GC35E/tidal.json`. Both studies are pending PI closure. Any suggested run needs the PI's go-ahead and an `explore/LOG.md` entry first. Follows [the 35C/D review](EX-GC35CD_claude_review.md).

## Summary

This is good news on both fronts, though less than the registered verdicts alone would suggest for 35E, and more modest than it may look for 35F.
- **35E:** projection removes the constraint blow-up. The remaining tidal error converges at second order and **stops growing** after t ≈ 12. The registered tidal gate failed on its absolute threshold at h = 0.25, but the behaviour is what a usable benchmark should show.
- **35F:** linear Y1 plus the scalar runs as a local event graph with no global barrier, and every legal schedule gives bitwise-identical results. As the note itself says, this holds by construction: immutable versioned inputs and a single writer per value guarantee order independence. So 35F shows that Y1 *can be expressed* as a causal graph. It is not yet a nontrivial causal-invariance test.

## EX-GC35E

Coarse/fine ratios of the projected tidal error, computed from `tidal.json`:

| t | coarse (h = .5) | fine (h = .25) | ratio | TT-only ratio |
|---:|---:|---:|---:|---:|
| 4 | 9.28% | 2.31% | 4.02 | 3.89 |
| 8 | 19.3% | 4.82% | 4.00 | 3.87 |
| 12 | 32.1% | 8.11% | 3.96 | 3.87 |
| 16 | 36.2% | 9.12% | 3.97 | 4.26 |
| 20 | 35.3% | 8.94% | 3.95 | 3.78 |

- The ratio is about 4 at every time, so the error is a clean O(h²) discretisation error. The error saturates (9.1% → 8.9% from t = 16 to 20) instead of growing like the free run (398%).
- The TT part, which the projection does not touch, converges the same way. My inference: the remaining error is ordinary lattice discretisation of the packet and its source, most likely accumulated dispersion phase error, not an artefact of projection. Extrapolating, h = 0.125 would give about 2.3%.
- **Caveat:** the velocity corrections are large (an impulse scale of about 23–40% of the tidal signal after t = 8). The projection is doing real work, because the next-order source is not conserved at finite h. That is acceptable for a benchmark, but projected runs should always report the correction size, and any claim should rest on gauge-invariant (TT or curvature) readouts.

My view: adopt **projected Y1** as the working reciprocal benchmark, with a stated O(h²) error. Its accuracy is controlled; it just isn't exact.

## EX-GC35F

- The scheduler overlaps depths and needs no barrier. The stale-read control breaks order independence. The common-slice constraint defect converges at second order in dt (×4.01). All of this is correct and useful.
- What it does not test:
  - The events still carry **integer depth labels**, so every site ticks at the same rate. Gravity's main effect, the lapse making clocks run at different rates, plays no part in the event structure.
  - Events never **compete or consume inputs**. Order independence is therefore automatic. The hard version of causal invariance (overlapping matches of a destructive rewrite, as in the parent ModelG engine and GC31B) is untouched.
  - Coupling is one-way only, and the box is a 5³ toy.

## Suggested next steps

### 1. Close out the Y1 benchmark (~1–2 core-h; or Richardson only, near zero)

- Either run 35E once at h = 0.125 to confirm the predicted ~2.3% and the saturation, or Richardson-extrapolate the existing coarse/fine pair and state the result.
- Then fix the reciprocal-run recipe in `docs/`: projection on, correction norms reported, TT/curvature readouts.
- Fold 35A–F into the white paper draft.

### 2. Skyrme action built natively on Y1 (from the 35C/D review, still open)

This is the prerequisite for any moving-knot run. Build the sigma-model and stabiliser terms on Y1's staggered placements so the static Ward identity holds exactly, in the style of 34O. Then the moving knot (a falling knot, then two attracting knots) can run at h = 0.2 with the projected evolution.

### 3. New direction: let gravity set the event structure (paper plus minutes of compute)

Make the lapse set **how many events each site has**, a local clock rate, instead of every site sharing integer depth. In a 1D or 2D toy on a prescribed static lapse field:
- each site fires at its own rate, with interpolation or a version rule at boundaries where rates differ;
- check that schedule independence still holds;
- check that a wave packet's frequency shift (redshift) and its deflection match the continuum lapse prediction.

If this works, the gravitational time dilation becomes a property of the causal graph itself, the event density, rather than a coefficient written into a global stencil. That is the link between the 34–35 gravity work and the programme's causal-invariance goal.

After that comes the harder test: a rule in which events consume their inputs and matches can compete, the true Wolfram-style causal-invariance question. It links back to the GC31B series in this repo and to ModelG's MG189 causal completion work.

### Order

Step 1 is cheap and should come first, so later runs have a fixed recipe. Steps 2 and 3 are independent and can run in parallel.
