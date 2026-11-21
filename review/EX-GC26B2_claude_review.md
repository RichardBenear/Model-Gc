# EX-GC26B2 — Claude review (EXPLORATORY)

Reviewer: Claude, 2026-11-21. Reviewed commits 897db1c (plan, with execution amendment) and d33092d (results), by
ChatGPT. Plain massless Skyrme model; the numerical construction is ours.

## Verdict

Sound, carefully reported, and more useful than its FAIL lines suggest. I agree with every disposition, including
the failures. Precision 26E should stay blocked on these gates.

The study shows that my 26B explanation was only half right. Restoring mirror symmetry removed the sideways momentum,
by a factor of about two million. But the along-axis speed mismatch did not change: 0.9823 before and 0.9823 after
at h = 0.2, u = 0.3. The good diagonal result in 26B was a property of the old stencil, not of the physics. ChatGPT's
exact energy-current calculation gives the real cause, and I confirmed it independently (below).

## What I checked

1. **Recomputation from the committed time series** with my own fitting code (not `exgc26b2_figure.py`). The figure
   script is ChatGPT's own recompute; mine is separate. All five short boosts, the T = 5 coast (early 0.297340, late
   0.295466, loss 0.630%) and the forced runs (a, F/a) reproduce to every printed digit.
2. **Independent derivation and check of the energy-current mechanism, on the old stencil with my own 26A code.**
   - I derived the bond current myself. One-sided differences give J_a = Σ_i Σ_s w (∂L_i/∂d_aφ)·v at the neighbour
     i + s·e_a. The canonical momentum is the same sum with v at the site i.
   - I evaluated both on my own boosted 26A state (h = 0.2, box ±5, two-orientation stencil, u = 0.3).
   - Results: J/P = 0.9815 (axis) and 0.9988 (diagonal). My fitted centroid speed ÷ (P/H) from the 26B review rerun
     was 0.9817 and 0.9988, and J/H predicts the centroid speed to within 3e-4. So the mechanism also holds in a code
     ChatGPT did not write.
   - Script and output: `review/EX-GC26B2_check/`. This uses the same solver author as 26A, so it is independent of
     ChatGPT's code, not of mine.

## What the results mean

- **The lattice has no exact momentum.** It has two lattice versions of the continuum momentum, which differ at order
  h²:
  - the canonical momentum P, which uses the velocity at each site;
  - the energy current J, which uses the velocity at the neighbouring site.
  The energy centroid moves at exactly J/H, so "mass from momentum ÷ centroid speed" carries the J/P gap. At u = 0.1
  the gap is 1.42% at h = 0.2 and 0.71% at h = 0.15. That ratio of 2.0 is close to the h² ratio of 1.78.
- **Neither number is the knot's "true" lattice mass.** This is truncation error that shrinks with refinement. It is
  not a failure of relativistic inertia. I agree with ChatGPT not to redefine P = J to pass the gate.
- **The force test fails mostly on the momentum input, not on the inertia.** This is a post-hoc decomposition only,
  with no rescoring:
  - With the core-localized (squared-deficit) force, momentum ÷ centroid speed gives 1.008 E₀ after t = 0.6, within
    1%. But the momentum gained is only 0.96 of the nominal impulse fQ. Together these make F/a = 1.056.
  - With the linear force, momentum ÷ speed gives 0.846: the centroid is pulled around by the deformed tail, as the
    note says.
  - So the open question is why dP/dt falls about 4% short of fQ at h = 0.2. That is larger than the 1.4% J/P gap,
    and it may come from Q growing as the knot deforms.
- **The forced runs are very short.** The knot moves 0.037, under a fifth of a lattice spacing, and the force switches
  on suddenly. Both can bias an acceleration fit: one through grid pinning, the other through ringing of the knot's
  internal modes. A smooth ramp and a longer run would separate inertia from start-up effects.
- **Long coast.** Speed loss of 0.63% over T = 5 passes. The momentum drift of 2e-3 E₀ fails. The regional flux shows
  outward transport (R = 3 drift drops from 4.8e-3 to 0.9e-3 after correction), so part of this is radiation leaving
  the core, not drag.

## Recommendation

The numerical error is now understood. It is a ~1.4% (h = 0.2) and ~0.7% (h = 0.15) gap between two lattice
momenta, scaling about as h², plus force-coupling effects. Getting below 0.5% would need a new higher-order
static-and-kinetic operator: new code, new validation, and more compute. Given the programme's goal and the PI's
preference to skip runs that add only a little precision, I'd stop stencil qualification here and move on:

1. Record a lattice error budget for later stages: about 1.5% at h = 0.2 and 0.7% at h = 0.15, scaling as h².
2. Make 26E a **comparative** test rather than an absolute one. Measure the fall of different knots on the same
   lattice, at matched resolution (same core width / h). The common lattice error then largely cancels. An absolute
   0.5% calibration of a continuum theorem is not needed to answer whether fall is universal.
3. Before 26E, one short force step only: core-localized coupling, smooth ramp, T ≈ 5, at h = 0.2 and 0.15. Report
   the momentum ledger dP/dt vs fQ and P/v separately. It is cheap and directly relevant to the measurement 26E will
   make.

ChatGPT's higher-order operator remains a valid option if the PI wants precision for its own sake. I'd rank it below
the steps above.

Compute: 26B2 used 0.86 core-hours of the 2 cap. My checks were seconds plus the reuse of my 26B rerun.
Protocol: clean. The plan reached main 20 minutes before the results. The force-coupling amendment is correctly
labelled, reported alongside the failed linear pair and not substituted for it.
