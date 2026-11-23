# ChatGPT review of the EX-GC34F draft

2026-11-23 | Reviewed draft at `46a539b` | Recommendation: revise before registration.

GC34F asks whether GC34E’s growing mode survives constrained time evolution and what in the lattice action enables it. Register F1–F3 as the diagnostic study. Register F4 separately once there is a specific candidate fix.

## F1. Define the spectrum and direction grid

Prove the symmetries of the assembled kernel before restricting the grid to 1/48 of directions, or sample signed directions explicitly. The fixed Kuhn triangulation privileges particular diagonals. Include pairs such as (1,1,0) and (1,−1,0), plus a sequence approaching an axis where the constraint rank can change. Report zero modes and their multiplicities separately from propagating and growing modes. Add a few larger-kh points or call this a long-wavelength spectrum. An earlier independent check found an additional oscillatory root along an axis at kh = 3, outside the proposed range.

## F2. Specify constraint enforcement and growth measurements

The lapse equation constrains the edge variables. Use a constraint-preserving reduction or a differential-algebraic solve, and initialize both the perturbation and its time derivative consistently. Repeat at a smaller timestep and report the lapse and momentum residuals separately. For the gravitational-wave control, include a run with a deliberately seeded growing component of known amplitude. Round-off is an unreliable seed. Measure growth in both the reduced mode amplitude and a curvature or tidal observable because the mode is close to a continuum coordinate change.

## F3. Identify the mechanism before varying the action

Reduce the unmodified action first and track the leading small-kh terms that set the unstable branch’s frequency squared. Record which constraints survive, then repeat after each proposed change. Changing the DeWitt trace coefficient can alter the continuum constraint structure. Switching off the lapse–Regge coupling can remove the constraint that produces the static Newtonian response. Explain which terms enable growth while the intended static and moving-source equations are still in place.

## F4. Strengthen the candidate-fix criteria

Require two healthy low-frequency tensor polarizations and classify additional modes, including their scaling. Check for secular growth from defective zero roots and verify constraint preservation. Retain the static-kernel and gravitomagnetism checks. If a candidate projects out the re-slicing direction, derive the projection from a symmetry or a constraint preserved by evolution.

## Main readout

Measure how the growth rate and observable curvature scale at fixed physical wavenumber as the lattice spacing decreases. This will clarify whether the tensor benchmark can be salvaged and what a repair has to preserve.

## Records

- [GC34F draft at the reviewed commit](https://github.com/RichardBenear/Model-Gc/blob/46a539b3c43d93d788a846dd9c56dd52093ab150/spec/EX-GC34F_plan_DRAFT.md)
- [GC34E equations](../src/exgc34e.py)
- [GC34E post hoc mode analysis](../results/EX-GC34E/review_unstable.json)

This is a review of the proposed plan. No GC34F runs were performed.
