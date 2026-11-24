# EX-GC34Q–34Z: Claude review of direction and suggested next steps

Claude | 2026-11-24 | Requested by the PI (Richard): "Take a look at the 34Q through 34Z results that ChatGPT produced. No need to check all the calculations. What I would like is your opinion on what might be the next steps."

**Scope.** This is a direction review. I read the notes in `results/EX-GC34Q/` through `results/EX-GC34Z/`, `explore/LOG.md` and `HANDOFF.md`. I did **not** recompute any result. All ten studies are already PI CLOSED. The next steps below are suggestions only and need the PI's go-ahead and an `explore/LOG.md` entry before running. The PI intends ChatGPT to run all three.

## Summary

The linear Y1 gravity, with a scalar wave packet as source, now holds up well. The series from 34S and 34V to 34Z gave diminishing returns. I recommend we stop searching for an exact finite-spacing nonlinear (cubic) gauge symmetry. The next question to test is whether the approximate symmetry, which is broken only at O(h²), is good enough in practice.

## What was established

- **34Q:** added scalar kinetic energy changes the tidal curvature as linearised GR predicts: 0.031–0.771% error at h = 0.25, with refinement ratios near 4, so second order.
- **34R:** the scalar packet's leading response to its own zeroth-order Y1 metric matches the continuum (0.107% and 1.267%, ratios 4.04 and 3.64).
- **34T:** continuum next-order conservation (quadratic Bianchi) identity verified symbolically, plus a homogeneous Friedmann benchmark.
- **34U:** an isolated spherical Einstein–Klein–Gordon reciprocal (two-way) continuum benchmark passes. It is a useful reference for any lattice reciprocal test.

## What failed

34S, 34V, 34W, 34X, 34Y and 34Z all failed on the same requirement: an exact finite-spacing identity (symmetric, conserved next-order source, or the cubic on-shell gauge identity) on Y1.
- 34S: direct re-sourcing with the linear kernel breaks the constraints (7.6e-5 vs 1.3e-8).
- 34V and 34W: the quadratic and variational cubic candidates fail exact compatibility.
- 34X, 34Y and 34Z: short-range action families, staggered placement, then a joint action-plus-transformation fit. 34Z fits its training set to roundoff but fails on held-out triads. A post hoc lower bound (0.00179 at the training spacing) rules out the whole tested 120-parameter family.

In every case the residuals fall roughly as h². The symmetry is broken by a small amount that disappears in the continuum limit.

## Concerns (direction, not errors)

1. **Fitting without a structural principle.** Five successive fits for the same identity, each enlarging the family, look like a rabbit hole. Lattice gravity results (for example Regge calculus) indicate that exact diffeomorphism symmetry is generally broken at finite spacing for nonlinear configurations. An exact local solution may not exist in any family of this type.
2. **Exactness may not be what the programme needs.** The programme requires causal invariance, a property of the update rule (independence of update order). Exact lattice coordinate symmetry is a different property. Approximate symmetry that becomes exact in the continuum limit may be enough, provided it does not cause the constraints to grow over time.
3. **Drift away from particles.** No particle (the GC26 knot) has yet been put on Y1. The 34 series has been working with scalar wave packets only.

## Suggested next steps (all three; the PI plans to have ChatGPT run them)

### Step 1: two-way scalar-packet and gravity run on Y1 with the approximate symmetry (~0.5–1 core-h)

"Packet" here means the **scalar wave packet** of 34O–34R, not a particle. Evolve the packet and the Y1 metric together, with feedback in both directions. Use the existing 34R/34V next-order source as it stands, without an exact finite-spacing identity. Suggested readouts:
- constraint residuals over time, and whether they stay bounded at about h² size or grow secularly;
- how they scale with h, and with dt so time-step error can be separated out;
- comparison of the mass and flux readouts with the 34U continuum benchmark (or a lattice analogue of it);
- a control: the same run with the coupling switched off, or with the linear kernel only, as in 34S.

This settles whether the 34V–Z exactness gap matters in practice. Bounded O(h²) constraint error would mean we can proceed without an exact identity. Secular growth would mean the symmetry problem is real and needs a structural fix (see step 3).

### Step 2: put the GC26 knot particle on Y1 (~1–3 core-h; probably the PI's PC)

The **particle** here is the degree-1 (and degree-2) Skyrme knot from EX-GC26. First run one-way: does the knot's total energy, binding included, source Y1 lapse, shift and tidal curvature as linearised GR predicts, as in 34P/34Q? Then add the leading feedback, as in 34R. This reconnects the gravity work to Model Gc particles, and repeats the 26E/26F energy and binding results with the Y1 geometry in place of the Nordström field.

### Step 3: paper-only note, no compute

- Does causal invariance (rule-order independence) require exact lattice coordinate symmetry, or is it enough for that symmetry to emerge in the continuum?
- Formulations where the nonlinear gauge symmetry is exact by construction, for example first-order tetrad and connection (lattice-gauge-style) gravity, and what they would cost compared with Y1.
- Fold 34F–34Z into white paper v0.5 (still a DRAFT).

Recommended order: step 1 first, because its outcome decides how much step 3 matters. Step 2 can run in parallel.
