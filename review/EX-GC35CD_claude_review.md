# EX-GC35C/D: Claude review and next steps, including new directions

Claude | 2026-11-24 | Requested by the PI (Richard): "Both 35C and 35D are on the repo for your review. Any new ideas?"

**Scope.** This is a direction review. I read `results/EX-GC35C/EX-GC35C_note.md` and `results/EX-GC35D/EX-GC35D_note.md` and did **not** recompute anything. A few readings below are my own inferences from the published tables; they are marked as such. Both studies are pending PI closure. Any suggested run needs the PI's go-ahead and an `explore/LOG.md` entry first. Follows [the 35A/B review](EX-GC35AB_claude_review.md).

## Summary

Both registered gates failed, and both notes report it plainly. Read together, though, the two failures look less alarming than "the approximate Y1 model does not work":
- **35C:** the growth is most likely the known drift of an unconstrained synchronous-gauge evolution driven by an O(h²) constraint error. It does not look like a lattice instability, so it should be fixable by controlling the constraints.
- **35D:** the knot's stress mismatch converges with refinement at about second order. The finite-spacing remainder looks like the lattice Derrick (virial) offset already seen in 34B. A matter action built natively on Y1, as 34O did for the scalar, could remove it exactly.

## EX-GC35C

What the data show:
- The power-law fits (exponent about 4.3–5.5) beat the exponential fits, and the fitted exponential rate *falls* from 1.6 to 0.27 across windows. That is not a fixed-rate unstable mode.
- 35C's vacuum analysis rules out a 34E-type exponential mode: Q²(Q+I) = 0, with a defective zero sector.
- The coarse/fine ratio stays about 3.3–4.2 at every time (t = 4: 3.3; t = 8: 3.7; t = 12: 4.1; t = 20: 4.2, from the table). So at any fixed t, the defect still converges at roughly second order. I worried in the 35A/B review that refinement might lower only the amplitude and not the rate. The data say refinement lowers the whole curve by about 4× per halving, which is what an O(h²) seed amplified by gauge drift looks like.
- The compensator is not the cause: switching it off changes the result by 1%.

**My reading (inference).** The vacuum note gives S̈ = 0 and L̈ = −|K|²S/2 for constraint-violating data. A source that violates the constraints by O(h²) feeds this defective chain: S grows like t², and L like t⁴. That matches the observed t^4–5. The error is then about h²·t⁴. That is gauge/formulation drift, the classic weakness of free synchronous-gauge evolution, not a geometric defect of Y1. In numerical relativity the standard cure is constraint damping or constrained evolution, not a better lattice.

## EX-GC35D

What the data show:
- All eight cases relaxed to the tighter tolerance.
- Wall size and tolerance barely matter (about 1%).
- The residual falls with h at an observed order of about 2.3–2.4: 2.96% → 1.65% → 0.87% for B1.

So my 35A/B guess that the knots were unrelaxed was **wrong**. The remaining mismatch is a resolution effect. Two signs:
- The integrated stress trace (E₄ − E₂, which is zero for an exact continuum soliton by Derrick's theorem) is 4–12% of E and falls with h.
- The degree integrals are still −0.96 and −1.93 at h = 0.15.

Extrapolating the trend, the 0.1% gate would need h ≈ 0.06 (about 200³ sites). That is too expensive for routine moving-knot runs.

## Suggested next steps

### 1. Constraint-controlled rerun of 35C (~1 core-h, cloud)

Keep Y1 and the 35C source unchanged. Remove the secular drift in one of two ways:
- **Projection:** after each step, project the second-order metric onto the per-mode constraint surface. Y1 is solved mode by mode and the 8-dimensional constraint kernel is already known, so this is cheap.
- **Damping:** add Z4-style constraint-damping terms with a damping rate κ set in advance.

Gate: over t ≤ 20 the defect stays bounded, converges as h², and the next-order tidal error vs continuum stays at the percent level. If projection works, it is also a ready tool for every later reciprocal and knot run.

### 2. A Skyrme action built natively on Y1 (derivation plus ~0.5–1 core-h)

Instead of pushing h down, apply the 34O trick to the knot. Place the sigma-model fields and their products on Y1's staggered sites, links and plaquettes so that the static Ward identity (stress divergence = Y1 constraint compatibility) holds exactly at finite h, as the product-form stress did for the scalar. The quadratic (sigma-model) term should go through as in 34O, since it is four scalars with a constraint. The quartic Skyrme term is the open part. A fallback is a different stabiliser with simpler products, for example a sextic baryon-current term or a potential. Gate: full-source residual at roundoff at h = 0.2. After that, the moving knot (the old step 3) can run at affordable resolution.

As a cheaper alternative, Richardson-extrapolate the 35D h = 0.25/0.2/0.15 series and accept a stated O(h²) source error. That is good enough for qualitative moving-knot tests, but not for exact gates.

### 3. New direction: test causal invariance of Y1 itself (paper first, then minutes of compute)

So far, all the gravity work (34A–35D) used synchronous global time steps. The programme's stated requirement, causal invariance, has not been tested on any of it. The linear Y1 scalar + gravity system is local and staggered, much like a leapfrog. Ask whether it can be written as an **asynchronous local event rule**: each site or link updates when its neighbours' inputs are ready, with no global clock. Then test whether different legal update orders give the same final field (exact for linear rules if the update is a fixed local map; worth checking). If they do, it connects the 34–35 gravity results to the parent ModelG causal-invariance goal. It also gives a principled place for constraint preservation, because constraints that hold on every local event cannot drift. If they don't, we learn which part of Y1 needs a global time slice, which is important to know before investing further.

This is a design note plus a small numerical check (a 2D or 3D toy box with randomised schedules, compared bit for bit or to roundoff). It returns to the programme's main goal rather than adding another increment.

### Order

Steps 1 and 3 are independent and can run in parallel. Step 2 is the prerequisite for any moving-knot run. If step 1's projection fails to bound the defect, Y1's reciprocal coupling would need rethinking before steps 2 and 3.
