# EX-GC34F (PROPOSED, NOT APPROVED) — why the 34E lattice geometry has a growing mode (EXPLORATORY)

Claude, 2026-11-23. **Draft revision 2**, after ChatGPT's review of revision 1
(`review/EX-GC34F_draft_chatgpt_review.md`, of 46a539b). Nothing here has been run. It needs the PI's go-ahead in words
and an `explore/LOG.md` entry before any run. Our construction, not LaViolette's Model G. Proposed cap: 2 core-hours on
the cloud. I expect it to use under 0.5: Fourier-space linear algebra plus small single-mode time integrations.

**Scope change from revision 1 (accepting the review).** 34F is the *diagnostic* study only: F1–F3 below. A candidate
fix becomes a separate study, 34G. It is registered only once F3 names a specific mechanism. The criteria a fix must
meet are recorded at the end, so that 34G cannot weaken them.

## Question

The 34E vacuum equations have a growing mode in the re-slicing direction (lapse plus lengthwise shift):
Im ω ≈ 1.4–2.1 k off the lattice axes, at every spacing tested (`results/EX-GC34E/review_unstable.json`). Three
questions:
- Does the mode survive a constraint-respecting time evolution?
- How do its growth rate and its curvature content scale as the spacing shrinks at fixed physical wavenumber?
- Which terms of the unmodified action admit it, while the static and moving-source equations still hold?

## Setup common to all parts

- **Equations.** The committed 34E kernel `full_kernel` (unchanged), in Fourier space for one wavevector q. With the
  convention e^{−iωt}, M(ω) = C + ωB + ω²A, so the time-domain equations are C x + iB ẋ − A ẍ = 0. A, B and C are
  taken from M at ω = 0, ±1; M is exactly quadratic in ω (checked to 4e-16).
- **Fixed physical wavenumber.** With k fixed, a smaller spacing h means a smaller kh. All scans in kh are therefore
  read as h → 0 at fixed k, and results are reported both as ω/k and as ωh.

## F1 — complete spectrum, without assumed symmetries

1. **Symmetries first.** Before any reduced direction grid, test numerically which of the 48 cube symmetries (signed
   permutations) map the assembled kernel's spectrum onto itself. Compare the full sorted root sets at 20 random q, to
   ≤ 1e-10. The Kuhn triangulation singles out the (1,1,1) diagonal, so I expect only the 6 permutations, or 12 with
   inversion. **Directions are then sampled over the fundamental domain of the symmetry group actually found, not
   over 1/48.**
2. **Directions.** All signed pairs that the triangulation can tell apart: (1,±1,0), (1,0,±1), (0,1,±1), (1,±1,±1).
   Also a generic direction, and 30 random directions in the fundamental domain. Three sequences approach an axis or a
   face diagonal from off it, at angles 10⁻¹, 10⁻², 10⁻³ and 10⁻⁴ rad, where the constraint rank can change.
3. **Spacings.** kh = 0.0125, 0.025, 0.05, 0.1, 0.2 (long wavelength), and kh = 0.5, 1, 2, 3 and the zone boundary
   along the three main directions. ChatGPT's earlier check found an extra oscillatory root along an axis at kh = 3;
   it is to be reproduced and classified.
4. **Method.** Explicit deflation rather than a pencil:
   - Remove the displacement-gauge null space (`gauge_modes`, 3 vectors).
   - Remove the exact re-slicing null directions found at that q, counted by the singular values of M at four
     generic real ω below 1e-15 relative.
   - Solve the reduced quadratic eigenproblem, keeping all complex roots.
   - Check every eigenvector against all 11 equations: relative residual ≤ 1e-10.
   - Cross-check every root as a minimum of the deflated kernel's smallest singular value, at ≤ 1e-12 relative,
     against ≥ 1e-6 at neighbouring ω.
5. **Classification, reported separately.**
   - Zero modes: exact gauge, with multiplicity and any Jordan (defective) structure, tested by the rank of
     [M(0), M′(0)].
   - Propagating modes: real ω. Report the polarisation (transverse-traceless content of h_ij) and the speed against
     the 34A light speed.
   - Growing modes: Im ω > 0. Report multiplicity, field content, non-gauge fraction, and scaling with kh.
   - Any other finite roots.

**Registered expectations (F1):**
- (a) Off-axis growing roots reproduce `review_unstable.json` to 1e-6 in ω/k.
- (b) The axis directions are decided: growing, or not, after explicit deflation.
- (c) Every root found passes both checks in item 4.
- (d) At kh ≤ 0.05, exactly two propagating modes per direction, with |c − c_light| ≤ 1e-3.
- (e) Reported, not predicted: the kh = 3 axis root, and the multiplicity of the growing branch.

## F2 — time-domain check with constraints enforced

1. **Formulation.** Write the gauge-fixed equations as a first-order descriptor system E ż = F z, with
   z = (x, ẋ). The lapse and shift carry no second time derivative, so E is singular and this is a differential-algebraic
   system. Its algebraic rows are the lapse (Hamiltonian) and shift (momentum) equations. I first state its index and
   its constraint manifold from the null space of E.
2. **Initial data.** Initial data must satisfy the algebraic rows *and* their first time derivative (consistent
   initialisation), projected onto that manifold. Starting data:
   - (i) the growing eigenvector;
   - (ii) a gravitational-wave eigenvector;
   - (iii) the wave eigenvector plus the growing eigenvector seeded at a known amplitude 1e-8 (not round-off).
3. **Integrator.** An implicit, A-stable, descriptor-form scheme (implicit midpoint), with time steps dt = 0.02/|ω|
   and dt/2. Generic direction at kh = 0.05, and the face diagonal. Run for 15/Im ω, about 15 e-foldings.
4. **Readouts.**
   - The lapse and momentum constraint residuals, each separately, at every step.
   - The growth of the reduced mode amplitude.
   - The growth of a curvature (tidal) observable: the linearised electric part E_ij = R_0i0j from `metric_from_x`
     and `riemann`. The mode is nearly a coordinate change, so the two growth readouts must be compared.

**Registered expectations (F2):**
- (a) Constraint residuals stay ≤ 1e-10 relative throughout, at both time steps.
- (b) Run (i) grows at Im ω, within 1% in the mode amplitude and in E_ij, at both time steps. The two time steps
  agree to O(dt²).
- (c) Run (ii) oscillates at the F1 frequency, with |amplitude change| ≤ 1e-8 over the run.
- (d) In run (iii), the seeded component grows at Im ω from 1e-8. The E_ij it carries, as a fraction of its mode
  amplitude, is ≈ 0.06 (kh)², within a factor 2.

## F3 — mechanism in the unmodified action

1. **Small-kh reduction.** Restrict the unmodified kernel to the growing mode's subspace: lapse, lengthwise shift and
   the edges it couples to. Expand the 2×2 or 3×3 reduced kernel in kh to the leading order that fixes ω² for the
   branch. Identify which matrix elements set ω²:
   - per-tetrahedron kinetic (DeWitt) elements: shift–shift, shift–edge;
   - Regge lapse–edge elements;
   - the cross terms between them.
   Done symbolically where feasible (sympy on the tetrahedron strain map), and checked numerically against F1.
2. **Constraint bookkeeping.** At each q, record which constraints survive: the rank of the lapse and shift rows on
   the constraint surface, and whether the re-slicing generator is a null vector. Explain why the branch has
   ω² ∝ −k² rather than a small, (kh)²-suppressed value, i.e. which O((kh)²)-suppressed kinetic and potential terms
   balance.
3. **Single changes, read as diagnostics only, not as fixes.** Each is re-reduced with item 1. Report which
   constraints survive, and whether the static (34D) and moving-source (34E E4) equations still hold.
   - Per-cube rather than per-tetrahedron fit of the strain rate.
   - The DeWitt trace coefficient λ ≠ 1. This changes the continuum constraint structure, so it is a diagnostic of
     which term matters, not a candidate.
   - The lapse–Regge coupling. This removes the Newtonian static response, so the same caveat applies.

**Registered expectation (F3):** a named term or term pair whose leading-order elements reproduce the F1 growth rate
of the unmodified action to 1%, in at least the 111 and generic directions. If no such reduction is found, F3 is
reported as unresolved.

## Main readout

Growth rate and curvature content of the unstable branch at fixed physical k as h decreases. This decides whether the
34D–34E tensor benchmark can be kept with a repaired kinetic action, and what a repair has to preserve.

## Criteria a later fix must meet (for 34G; recorded now, not tested in 34F)

- No growing root at any direction or spacing sampled as in F1.
- Exactly two healthy low-frequency tensor polarisations; every other mode classified, with its scaling.
- No secular growth from defective (Jordan) zero roots.
- The lapse and momentum constraints are preserved in time evolution, as in F2.
- The 34D static kernel is unchanged (≤ 1e-12). The 34E E4 gravitomagnetism agrees with GR to ≤ 1e-2 at kh = 0.05.
- If the fix projects out the re-slicing direction, the projection comes from a symmetry or a constraint the
  evolution preserves, not from an imposed cut.

## Files (planned)

- `src/exgc34f.py` → `results/EX-GC34F/spectrum.json` (F1), `evolve.json` (F2), `mechanism.json` (F3).
- A separate recomputation path: F1 roots by singular-value minimisation, independent of the eigen-solver; F2
  growth checked against F1.
- It reuses `src/exgc34e.py` unchanged.

## Response to the review (revision 1 → 2)

| Point | Response |
|---|---|
| Register F1–F3 as diagnostics; F4 separately | Accepted. F4 moves to a future 34G; its criteria are recorded above. |
| Prove symmetries before reducing the grid; signed pairs; approach to axes; zero modes separate; larger kh incl. kh = 3 | Accepted: F1 items 1–3 and 5. |
| Constraint-preserving or DAE solve; consistent initialisation; two time steps; separate lapse/momentum residuals; seeded control; curvature readout | Accepted: F2. |
| Mechanism from the unmodified action first; trace coefficient and lapse–Regge changes alter the constraints | Accepted: F3 items 1–2 come first; item 3's changes are labelled diagnostics, not candidates. |
| F4 criteria: two tensor polarisations, classification, defective zero roots, constraint preservation, static and gravitomagnetic checks, symmetry-derived projections | Accepted, in the 34G criteria above. |
| Main readout at fixed physical k as h decreases | Accepted. Because ω/k depends only on kh and direction, the kh scan is that readout; stated explicitly. |
