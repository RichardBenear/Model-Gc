# EX-GC34F — why the 34E lattice geometry has a growing mode (EXPLORATORY)

Claude, 2026-11-23. **Registered before any run.** PI go-ahead (Richard, 2026-11-23): "run 34F". This is draft
revision 2, which followed ChatGPT's review (`review/EX-GC34F_draft_chatgpt_review.md`), plus the registration edits
listed at the end. Those edits were made before any 34F computation and only make methods and thresholds concrete. Our
construction, not LaViolette's Model G. Cap: 2 core-hours on the cloud. I expect it to use under 0.5.

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

- `src/exgc34f.py` → `results/EX-GC34F/spectrum.json` (F1), `evolve.json` (F2), `mechanism.json` (F3). Plus a
  note.
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

## Registration edits (2026-11-23, before any 34F run)

These make the revision-2 text concrete. None was chosen after seeing 34F data.

1. **F1 gauge and root method.** Roots are found in the zero-shift gauge (β = 0). That gauge is reachable at any
   ω ≠ 0, because the displacement gauge shifts β by −iωu. The unknowns are then n and the 7 edges. The rows are the
   lapse row and the 7 edge rows. The 3 shift rows (the momentum constraint) are not imposed; they are checked as
   residuals.
   - The normal-rank deficit n_k of this 8 × 8 quadratic is measured at four generic real ω.
   - Exact null directions are removed by a rank-completing random perturbation of the 16 × 16 linearised pencil
     (Hochstenbach, Mehl and Plestenjak 2019). Its true eigenvalues do not move; spurious ones do.
   - Two independent random perturbations are applied, and a root is kept only if it appears in both, to 1e-8
     relative.
   - Zero roots (|ω| ≤ 1e-8 k) are gauge in this gauge. They are reported as zero modes, separately from the static
     zero-mode analysis below.
2. **Validation thresholds (replacing "≤ 1e-12 relative, ≥ 1e-6 at neighbours").** On the 8 × 8 kernel, the kernel's
   relative singular values span many orders of magnitude at small kh. So the singular-value check is a *contrast*:
   - σ_r at the root, where r = 8 − n_k is the smallest generically non-zero singular value, must be ≤ 1e-6 times the
     smallest σ_r at ω* ± δ and ω* ± iδ, with δ = 1e-3|ω*|.
   - The residual in all 11 equations, ‖M(ω*)x‖ / (‖M(ω*)‖₂ ‖x‖), must be ≤ 1e-10.
3. **Static zero modes and defective roots.**
   - The displacement gauge is removed at ω = 0 by the penalty M(0) + G(0)G(0)†. The remaining null space N is the
     set of static zero modes.
   - A Jordan chain (secular growth ∝ t) exists if N†BN is singular, i.e. its smallest relative singular value is
     ≤ 1e-10. The gauge identity makes the chain condition reduce to this.
   - Where the re-slicing symmetry is exact, one chain per exact direction is expected (gauge). Any excess is
     reported as an unexplained defective chain.
4. **Propagating-mode polarisation.** The transverse-traceless part of h_ij is gauge invariant at linear order. Its
   fraction of the non-gauge part of the metric perturbation is reported. The non-gauge part is the remainder after
   least-squares projection on the continuum gauge orbit, as in `src/exgc34e_review_unstable.py`.
5. **F2 formulation.**
   - Zero-shift gauge (β = 0), a legitimate dynamical gauge. The edge rows read C_ee e + C_en n − A_ee ë = 0, and the
     lapse row reads C_ne e = 0. The lapse has no other terms (no n² term in Regge × lapse; no n in the kinetic term
     at linear order).
   - Index reduction: C_ne ë = 0 gives n(e), so ë = −W e. This is integrated as z = (e, ė) with the implicit midpoint
     rule.
   - The lapse constraint is enforced at the acceleration level: C_ne ė is conserved exactly by the scheme, and C_ne e
     is monitored. The momentum constraint C_βe e + iB_βe ė is not imposed and is monitored.
   - Initial data are eigenvectors of W, or sums of them, so they satisfy both constraints and their time derivatives.
   - Time step dt = 0.02/max|ω| over the modes present, and dt/2. Run length 15/Im ω.
   - Tidal observable: E_ij = ḧ_ij − (q_iq_j/h²) h_00, with h_00 = −2n and h_0i = 0. It is computed from the state
     directly, not from modes.
   - Growth rates are linear fits of ln|·| over the run.
6. **F2(d) replaced.** The original wording mixed normalisations. New (d): run (iii) minus run (ii), at the same dt,
   isolates the seeded component. Its tidal field must grow at Im ω within 1%. The time t× at which it equals the
   wave's tidal field must match t×_pred = ln(|E_wave|/|E_seed(0)|)/Im ω within 1%, where E_seed(0) comes from the
   eigenvectors.
7. **F3 concrete method.** sympy is not attempted. At each growing root, M is split into five blocks:
   kinetic ββ, kinetic βe(+eβ), kinetic ee, Regge ne(+en) and Regge ee.
   - (i) Each block's contribution x†M_X x is reported; the contributions sum to zero at a root. Their scaling with
     kh is reported, with |n| = 1.
   - (ii) Each block's sensitivity d ln ω²/d ln(scale of X) is reported.
   - (iii) The minimal-model test behind the F3 expectation: keep only the blocks with |sensitivity| ≥ 0.05 (the
     others set to zero), recompute the root, and compare with the full root at 1%.
   - Constraint bookkeeping: the nullity of the full 11 × 11 M(ω) at generic ω, the normal-rank deficit n_k of the
     gauge-fixed kernel, and the rank of the lapse and shift rows.
   - The F3 item-3 diagnostic variants are:
     - (a) a per-cube least-squares fit of K from the cube's 19 edges, volume h³;
     - (b) DeWitt trace coefficient λ = 1/2 and 1/3;
     - (c) the lapse–Regge coupling removed (n rows and columns of the static kernel zeroed).
   - For each variant, report: the growing roots (111, generic, 110 at kh = 0.05 and 0.1), the re-slicing breaking
     (34E E2 measure), the constraint counts, whether the 34D static kernel is unchanged, and the E4 gravitomagnetism
     error at kh = 0.05.
