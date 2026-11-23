# EX-GC34H — constraint-reduced analysis of the cube-mean kinetic term B (EXPLORATORY)

Claude, 2026-11-23. Pre-run edit (before any computation): the zero-root handling in H3 and H4 is now defined by a
zero cluster, because defective zero roots split numerically. PI go-ahead (Richard, 2026-11-23): "go ahead with the B run first." This plan and the
`explore/LOG.md` entry are pushed before any 34H computation. Cap: 2 core-hours on the cloud; I expect under 0.1. Our
construction, not LaViolette's Model G. It follows priority 1 of ChatGPT's review
(`review/EX-GC34F-34G_chatgpt_review.md`).

## Question

In 34G, candidate B (per-cube mean of the six tetrahedra's K, V = h³) had no genuine growing root and no third
low-frequency mode, post hoc. It failed three things:
- a static Jordan chain in every direction;
- a strongly anisotropic polarisation;
- an inconclusive time-domain look.

The suspected cause is one edge combination per vertex with no kinetic energy. 34H asks four things:
- (a) What is that combination, and what constraint does it obey?
- (b) Once it and the lapse are eliminated, is the remaining system a regular second-order system with exactly two
  physical waves and no growth?
- (c) Is the zero-frequency (drift) sector pure gauge, meaning no spacetime curvature?
- (d) Does a constraint-reduced time evolution with a direct tidal observable behave?

The wave-speed error is a separate gate. The 34G result is not re-scored.

## Setup

B exactly as in `src/exgc34g.py` (`kin_B`), with the unchanged 34D lapse × Regge potential.

**Variables.** Zero-shift gauge (β = 0), unknowns x = (n, 7 edges). There, the kinetic term is ω² A_ee only, so the
time-domain equations are C x − A ẍ = 0 on the (n, edge) rows. The momentum rows are β: C_βx + iB_βe ė = 0.

**Reduction.**
- u(q) is the unit null vector of A_ee. Q(q) (7×6) is an orthonormal basis of its complement.
- Write the edges as Q y + a u.
- The lapse row and the u-projected edge row contain no second time derivatives. They form a 2×2 linear system for
  (n, a) in terms of y.
- If that block is non-singular, the remainder is the regular system ÿ = F y with
  F = (Q†AQ)⁻¹ Q†C x(y), where x(y) = (n(y), Q y + a(y) u).
- **Fallback, fixed now.** Where the 2×2 block is singular (σ_min/σ_max < 1e-10), the case is reported as singular
  and not reduced.

**Cases.** The 34G direction set, i.e. the 34F set:
- axes, signed face and body diagonals, generic, near-axis, the 30 random directions, and the approach sequences;
- kh = 0.0125, 0.025, 0.05, 0.1, 0.2;
- along 100, 110 and 111 also kh = 0.5, 1, 2, 3.

The light-speed reference is `exgc34e.light_speed(q)`.

## Criteria (thresholds fixed now)

- **H1: a single non-kinetic edge combination.** In every case A_ee has exactly one singular value below 1e-12 of its
  largest, and the next is above 1e-6.
  - Reported: u at small kh, as a combination of the 7 edge types, and its limit as kh → 0.
  - Reported: the coupling of u to the lapse, |C_n,e u| relative to |C_n,e|.
- **H2: the reduction is regular.** The 2×2 (n, a) block is non-singular in every case (σ_min/σ_max ≥ 1e-10).
- **H3: mode count and no growth.** Eigenvalues of F give ω² = −λ.
  - A reduced solution is physical if it also satisfies the momentum rows: relative residual ≤ 1e-8, using the full
    11-variable M(ω) with the reconstructed x.
  - The zero cluster is the eigenvalues with |λ| ≤ 1e-6 ‖F‖₂. Defective zero eigenvalues split numerically, so the
    cluster is analysed only through H4. Its size is reported against the generalised-null-space dimension of L.
  - Pass if, in every case:
    - no eigenvalue outside the zero cluster has |Im ω| > 1e-6 |ω|;
    - at kh ≤ 0.05, exactly two physical roots with ω > 0 and 0 < ωh < 0.5, counted with multiplicity;
    - every other physical root outside the zero cluster has ωh ≥ 0.5.
  - Internal check: each physical root satisfies the full 11-variable equation with residual ≤ 1e-10.
- **H4: the zero-frequency sector is pure gauge.**
  - Take the invariant subspace Z of the first-order form L = [[0, I], [F, 0]] for its eigenvalues with
    |μ| ≤ (1e-6 ‖F‖₂)^½, from an ordered Schur form. Solutions there are z(t) = Z exp(N t) c with N = Z†LZ.
  - Keep the solutions that satisfy the momentum rows at t = 0, 1/k and 2/k (singular value ≤ 1e-8 of the stacked
    map's scale).
  - For each kept solution, form the full spacetime metric (n, β = 0, edges via the 34D map) and its linearised
    Riemann tensor at t = 0 and t = 1/k.
  - The measure is r = max_t |R| / (k² max_t |h|).
  - The control r_disp is the same measure for the three exact displacement-gauge modes at the same q.
  - Pass if every kept solution has r ≤ max(10 r_disp, 1e-10) in every case.
  - Reported: the kh scaling of r and r_disp, and r for a physical wave (expected O(1)).
- **H5: reduced time evolution.** Generic and (1,1,0) directions at kh = 0.05.
  - Run the regular system ÿ = F y with the implicit midpoint rule, steps dt = 0.02/ω_max and dt/2, for 40 periods of
    the slower physical wave.
  - Initial data: that wave (unit amplitude), plus every other eigen- or generalised-eigen component of L seeded at
    1e-8. The seeded components include the zero sector, so they may drift linearly.
  - The observable is the tidal field E_ij = ½(ḧ_ij − q_iq_j h00/h²), from the reconstructed n and edges.
  - Pass if:
    - the momentum-row residual is ≤ 1e-8 relative throughout;
    - max_t |E(t) − E_wave(t)| ≤ 1e-6 max_t |E_wave(t)|, where E_wave is the same run with the seeds off;
    - the two steps agree on that measure to within a factor 2.
- **H6: wave speed (separate gate).** At kh ≤ 0.05 in every direction, both physical waves satisfy
  |ω/k − c| ≤ 1e-3.
  - Reported: (ω/k − c) against kh, a fit of a(dir)(kh)², and the overlap of each polarisation with u.
- **C1: control.** The same classifier, applied to the per-tetrahedron term T in the 34F reduction (n from the twice
  differentiated lapse constraint, `exgc34f.reduced_system`), must flag growth in the generic direction at kh = 0.05.
  - Its growth rate must match the confirmed 34F/34G root to ≤ 1e-6 relative.
  - If C1 fails, H3 is not scored.
- **S1: separate path.** Every non-zero physical root of the reduced system must match a confirmed root of B in 34G's
  `candidates.json` to ≤ 1e-6 relative, at the cases recorded there.
  - Shared dependencies: the B kernel and the 34D potential.

**Verdict.**
- **B is a viable kinetic term at linear order** if H1–H5 pass, C1 passes and S1 passes.
- H6 is reported as its own result. A pass on H1–H5 with an H6 failure means "stable, but not isotropic".

**My expectations, stated before running:**
- H1 passes.
- H2: I guess it passes except possibly on the axes, but I'm not confident.
- H3 passes (consistent with 34G's post hoc reading).
- H4: no firm prediction. I lean towards gauge up to lattice error.
- H5 passes if H2–H4 do.
- H6 fails (34G: one polarisation 5.8% fast at kh = 0.05).

## Scope

Linear, flat background, continuous time. A pass does not make B physically preferred. Non-linear and coupled-matter
tests would come later.

## Files

`src/exgc34h.py` → `results/EX-GC34H/reduced.json` (H1–H4, H6, C1, S1) and `evolve.json` (H5). It reuses
`src/exgc34e.py`, `exgc34f.py` and `exgc34g.py` unchanged.
