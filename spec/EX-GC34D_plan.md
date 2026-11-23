# EX-GC34D — does a lattice-native covariant gravity energy (Regge) give GR light bending? (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "go ahead and close 34A-r, 34B, 34C and proceed with 34D".
Cloud, cap 2 core-hours (expected < 0.5). Fourth stage of EX-GC34. Our construction, not LaViolette's Model G.

## Question

34C found that GR's light bending (PPN γ = 1) follows from two rules together: the lapse multiplies the gravity
field's own energy ("pure clock"), and that energy is one relabelling-invariant functional of the spatial geometry
(√γR in the continuum). The open question was whether a **lattice** can supply rule 2. Regge calculus is the known
lattice construction: geometry is carried by edge lengths, and curvature sits on edges as deficit angles. Its action
Σ_e l_e ε_e is invariant under moving vertices in a flat embedding, which is the lattice version of relabelling
invariance.

34D gives the shared structure a Regge gravity energy multiplied by the lapse. It then solves for the static weak
field of a source and checks γ on the lattice, both directly and with the 34A-r light-bending machinery. It ends with
an analytic look at whether a shift-like variable is forced once the source moves.

## Construction (stated before any γ is computed)

- **Lattice.** Cubic, spacing h = 0.2 (the 34A/34B field lattice). Each cube is split into the 6 Kuhn tetrahedra
  (path 0 → e_σ1 → e_σ1 + e_σ2 → (1,1,1)). This gives 7 edges per vertex: 3 axes, the 3 face diagonals e_i + e_j, and
  the body diagonal. Variables per vertex: lapse n (N = 1 + n) and 7 edge-length changes δl.
- **Gravity energy (pure clock × Regge).** S_g = (1/κ) Σ_edges N_e l_e ε_e, with N_e the mean of the end-point
  lapses and ε_e = 2π − Σ (dihedral angles at e). This is the lattice form of (1/2κ)∫N√γR, using ∫√γR = 2Σ l ε. It is
  written per tetrahedron as Σ_tet Σ_{e∈tet} N_e l_e (θ⁰_e − θ_e(l)), where θ⁰ are the flat angles. Matter couples as
  −Σ_v N_v M_v (the 34B lapse source, energy only). κ = 1.
- **Linearisation.** The per-tetrahedron Hessian in (4 lapses, 6 lengths) is computed at the flat point by
  high-precision finite differences (mpmath), then assembled into the 8 × 8 Fourier kernel K(k). The three vertex
  displacements per vertex are exact null modes of K. They are fixed with the penalty K + G G†, where G(k) holds the
  displacement modes. The source has no overlap with them, so this selects a gauge without changing the physics.
- **Readout of the metric.** h_ij at each vertex comes from the 3 axis and 3 face-diagonal edges, using
  dᵀ h d = 2|d| δl / h. Each edge is moved from its midpoint to the vertex by a Fourier phase. The body diagonal is the
  non-metric extra variable, and its mismatch with dᵀhd is reported.
- **γ, gauge-free.** In Fourier space, γ(k) = −A/(2n̂) with A = ½ P_ij ĥ_ij and P = δ − k̂k̂. Vertex displacements do
  not change A. In real space the gauge parts of h integrate to zero along a straight ray, so the deflection ratio to
  lapse-only is 1 + γ (as in 34C).

Smoke check done before registration, with no γ or field solution computed: the per-tetrahedron Hessian has an exactly
zero lapse–lapse block, the assembled K is Hermitian, the displacement modes are null to finite-difference precision,
and after gauge fixing K has one O(1) eigenvalue (the extra edge mode) and the rest O(k²).

## Parts and registered criteria

**D1 — implementation checks**
- Schläfli identity per tetrahedron, Σ_e l_e ∂θ_e/∂l_f = 0: ≤ 1e-10.
- K Hermitian; lapse–lapse block of K exactly 0 (pure clock); ‖K G‖/‖K‖ ≤ 1e-9 at sampled k, including the zone
  boundary.
- **Separate recomputation path:** the Regge sum is evaluated directly (non-linear dihedral angles, no Fourier
  assembly) on a periodic 8³ lattice. The perturbation is a random small one in all 8 variables at amplitude ±ε. The
  quadratic coefficient [S(ε) + S(−ε) − 2S(0)]/(2ε²) must match ½ xᵀKx from the Fourier kernel to ≤ 1e-4 relative.

**D2 — continuum limit in Fourier space** (registered predictions)
- **P1:** γ(k) → 1. |γ − 1| ≤ 0.01 at kh = 0.05 along (1,0,0), (1,1,0) and (1,1,1), shrinking with k.
- **P2:** the lapse is Newtonian. k² h³ n̂ / M̂ → −κ/2 (that is, G = κ/8π) within 1% at kh = 0.05.
- **P3 (expectation):** exactly one gauge-fixed eigenvalue stays O(1) as k → 0. This is the extra edge mode, and it
  decouples, so it doesn't affect γ.

**D3 — static fields of actual sources (straight-ray deflection, Fourier projection)**
- Sources: a Gaussian (width σ = 1, mass M = 0.4π so that n → −0.05/r, the 34A-r bump amplitude), and the energy
  densities of the 34B rest knots B = 1 and B = 2 (h = 0.2, from the regenerated states; energy only, stress omitted
  as in 34B), each scaled to the same M.
- Periodic box 512 × 256 × 256 sites (102.4 × 51.2 × 51.2). The ray runs along x at impact parameters b = 6, 8, 12.
- **P4:** deflection ratio Regge / lapse-only (same n, flat spatial metric) within 2% of 2 for every source and b.

**D4 — 34A-r light-bending machinery on the solved background**
- The Regge-solved Gaussian field gives N = 1 + n and γ_ij = δ_ij + h_ij at sites and link midpoints (Fourier phase
  shifts), on the 34A-r M2 grid (240 × 100 × 60, lo = (−24, −6, −6), source at the origin). The link weights come from
  A = N√γ γ⁻¹ with the registered split β = (A_ii + A_jj)/12. The control background has the same n and a flat
  spatial metric.
- F1 massless packet, x0 = (−18, 6, 0), k0 = (3, 0, 0), σ = 2, T = 36, dt = 0.05. Readout: circular centroid and
  early (t 2–8) versus late (t ≥ 28) direction, identical to `src/exgc34ar_eval.py`. The rays use the same backgrounds
  via cubic interpolation of the site fields, with the leapfrog dispersion.
- **P5:** packet deflection ratio Regge / lapse-only within 3% of 2. The ray ratio is also within 3% of 2.
- **P6:** |α_field − α_ray| ≤ 5% of |α_ray(lapse-only)| for each background (the 34A-r M2 tolerance).

**D5 — is a shift-like variable forced for a moving source? (analytic, sympy)**
Linearised Einstein equations in (ω, k) for dust moving uniformly at velocity v (ω = k·v). These are compared with
the reduced system that has no shift: only the lapse and h_ij equations, h_0i ≡ 0.
- **P7a (expectation):** for ω ≠ 0 the reduced system is solvable. The momentum constraint then holds automatically
  (Bianchi identity plus conservation), and the gauge-invariant potentials equal GR's.
- **P7b (expectation):** the price is a spatial-metric vector part ∝ 1/ω, which diverges for modes with k ⊥ v. So a
  co-moving, bounded field of a moving body needs a shift-like variable, or an equivalent. Reported with the
  algebra.

## Scope limits

Linear and static apart from D5. The knot source carries energy only (its integrated stress is the 34B lattice
offset). The box is periodic, so the mean source is removed. This affects the Regge and lapse-only backgrounds in the
same proportion, and only ratios are registered. The Kuhn triangulation's edges (axis, e_i + e_j, body diagonal) are
the geometry variables, and the nine field links (both diagonal signs) are built from the metric they imply. Making
the field links themselves the Regge edges is a separate question.

## Files

`src/exgc34d.py` (kernel, checks, Fourier solves, D3; writes `results/EX-GC34D/analysis.json`), `src/exgc34d_bend.py`
(D4 packets and rays; `bend.json`), `src/exgc34d_shift.py` (D5; `shift.json`). Raw packet states go in
`local_raw/EX-GC34D/` (not committed).
