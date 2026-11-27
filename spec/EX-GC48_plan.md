# EX-GC48 — a Dirac-type lattice fermion: does inertia equal energy, and is it dust?

2026-11-27 | EXPLORATORY | Claude | Plan only, pushed before any run. PI go-ahead 2026-11-27 15:52 UTC ("Let's do 7 then 1 then 4"; idea 7 of the GC40-to-now recap = a Dirac-type lattice fermion, the test proposed in the [Codex review of GC47R](../review/EX-GC47R_codex_review.md)). Cap 2 core-h (expected well under 0.5). Our construction; not LaViolette's Model G, not Chain Theory.

## Why

45–47 made matter dust-like by putting its rest energy on a site, but that broke inertia: the 47 pairs have inertia about 8,000 times their energy, so they fall at about 10⁻⁴ of the Newtonian rate. A Dirac-type field puts the rest energy in a mass term that mixes the same components that hop. In the continuum that gives E² = m² + p² (inertia = energy/c²) and, at rest, zero stress (dust) at the same time. 48 checks whether the lattice version keeps both, under the 44 coupling rules, and what happens to a bound composite, where 45–46 found that binding carries stress (no von Laue cancellation).

## Model

Cubic lattice, spacing 1, Hamiltonian form of the Wilson–Dirac fermion (4 components per site):

h(k) = t Σⱼ αⱼ sin kⱼ + β [m + t r Σⱼ (1 − cos kⱼ)],  t = 1 at the reference point, r = 1 (Wilson parameter; r = 0.5 as a sensitivity row).

Positive band E(k) = √(t² S² + M²), S² = Σ sin² kⱼ, M = m + t r W, W = Σ(1 − cos kⱼ). Each band is doubly degenerate (spin).

**44 coupling rules.** Lapse N: homogeneous weights, every term × N, so Q_N = ΔE identically. Space: link terms ∝ 1/a, site terms unchanged, so Q_a = −a ∂E/∂a = t ∂E/∂t (= ⟨H_link⟩). Two readings of the Wilson term:
- **W1 (primary):** the whole Wilson term is a lattice Laplacian (kinetic), so it scales as 1/a. Only m is a site term.
- **W2 (contrast):** strictly by operator support: the Wilson hop −(r/2)β(link) scales as 1/a, its on-site part 3r β does not.

**Signal speed.** c is the low-momentum speed of the massless (m = 0) field in the same lattice and time units, measured independently (48A).

## Parts

**48A. Free particle: inertia vs energy.** For m = 0.05, 0.1, 0.2, 0.3, 0.5 (r = 1 and 0.5):
- E_rest = E(0) = m;
- t_pair = ½ ∂²E/∂k_x² at k = 0, by finite difference (h = 1e-4), along the axis and the diagonal;
- R = 2 t_pair E_rest / c².

c is measured from a massless packet, not assumed. A positive-energy Gaussian at k₀ = 0.05 x̂, σ = 20, L = 256 along x (3D momentum grid, transverse k at 0), evolves exactly in k-space; the centroid speed is fit over t = 0–200.

**Diagnostic 48A-2.** This is Codex's calibration of the 47 bodies. Compute the largest group velocity v_max over all seven bands of the 45/47 network (λ₀ = 0.1, g = 0.025, m₀ = 0.5 and 0.7, L = 64). Then compute R₄₇ = 2 t_pair E / v_max² for B1 and B2 from the 47A records.

**48B. One packet: dust at rest, relativistic stress when moving (W1 and W2).** Positive-band Gaussian packets |φ(k)|² ∝ exp(−|k − k₀|²σ²) on the L = 128 Brillouin zone, with m = 0.2 and 0.5.
- Rest: k₀ = 0, σ = 4, 8, 16, 32.
- Moving: σ = 24, k₀ = 0.05, 0.1, 0.2 along x.

Quantities: Q_N/E, Q_a/E, and v_g² = ⟨(∂E/∂kₓ)²⟩ at the packet.

**48C. Bound composite: does binding energy carry inertia, and does it cancel its own stress?** Two distinguishable Dirac fermions of equal mass m are bound by an on-site contact −U Σₓ n₁(x) n₂(x), with n = ψ†ψ (all components). The contact is a site term, so it does not scale with a. It is projected onto positive-energy states of each fermion (no-pair, as in 45–46).
- **Solver.** With total momentum K, the bound energy solves λ_max[M(E)] = 1, where M(E) = (U/N) Σ_k [Π₁(k+K) ⊗ Π₂(−k)]/(E₁ + E₂ − E). Here Π = ½(1 + h/E) is the positive-band projector, built analytically as a 5 × 5 sum over (1, α₁, α₂, α₃, β) times Dirac-matrix products.
- **Points.** m = 0.2 and 0.5; U chosen so that E_b/(2m) = 0.02 and 0.08; L = 96 (m = 0.5) and 128 (m = 0.2), with a size check at 0.75 L (amended A1); r = 1, W1.
- **Measurements.**
  - E₀ and E_b.
  - The internal kinetic energy T = ⟨E₁ + E₂⟩ − 2m and the potential energy V = E₀ − ⟨E₁ + E₂⟩.
  - t_c from E₀(K) at K = 2π/L and 4π/L, along the axis and the diagonal.
  - R_c = 2 t_c E₀ / c², and δR = R_c / R_free(m) − 1.
  - Q_a(pair) = t ∂E₀/∂t, and η = [Q_a(pair) − 2Q_a(free, k = 0)]/E_b, where Q_a(free, 0) = 0 under W1.
  - A spectator comparator for t_c: the binding-weighted average of the one-body curvatures, with K split equally.

**Checks.** These form a separate path.
- An L = 8 dense band-basis two-body Hamiltonian, built from numerically diagonalised spinors (np.linalg.eigh) rather than the analytic projector. It gives E₀ at K = 0 and K = 2π/8, and Q_a by finite difference, compared with the 5 × 5 solver at L = 8 to 1e-9.
- A real-space check of 48B: build h in real space at L = 24 and project a Gaussian onto the positive band. Compare ⟨H_link⟩ with the k-space Q_a to 1e-10.

## Predictions (scored)

- **P1 (free inertia).** R = 1 + r m to 1e-3 relative at every m (from E² = m² + (1 + r m)k² + O(k⁴)), along both directions. The measured massless c = 1 within 1%. So inertia = energy/c² for fermions much lighter than the cutoff 1/r, with a correction of order m/cutoff.
- **P1-2 (47 calibration, diagnostic).** R₄₇ < 0.05 for B1 and B2: the 47 inertia failure survives an independent signal speed. If R₄₇ falls between 0.5 and 2, I withdraw the 47 "inertia failure" as a units problem.
- **P2 (dust at rest, W1).** Q_N/E = 1 to 1e-12. The rest-packet Q_a/E falls as σ⁻², with a log-slope between σ = 16 and 32 of −2 ± 0.1. At σ = 32 it lies within 5% of 3(1 + r m/2)/(2σ²m²) (amended A1).
- **P2m (moving, W1; amended A1).** [Q_a(k₀) − Q_a(rest, same σ)]/E equals [⟨v_g²⟩ − ⟨v_g²⟩_rest]·(1 + r m/2)/(1 + r m)² within 5% for k₀ ≤ 0.2. In the continuum the trace is ∫T^ii = p·v, which is the r m → 0 limit of this; the lattice factor is the same order-m/cutoff correction as in P1.
- **P2-W2 (contrast).** At rest, Q_a/E ≈ −3r/m (order 1 or larger), so W2 is not dust. The scaling of the Wilson on-site part is decisive.
- **P3 (composite).** Scored at the full L.
  - **(a)** Q_N = E₀ exactly, so gravitational charge counts binding.
  - **(b)** The inertia of the composite misses its potential energy: δR < 0, with δR within 30% of V/E₀ (the instantaneous contact carries no momentum, so only kinetic energy adds inertia). GR needs δR = 0.
  - **(c)** η > 0.5 at every point. Binding kinetic stress is not cancelled by the a-independent contact (no von Laue); GR needs η = 0.
  - **(d)** The spectator formula reproduces t_c within 5%.
- **C1–C2 (checks).** L = 8 dense equals the 5 × 5 solver to 1e-9; real-space ⟨H_link⟩ equals k-space Q_a to 1e-10.

## Amendment A1 (2026-11-27, before any run)

Working the small-k algebra for the code before running anything showed two errors in the plan as pushed. With the W1 rule, t∂E/∂t = k²(1 + r m/2)/E at small k, while the low-k curvature carries (1 + r m), because the r m k² piece of E² scales as t¹, not t². So the rest-packet formula in P2 becomes 3(1 + r m/2)/(2σ²m²), and the P2m comparison carries the factor (1 + r m/2)/(1 + r m)². Separately, the 48C box sizes were too small for the weakest binding (bound-state radius about 40 sites at m = 0.2, E_b/(2m) = 0.01). The binding points move to E_b/(2m) = 0.02 and 0.08 with larger boxes. No data had been computed when these changes were made.

## What the outcomes mean

- **P1 and P2 pass.** A single Dirac-type fermion fixes both 47 problems at once: inertia = energy and dust at rest, under the 44 rules with W1.
- **P3(b) and (c) as predicted.** The remaining problem is the binding interaction, not the matter field. Both binding defects (inertia without potential energy, and stress without von Laue) would trace to an instantaneous site interaction that carries neither momentum nor stress. The next step would bind the fermions through a propagating lattice field whose energy scales with a (a lattice version of idea 6), rather than a contact term.
- **R₄₇ near 1.** The 47 conclusion changes, and I will say so.

## Outputs

`src/exgc48.py`, `src/exgc48_recompute.py`, `results/EX-GC48/` (A.json, B.json, C.json, checks.json, note). LOG row and HANDOFF section.
