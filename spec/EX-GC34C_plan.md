# EX-GC34C — what selects the gravity field's own equations, and hence light bending? (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "go ahead with 34C". Cloud, cap 2 core-hours (expected
< 0.1: analysis plus Fourier-space solves). Third stage of EX-GC34. Our construction, not LaViolette's Model G.

## Question

34B showed that matter on the shared structure sources the lapse n with its energy density ρ and the spatial metric
h_ij with its stress T^{ij}, and that ∫T^{ij} = 0 for a body at rest. Light bending (PPN γ) then depends on the
gravity field's own action L_g[n, h]. Which physically motivated rule for L_g gives which γ? Can a rule that belongs
to the shared-structure idea select γ = 1 (GR), rather than imposing GR?

## Candidate rule (stated before the analysis)

**"Pure clock": the lapse is a clock rate for gravity too.** Every matter field's energy is (static energy) × N +
(kinetic)/N, with no energy stored in gradients of N. If the gravity field obeys the same rule, then N appears in
L_g only multiplicatively. In the static weak field, L_g = −(1/κ)[(1 + n) Λ(h) + Q(h)], with Λ linear and Q quadratic
in ∂h, and no (∂n)² term. The n-equation then becomes a *constraint* that ties the spatial metric to the energy
density. With an independent (∂n)² term it would instead be a separate Poisson equation for n.

## Actions compared (static, linear, κ = 1, source ρ = M × normalised Gaussian of width 0.5, T^{ij} = 0)

| id | action | origin |
|---|---|---|
| A0 | −½(∂n)² − Q₀(h) with Q₀ = ¼ Σ(∂h_ij)², no n–h coupling | "each part obeys its own equation" (34B's simplest) |
| A1 | Nordström: n = φ, h_ij = 2φδ_ij, −c(∂φ)² | conformally flat scalar gravity |
| A2 | pure clock, general: Λ = a ∂_i∂_j h_ij + b ∇²h_kk; Q = c₁ ∂_k h_ij ∂_k h_ij + c₂ ∂_i h_ij ∂_k h_kj + c₃ ∂_i h_ij ∂_j h_kk + c₄ ∂_k h_ll ∂_k h_mm | the rule above, with nothing else imposed |
| A3 | A2 at the coefficients of the static Fierz–Pauli action (linearised GR): a = −b, and (c₁..c₄) ∝ (−½, 1, −1, ½) in FP normalisation | the rule plus spatial covariance (invariance under h_ij → h_ij + ∂_iξ_j + ∂_jξ_i) |

The matter coupling is the same in every case and comes from 34B: −nρ + ½h_ij T^{ij}.

## Method (`src/exgc34c.py`)

1. **Symbolic, Fourier space (sympy).** Decompose the static scalar sector as h_ij = A δ_ij + B k̂_i k̂_j. Derive
   the field equations from each action and solve for n̂, Â, B̂ in terms of ρ̂.
2. **γ in a gauge-free way.** Light deflection depends on n and A only, because the B part integrates to zero along
   a straight ray. Define γ = −Â/(2n̂) at k → 0. In isotropic gauge this is the PPN γ, and deflection ∝ (1 + γ).
3. **General A2.** Find the conditions on (a, b, c₁..c₄) for a solvable static problem and for γ = 1. Report whether
   the pure-clock rule alone fixes γ.
4. **Numerical check (separate path).** For A0, A1, A3 (and one generic A2 point), build n and h_ij on a 3D grid by
   inverse FFT of the solved spectra (128³, box 64). Compute the straight-ray deflection integral
   α(b) = ∫ ∇_⊥(−n + ½ h_ll) dl numerically at impact parameters b = 8, 12, 16. Compare with (1 + γ) × the
   lapse-only value from the same n.

## Registered predictions

- **P1:** A0 gives γ = 0 (lapse-only bending).
- **P2:** A1 gives γ = −1 (no bending of light).
- **P3:** A3 gives γ = 1, with n = −κM/(8πr) far away. In those units this is Newtonian with G = κ/8π.
- **P4 (my expectation, could fail):** the pure-clock rule alone does *not* fix γ. Within A2, γ depends on the
  coefficients, and γ = 1 needs the covariance relations (or a weaker condition the analysis will identify).
- **Numerical check:** the deflection integrals agree with (1 + γ) × lapse-only to within 2% for each case at
  b = 12 and 16 (b = 8 may feel the Gaussian core).

## Scope limits

Static, linear, weak field only. This decides γ, not the full dynamics. No shift vector exists in the structure, so
the radiative and moving-source sector (gravitomagnetism, gravitational waves) is outside this stage. If P4 holds,
the next question (not part of 34C) is whether something native to the lattice can supply spatial covariance, for
example invariance of the gravity energy under relabelling or re-meshing of sites.
