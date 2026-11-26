# EX-GC39A — entanglement equilibrium on the lattice (does matter set the geometry?)

2026-11-26 | EXPLORATORY | Claude | PI go-ahead in thread ("Go ahead and do 39A, 39B1/2, and 39C"; "The core-hours are not a concern", 2026-11-26 04:53 UTC). Written before any run. Cap 4 core-h. Design: [EX-GC39_design.md](EX-GC39_design.md). This is our construction, not LaViolette's Model G.

## Question
Jacobson's equilibrium: for every small region at fixed volume, δS_geometry + δS_matter = 0. In the continuum, with δS_geo = η δA (η = the area-law coefficient) and δS_matter = δ⟨K⟩ (first law), this gives ∇²ψ = −4πGρ with G = 1/(4η). Does the same chain work on our lattice?

## Model
- Free scalar on a periodic 16³ cubic lattice, h = 1, mass m = 0.3 (main) and 0.5 (repeat). Gaussian ground states.
- Metric coupling with spatial metric γ_ab = (1 + 2ψ)δ_ab:
  H = Σ_v [½π_v²/√γ_v + ½√γ_v m²φ_v²] + Σ_links ½(√γ γ^{aa})_link (φ_{v+a} − φ_v)²,
  with √γ = (1 + 2ψ)^{3/2}, √γγ^{aa} = (1 + 2ψ)^{1/2}, and the link value taken from the mean ψ of its two ends.
- Ground state: X = ½D^{1/2}S^{-1/2}D^{1/2}, P = ½D^{-1/2}S^{1/2}D^{-1/2}, with D = diag(1/√γ) and S = D^{1/2}VD^{1/2}, from dense eigendecomposition in double precision. Entropy from symplectic eigenvalues (38A `sym_entropy`).
- Blocks: cubes of ℓ = 2, 3, 4, 5 sites, centred on the curvature/matter centre.

## Readouts
1. **Geometry response (η_eff).**
   - ψ = −ε Σ_a cos(2π(x_a − c_a)/16), ε = ±2e-3. Near the centre this gives ∇²ψ = 3ε(2π/16)².
   - Fixed-volume variation: δS_fixV = δS_B[ψ] − δS_B[uniform ψ = ⟨ψ⟩_B].
   - Continuum comparator for a cube of side ℓ with quadratic ψ: δA_fixV = ℓ⁴∇²ψ/3 (derivation in the note), plus the exact first-order value for the actual cos profile, computed on the cube's faces and volume.
   - η_eff(ℓ) = δS_fixV/δA_fixV, against 38A's η = S/A for this mass (z cut: 0.022005 at m = 0.3, 0.019609 at m = 0.5).
2. **Matter response (k_B).**
   - Low-energy "matter" is a first-order, inhomogeneous Gaussian perturbation of the vacuum: a two-mode squeezing that pairs modes k and q − k for |k| ≤ k0 = 2π·2/16, with q = 2π/16 ẑ.
   - It changes X and P by δX, δP of order r = 1e-3, and gives a local energy density δρ(x) ∝ cos(q·(x − c)) at first order.
   - Matter term: δS_B (first law: equal to δ⟨K_B⟩ at first order). k_B = δS_B/δρ(c), with δρ evaluated as the block average of ⟨ℋ_v⟩ (site energy density with half of each link).
3. **Equilibrium G.**
   - η_eff δA_fixV + δS_matter = 0 gives ∇²ψ = −(3k_B/(η_eff ℓ⁴))ρ.
   - So G_eq(ℓ) = 3k_B/(4π η_eff ℓ⁴), with the exact-profile comparator in place of ℓ⁴/3 where it differs.
4. **Coherent matter.** For a displaced (coherent) state, check δS_B = 0 to roundoff while δ⟨ℋ⟩ ≠ 0.

## Predictions
- **P1 (sign):** η_eff > 0 and k_B > 0, so G_eq > 0 (attractive; ψ peaks at the matter). Expectation: uncertain.
- **P2:** η_eff/η = 1 within 10% for ℓ = 3–5. Expectation, recorded: **fails, with η_eff ≪ η.** On a lattice, the cutoff is the site count, so the area-law term counts cut links and is blind to the metric weights; only subleading terms respond.
- **P3:** G_eq(ℓ) independent of ℓ within 20% over ℓ = 3–5 (Jacobson's condition must hold for every small region). Expectation: uncertain.
- **P4 (descriptive):** G_eq against 1/(4η) (38A) and against 37A's G_ind (cross 21.0 at m = 0.1, 28.7 at m = 0.5).
- **P5:** coherent matter gives δS_B = 0 (|δS| < 1e-10), whatever its energy. This is structural. It means classical-field energy cannot gravitate through δS alone.

## Checks
- **G1:** ψ = 0 reproduces the translation-invariant (FFT) correlators to 1e-12, and block entropies agree with 38A's chain-separated method for a slab.
- **G2:** first-order linearity. δS from ±ε (and ±r) agree in magnitude to 1%, so the response is linear.
- **G3, separate path for k_B:** the first-law route. δ⟨K_B⟩ = ½Tr(ℳδX) + ½Tr(𝒩δP), with K_B computed in arbitrary precision from the high-precision flat correlators (k-sum in arb) for ℓ = 2, 3, must equal δS_B to 1%.

## Outputs
`src/exgc39a.py`, `results/EX-GC39A/`.
