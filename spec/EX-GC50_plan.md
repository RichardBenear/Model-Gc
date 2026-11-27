# EX-GC50 — is there a natural reason for the vacuum subtraction?

2026-11-27 | EXPLORATORY | Claude | Plan only, pushed before any run. PI go-ahead 2026-11-27 15:52 UTC ("Let's do 7 then 1 then 4"). Idea 4 of the GC40-to-now recap is a natural reason for 42's vacuum subtraction, in place of a counterterm. Cap 2 core-h, expected under 0.5. Our construction; not LaViolette's Model G, not Chain Theory.

## Why

In 42 and 43 the vacuum's own response to the lapse, χ(q) = δ²E_vac/δN δN, was about −6 at q → 0 (χ(0) ≈ 2ε_link). That turned the Poisson equation into a Helmholtz equation, so 1/r survived only with a local counterterm (43's prescription L, used again in 44 and 47).

Both 42 and 43 used 37B's event weights: links N_vN_w (degree 2 in N) and sites N_v. 44 then introduced the lapse-homogeneous rule, links N̄ = (N_v + N_w)/2 and sites N_v. Under that rule every term of H is linear in N and H(λN) = λH(N). But 44 and 47 kept 43's 37B-weight kernel for the screening factor.

Two exact consequences of homogeneity, argued here before any run:

1. **Euler.** E_vac[λN] = λE_vac[N], so Σ_y N_y ∂²E_vac/∂N_x∂N_y = 0. At N ≡ 1 this gives χ_NN(q = 0) = 0. A uniform change of lapse is only a change of time unit, so it cannot cost vacuum energy at second order.
2. **Concavity.** E_vac is the minimum over states of expectation values that are linear in N, so it is concave. Hence χ_NN(q) ≤ 0 at every q.

If both hold, the clock rule itself is the natural reason for the subtraction: there is no lapse mass term to cancel. What is left is an O(q²) vacuum stiffness, which changes G (anti-screening, by concavity), and a uniform first-order source, the vacuum energy density, which is the cosmological constant and is not addressed here.

The spatial factor a is different. Hopping goes as 1/a and sites do not depend on a, so H is not homogeneous in a. Euler then gives a mixed term χ_Na(0) = −ε_link (minus the vacuum link energy per site) and a nonzero χ_aa(0). In 44's fields (Φ = N − 1, Ψ = −ln a) these would put vacuum mass terms of order the vacuum link energy (its "pressure") into the field equations. The candidate natural reason there is boson–fermion pairing: a partner boson with the same spectrum cancels E_vac identically for every N(x) and a(x).

## Models

- **Network.** The 40A pair network: 6 link modes per site, λ₀ = 0.1 and 0.2, half filling (3 of 6 bands occupied), L = 32.
- **Dirac sea.** The 48 Wilson–Dirac fermion (r = 1, m = 0.2 and 0.5). The Wilson term is written as a sum of bond operators (r/2)β(c_v − c_w)†(c_v − c_w), so under W1 each bond, including its on-site pieces, scales with that bond's 1/a. The negative bands are filled.
- **Lapse weights.** Homogeneous (H): bonds N̄, sites N. 37B (anchor only): bonds N_vN_w, sites N_v.
- **Space.** Bond operators scale as 1/ā, with ā = (a_v + a_w)/2. Sites are unchanged.

## Method

The 2 × 2 vacuum kernel χ_ij(q), i, j ∈ {N, a}, is defined by E_vac = E₀ + ½ Σ δᵀχδ. It is computed from second-order perturbation theory for the filled Fermi sea at q = 2πm/L, with m = 1…8 along the axis and m = (1,1,1) and (2,2,2) along the diagonal. It has two parts:
- **Paramagnetic.** −Σ_k Tr[P_emp(k+q) O_i P_occ(k) O_j†]/Δ, with Bloch transition operators. The bond weight is e^{ik·d}(1 + e^{iq·d})/2 for H and twice that for 37B.
- **Diamagnetic.** This comes from the second-order dependence of the weights: 37B, NN: 2Σ_a e_a cos q_a. H, Na: −½Σ_a e_a(1 + cos q_a). H, aa: Σ_a e_a(1 + cos q_a). Here e_a is the vacuum bond energy per site along axis a.

The uniform reference values come from Euler: χ_NN(0) = 0, χ_Na(0) = −ε_link, and χ_aa(0) = ∂²E_vac(a)/∂a² by a uniform finite difference.

**Separate path (C1).** For the network at λ₀ = 0.1, apply δN = ε cos qx and/or δa = η cos qx on an (L, 1, 1) supercell and diagonalise exactly, with transverse k on a 32² grid. Take χ from the energy second differences, with Richardson extrapolation in ε (ε = 0.01, 0.02), at m = 1, 2 and 4, for both weightings.

**B–F check (C2).** On an L = 4 real-space Wilson–Dirac lattice with random smooth N(x) ∈ [0.8, 1.2] and a(x) ∈ [0.8, 1.2] (5 draws), compare:
- the filled-sea energy E_F = Σ_{ε<0} ε of h_{N,a};
- the zero-point energy E_B = ½ Σ ω of one real bosonic oscillator per fermionic mode, with stiffness h_{N,a}² (ω = |ε|).

## Predictions (scored)

- **P0 (anchor).** With 37B weights, the perturbative χ(q) reproduces 43's recorded χ(q) (λ₀ = 0.1, 0.2; m = 1…8, L = 32) within 1e-5 relative.
- **P1 (no lapse mass under H).** For the network and the Dirac sea, χ_NN(q) ≤ 0 at every computed q. A fit χ_NN = χ₀ + c₂q̂² + c₄q̂⁴ over m = 1…4 has |χ₀| < 1e-6, and c₂ < 0. At m = 1, |χ_NN| is below 1% of the 37B-weight value.
- **P2 (spatial sector keeps vacuum terms).** χ_Na(q → 0) extrapolates to −ε_link within 1e-5 relative, and χ_aa(q → 0) to the uniform ∂²E/∂a² within 1e-4. Both are of order the vacuum link energy, nonzero for the network and the Dirac sea.
- **P3 (Dirac sea).** P1 and P2 hold for m = 0.2 and 0.5.
- **C1.** The supercell energy differences reproduce the perturbative χ_NN, χ_Na and χ_aa at m = 1, 2 and 4, for both weightings, within 1e-4 relative, or 1e-6 absolute where the value is near 0.
- **C2.** E_F + E_B = 0 to 1e-10 relative for every draw. This is algebraic: chiral symmetry makes the spectrum of h_{N,a} symmetric. So with exact partners every vacuum kernel vanishes identically.

## What the outcomes mean

- **P1 passes.** The counterterm of 42/43 was an artefact of the 37B clock weights (degree 2 in N). Under the 44 clock rule the vacuum has no lapse mass, for a structural reason (time-unit homogeneity), and 1/r needs no subtraction. The leftover q² term renormalises G. I will also give the G_eff that 47 should have used, 1/(1 + 4πG c₂_H) in place of 43's 1/(1 + 4πG·1.001), and say whether it matters.
- **P2 passes.** The spatial/stress sector still carries vacuum terms of order the vacuum link energy. These are the lattice form of the cosmological constant problem: vacuum energy and pressure gravitate. The only natural cancellation found here is exact boson–fermion pairing (C2), which needs a boson partner with exactly the same spectrum as each fermion, i.e. unbroken lattice supersymmetry. The 40A network has no such partner, because its spectrum is not symmetric.
- **P1 fails.** Homogeneity is not enough, and I will report what breaks it.

## Outputs

`src/exgc50.py`, `src/exgc50_recompute.py`, `results/EX-GC50/` (kernels.json, checks.json, note). LOG row and HANDOFF section.
