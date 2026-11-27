# EX-GC50 — is there a natural reason for the vacuum subtraction?

2026-11-27 | EXPLORATORY | Claude | PI CLOSED 2026-11-27 | [Plan](../../spec/EX-GC50_plan.md) | About 0.3 core-h on the cloud container (Linux, Python 3.13, numpy 2.x). This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**The clock side has a natural reason, and it is exact. The spatial side does not.**

**The clock side.** Under 44's clock rule, every term of H is multiplied by the local lapse (bonds N̄, sites N). A uniform change of lapse is then only a change of time unit. The vacuum's lapse response has no mass term: χ_NN(q → 0) = 0, and χ_NN ≤ 0 at every q.
- For the pair network, the leftover q² stiffness is −0.0025 (λ₀ = 0.1) and −0.010 (λ₀ = 0.2). That is 400 and 100 times smaller than the +1.00 that 43's kernel gave.
- For the Dirac sea, the q² term goes to zero and only a q⁴ term is left. That is what a relativistic vacuum should do.
- So the clock-sector vacuum terms neither screen G nor need a counterterm.

**The spatial side.** The counterterm that 42–47 needed sits entirely in the spatial ("stress") sector. Stretching a link changes the vacuum's bonding energy. That gives a mixed term χ_Na(0) = −ε_link and a stiffness term χ_aa(0), both of order the vacuum link energy (≈ 3 and ≈ 6 per site for the network, ≈ 6.5 and ≈ 13 for the Dirac sea). These are the lattice form of the cosmological-constant problem: vacuum energy and pressure gravitate.

**43's χ(0) was this term all along.** 43's χ(0) with the 37B weights equals the spatial χ_aa(0) exactly (−5.997975 at λ₀ = 0.1). That follows algebraically, because the 37B coincidence weight N_vN_w is the homogeneous clock weight times a = 1/N̄. So 37B's way of getting light's factor of 2 also pulled the vacuum's stress response into the clock.

**The only cancellation found is exact boson–fermion pairing.** Give each Dirac mode a bosonic partner with the same spectrum (stiffness h²). The total vacuum energy is then zero to 1e-16 for random lapse and strain fields. So every vacuum kernel vanishes, the spatial ones included. That amounts to unbroken lattice supersymmetry. The pair network cannot have such a partner, because its spectrum is not symmetric about zero.

**A consequence for 47.** 47's screening factor 1/(1 + 4πG·1.001) came from 43's kernel. Under the 44 clock rule the lapse-sector factor is about 1, so G_eff is about 1% larger. 47's conclusions (a/a_N ≈ 10⁻⁴) are unchanged.

## Scorecard

| Prediction | Result | Score |
|---|---|---|
| P0 anchor: 37B-weight χ(q) reproduces 43's record within 1e-5 | Within 1.4e-10 at λ₀ = 0.1 and 0.2, for m = 1–8 | PASS |
| P1 network, homogeneous: χ_NN ≤ 0 at every q, \|χ₀\| < 1e-6, c₂ < 0, \|χ_NN(m=1)\| < 1% of the 37B value | Largest χ_NN = −9.7e-5 and −4.0e-4; χ₀ = 1.1e-8 and −3.2e-7; c₂ = −0.0025 and −0.0105; m = 1 ratio 1.6e-5 and 6.8e-5 | PASS |
| P2 network: χ_Na(0) → −ε_link within 1e-5; χ_aa(0) → uniform ∂²E/∂a² within 1e-4 | 4e-9 and 3e-8 (Na); 4e-8 and 8e-8 (aa) | PASS |
| P3 Dirac sea (m = 0.2, 0.5), P1 clauses | χ_NN ≤ 0 everywhere (PASS); c₂ < 0 (−5e-4, −3e-4; PASS); \|χ₀\| = 2.2e-5 and 1.3e-5 > 1e-6 (FAIL) | **P1 χ₀ clause FAIL** |
| P3 Dirac sea, P2 clauses | 2.5e-6 and 1.4e-6 (Na); 1.1e-6 and 5.9e-7 (aa) | PASS |
| C1 supercell vs perturbation, 1e-4 relative or 1e-6 absolute | Largest difference 1.5e-7 absolute (aa), relative ≤ 4e-6 (NN at m = 1, value 1e-4) | PASS |
| C2 boson–fermion sum zero to 1e-10 | \|E_F + E_B\|/\|E_F\| ≤ 1.3e-16 for 5 random (N, a) draws; spectral asymmetry ≤ 3e-14 | PASS |

**Why the Dirac χ₀ clause failed.** The L = 32 fit, over q = 0.2–0.8, sits outside the Dirac sea's analytic small-q region, which needs q ≪ m. Its χ_NN goes as q⁴ there, not q², so a quadratic-in-q̂² fit leaves an intercept of the same size as the data.

At L = 96, post hoc and not scored ([dirac_smallq_posthoc.json](dirac_smallq_posthoc.json)): χ_NN/q̂² falls in proportion to q̂² toward zero, and the same fit gives χ₀ = 6e-8 (m = 0.2) and 2e-8 (m = 0.5), with c₂ = −1.2e-5 and −4.8e-6. So the q² coefficient itself tends to zero. That is the Lorentz-invariant behaviour: the static energy–energy correlator of a relativistic vacuum starts at q⁴.

## Evidence

**Kernels** ([kernels.json](kernels.json)). These come from second-order perturbation theory for the filled sea, at L = 32 and q = 2πm/L. "Network" is the pair network; "37B" means the 37B weights.

| Model | ε_link | χ_NN(m = 1), homogeneous | χ_NN(m = 1), 37B weights | c₂ (clock), homogeneous | c₂ (clock), 37B weights | χ_Na(0) | χ_aa(0) |
|---|---|---|---|---|---|---|---|
| Network λ₀ = 0.1 | −2.9635 | −9.7e-5 | −5.9595 | −0.0025 | +1.0018 | +2.9635 | −5.9980 |
| Network λ₀ = 0.2 | −2.8656 | −4.0e-4 | −5.9311 | −0.0105 | +0.9987 | +2.8656 | −5.9713 |
| Dirac m = 0.2 | −6.5534 | −1.4e-5 | — | → 0 (q⁴) | — | +6.5534 | −13.113 |
| Dirac m = 0.5 | −6.5395 | −9.8e-6 | — | → 0 (q⁴) | — | +6.5395 | −13.107 |

**What the spatial terms would do.** In 44's fields (Φ = N − 1, Ψ = −ln a), the vacuum adds −χ_Na Ψ to ρ_N. The Ψ equation then becomes (∇² + 4πG χ_Na)Ψ = 4πGρ, a Helmholtz equation with an oscillation length of about 2π/√(4πG·3), roughly 36 sites at 4πG = 0.01. There are similar cross terms in the Φ equation. Without a counterterm or a partner cancellation, the far field is not 1/r. This is the problem 42 first ran into, now located in the right sector.

**Checks** ([checks.json](checks.json)).
- C1 diagonalises (32, 1, 1) supercells exactly, with cosine lapse or strain profiles, and takes Richardson energy differences. It agrees with the perturbative kernels for all four entries (NN, Na, aa, and 37B NN) at m = 1, 2 and 4.
- C2 builds the real-space Wilson–Dirac h_{N,a} at L = 4 with random smooth N(x) ∈ [0.8, 1.2] and a(bond) ∈ [0.8, 1.2]. Chiral symmetry (Γ₅ anticommutes with every α and β term) keeps the spectrum symmetric. So the filled-sea energy equals minus the zero-point energy of the partner oscillators, to rounding.

## What this does and does not show

**What it shows:**
- The clock half of 42's problem was an artefact of the 37B weights. The 44 clock rule (time-unit homogeneity) removes it exactly, and leaves G essentially unrenormalised by the vacuum.
- The spatial half is real. It is the vacuum's link (bonding) energy responding to stretching, of the size of the vacuum energy itself.

**What it does not show:**
- Why that spatial term should be absent. A natural cancellation exists only with exact boson–fermion partners.
- That the uniform vacuum energy density (the cosmological constant itself) is addressed. It is not.
- Anything about broken partner symmetry. With unequal partner masses the cancellation would leave terms of order the mass splitting, which was not tested.

**Suggested next (not approved):**
- Test whether the spatial field a can be defined so that the vacuum's link energy does not couple to it at q → 0. For example, a could couple only to departures of the link energy from its local vacuum value, which is homogeneity again, now for the spatial factor. Then check what that does to 44's stress charge Q_a and to light's factor of 2.
- Or test a softly broken partner (bosons with a different mass) to see how the leftover term scales.
