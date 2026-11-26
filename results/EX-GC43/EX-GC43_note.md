# EX-GC43 — what does the lapse couple to? One energy functional, positive-energy excitations, and the vacuum term

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC43_plan.md), written and pushed before the PI's run go-ahead ("go ahead", 2026-11-26 18:23 UTC) and before any run. It follows the [Codex review of GC42](../../review/EX-GC42_codex_review.md). PI closure pending. Cloud container (1 core), Python 3.13, numpy 2.5. About 0.14 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

**Under one consistent functional, gravitational mass is not energy, and it is not universal.** All six predictions and all four checks pass.
- **Free excitations attract, at about twice their energy.** Every positive-energy packet of the unchanged network has positive charge. But Q/ΔE runs from 1.81 to 2.25 for single particles, and 2.40–3.18 for a particle–hole pair. No state anywhere in the zone has ratio 1.
  - The reason is that this network's excitation energy is mostly link (hopping) energy, and link energy counts twice. That is the Tolman factor of light.
- **Hubs gravitate at 0.10–0.41 of their (negative) energy.** That is the same number as their redshift exponent, to 2e-6. So 42A's equivalence-principle failure and this non-universality are one defect.
- **The vacuum term does not cancel by itself.**
  - The raw variational vacuum has χ(0) ≈ −6.0, so the constraint is Helmholtz, not Poisson. The far field is boundary-set: n(r = 10) changes sign as the box grows.
  - Only the local counterterm (or subtracting the vacuum non-locally) restores 1/r. G is then screened by 1/(1 + 4πG c₂χ), with c₂χ ≈ 1.00: twice 42's value, because link energy counts twice here too.

## 43D — positive-energy excitations (the main test)

Packets are Gaussian superpositions of the Bloch states of one band, with no change to H. In this network the two lowest conduction bands (indices 3 and 4) are flat, at e = 1 − λ₀. The top band (index 5) is dispersive, from 1 − λ₀ to 1 + 5λ₀. So the packet shape matters only in the top band.

| Packet | λ₀ = 0.1: ΔE / Q / Q÷ΔE | λ₀ = 0.2: ΔE / Q / Q÷ΔE |
|---|---|---|
| Conduction (flat), any centre, σ = 1, 2, 4 | 0.900 / 1.900 / **2.111** | 0.800 / 1.800 / **2.250** |
| Top band, σ = 1 | 0.987 / 1.962 / **1.988** | 1.037 / 1.882 / **1.815** |
| Top band, σ = 2 | 0.925 / 1.916 / **2.070** | 0.879 / 1.795 / **2.042** |
| Top band, σ = 4 | 0.907 / 1.904 / **2.100** | 0.822 / 1.792 / **2.180** |
| Particle–hole pair, σ = 2 | 1.425 / 3.416 / **2.397** | 0.879 / 2.795 / **3.179** |

(Top-band packets centred at k = 0 and at the band minimum coincide, because the minimum is at k = 0.)

**Band map (P5).** The ratio 1 + e_link/e over the whole zone, excluding |e| < 1e-3:
- Conduction bands: 1.67–2.11 (λ₀ = 0.1) and 1.50–2.25 (λ₀ = 0.2).
- Valence bands: 1.83–3.0 at λ₀ = 0.1. At λ₀ = 0.2 the top valence band reaches e ≈ 0 at its maximum, so the ratio there grows without bound. Here e < 0, so a hole's ratio is set by these values too.
- No state has ratio within 1e-3 of 1. About 120 valence states have ratio 2 (zero site energy). No conduction state does.

**Scoring:**
- **P4: PASS.** Every packet has Q > 0 and ΔE > 0. Conduction packets have Q/ΔE > 1. The spread at fixed λ₀ is 5.9% (λ₀ = 0.1, single particles only) and 19% (including the pair). At λ₀ = 0.2 it is 20% and 62%. All exceed the 5% threshold.
- **P5 (descriptive).** No energy-only state exists. Excitations sit near the light-like value 2, but not at it.

## 43C — hub charges

| λ₀ | Hub | ΔE₀ | Q = Σ_v δq_v | Q / ΔE₀ | Q vs uniform dΔE₀/dN |
|---|---|---:|---:|---:|---:|
| 0.1 | 1 site, Δλ = 0.1 / 0.2 / 0.4 | −0.107 / −0.271 / −0.722 | −0.0104 / −0.0472 / −0.2250 | 0.098 / 0.174 / 0.312 | ≤ 1.8e-6 |
| 0.1 | 7 sites, Δλ = 0.1 / 0.2 / 0.4 | −0.749 / −1.905 / −5.084 | −0.0712 / −0.3282 / −1.6164 | 0.095 / 0.172 / 0.318 | ≤ 1.9e-6 |
| 0.2 | 1 site | −0.166 / −0.377 / −0.894 | −0.0356 / −0.1072 / −0.3503 | 0.214 / 0.284 / 0.392 | ≤ 9e-7 |
| 0.2 | 7 sites | −1.166 / −2.646 / −6.260 | −0.2508 / −0.7701 / −2.5439 | 0.215 / 0.291 / 0.406 | ≤ 6e-7 |

**P3: PASS.**
- Q/ΔE₀ = d ln ΔE₀/d ln N, to 2e-6.
- The spread over hubs is 0.095–0.318 and 0.214–0.406, far above 10%.
- Q < 0 for every hub.

## 43A and 43B — vacuum kernel and far field

| λ₀ | ε_vac | σ_vac = ε′ | χ(0) = ε″ | 2ε_link | χ(0) extrapolated from m = 1, 2 | c₂χ | χ at m = 8 (q = π/2) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | −3.0372 | −6.0007 | −5.9980 | −5.9271 | −5.9980 | 1.0010 | −3.999 |
| 0.2 | −3.1446 | −6.0102 | −5.9713 | −5.7311 | −5.9713 | 0.9962 | −3.990 |

**Far field of a unit point charge** (128³ periodic, fit r = 4–20):

| λ₀ | 4πG | L: G_eff/G (prediction 1/(1 + 4πG c₂χ)) | L: residual | V: 1/r residual | V ÷ L residual |
|---|---:|---|---:|---:|---:|
| 0.1 | 0.01 | 0.99010 (0.99009) | 1.1e-4 | 157 | 1.5e6 |
| 0.1 | 0.05 | 0.95238 (0.95233) | 5.2e-4 | 18.6 | 3.6e4 |
| 0.2 | 0.01 | 0.99013 (0.99014) | 1.0e-4 | 9.0 | 8.6e4 |
| 0.2 | 0.05 | 0.95253 (0.95256) | 5.1e-4 | 2.2 | 4.4e3 |

**Box resonance under V** (n at r = 10, unit charge, λ₀ = 0.1):

| 4πG | Lg = 96 / 112 / 128 / 144 / 160 | Dirichlet ball, R = 30 / 40 / 50 / 60 (k_χR/π) |
|---:|---|---|
| 0.01 | −2.0e-4 / −5.2e-5 / −1.7e-5 / +6.8e-5 / +9.9e-5 | +9.1e-5 (2.33) / +2.1e-4 (3.10) / −6.7e-5 (3.88) / +3.3e-5 (4.65) |
| 0.05 | +1.0e-3 / −2.6e-4 / −1.8e-4 / +1.3e-3 / +1.9e-3 | −1.1e-3 (5.10) / +2.1e-4 (6.80) / −2.2e-4 (8.51) / −6.3e-4 (10.21) |

The λ₀ = 0.2 values are similar (`far.json`). At 4πG = 0.05 they keep one sign but vary 5.5×.

**Scoring:**
- **P1: PASS.** χ(0) < 0 at both λ₀, and χ(0) ≤ 2ε_link (−5.998 ≤ −5.927 and −5.971 ≤ −5.731). V gives no 1/r: its residual is 4e3 to 1.5e6 times larger than L's.
- **P2: PASS.** L is 1/r with residual ≤ 5.2e-4, and G_eff/G matches 1/(1 + 4πG c₂χ) to 6e-5. c₂χ > 0, so the vacuum screens G.
- **P6: PASS (descriptive).** Across the box sweep, n(r = 10) changes sign at three of the four (λ₀, G) settings, and varies 5.5× at the fourth. The Dirichlet amplitude is largest next to k_χR = 3π (R = 40 at 4πG = 0.01) and changes sign between radii. So the raw-vacuum far field is set by the boundary.
  - Codex's point stands: this is a static boundary-value statement. A dynamical instability would need an evolution law, which is not tested here.

## Checks

- **C1 (Hellmann–Feynman, dense L = 8): PASS.** q_v from P equals the finite difference of E to ≤ 3.7e-10 relative, at two sites, for the vacuum and for a Δλ = 0.2 hub.
- **C2: PASS.** χ(q → 0), extrapolated, equals the uniform χ(0) to 1.2e-6 (λ₀ = 0.1) and 3.9e-6 (λ₀ = 0.2).
- **C3: PASS.** The real-space ⟨∂H/∂N⟩ of a σ = 2 conduction packet at L = 16 equals the band formula to 0, and its mean over 200 times equals it to 2e-16. This is trivially exact: a one-band packet has no interband terms, so Q cannot change in time.
- **C4: PASS.** The Woodbury Σ_v δq_v equals the uniform dΔE₀/dN to ≤ 1.9e-6 (gate 1e-4), and ≤ 9e-7 at λ₀ = 0.2 (gate 1e-3).
- **Separate recomputation** (`src/exgc43_recompute.py`). It rebuilds the Bloch H independently and takes E_link = t ∂e/∂t from eigenvalue finite differences, with no eigenvectors. It reproduces every packet ratio to ≤ 3e-11. The shared dependencies are the packet definition and the hopping/on-site split.
- **Deviation.** C3 was planned at L = 16 "dense". It was done with a real-space shift operator at L = 16 rather than a dense matrix, which would be 24,576² entries.

## What it means (inference)

- **The 42 result does not survive a consistent coupling.**
  - When the lapse equation and the matter coupling come from one functional, the charge is Q = ΔE + ΔE_link.
  - Pinned hubs gravitate at 10–40% of their energy, and free excitations at 150–320% (single-particle states 150–225%; the pair up to 318%). No class of states has Q = ΔE.
  - So in this substrate, objects with equal energy fall in each other's field with strengths that differ by factors of order one. An equivalence-principle (Eötvös) failure at order 1, not 1e-13.
- **Why it fails is instructive.** In GR the same bookkeeping happens: active mass density is ρ + 3p, and light has Q = 2E (Tolman). Massive bodies still get Q = E because their internal stresses integrate to zero (von Laue), and that rests on the continuum's dilation and translation symmetry: a stationary self-bound body has zero total stress.
  - This network has no lengths to dilate. Its free fermions cannot self-bind, and its hubs are held by an external coupling.
  - So nothing forces the excess link energy to vanish. The 37B weights give light the right factor 2, but then they give matter roughly the same factor.
- **The vacuum term needs a counterterm.** χ(0) ≈ 2ε_link is set by the bonding energy of the pairs. Nothing in the model cancels it. The local counterterm works, and is the lattice form of tuning the cosmological constant.
- **What would rescue it.** We need a substrate where bound composites exist and obey a lattice von Laue theorem: excess link energy is zero for a stationary self-bound state. Two places to look:
  1. An interacting network that can bind particles: does ⟨H_link⟩_excess → 0 for bound states?
  2. A rule where the lapse weights are chosen so that the conjugate charge equals the energy for any stationary state. The 37B weights do not do this.

  Either is a new design question. I would not scale up this Gaussian network further.

## Files

[Code](../../src/exgc43.py) · [Recompute](../../src/exgc43_recompute.py) · `packets.json` (43D), `hubs.json` (43C), `vacuum.json` (43A kernel), `far.json` (43A/B), `checks.json` (C1, C3), `recompute.json`. `src/exgc42.py` gained an extra returned field (`Slink`, the half-link energy per site); its 42 outputs are unchanged.
