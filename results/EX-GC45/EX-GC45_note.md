# EX-GC45 — heavy on-site matter: can a body's spatial charge be small, and does its binding energy gravitate as energy?

2026-11-27 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC45_plan.md), pushed before the PI's run go-ahead ("ok, go ahead with 45", 2026-11-27 01:37 UTC) and before any run. Follows the [Codex review of GC44](../../review/EX-GC44_codex_review.md). PI closure pending. Cloud container (1 core), Python 3.13, numpy 2.5. About 0.25 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

**Heavy matter can be made dust-like, but its binding energy does not gravitate as energy.**
- **One heavy carrier is nearly dust (45B, P2 PASS).** A rest mode held on a site, reached only through a weak hybridisation g, has a spatial charge Q_a/ΔE of 0.007–0.009 at g = 0.025. It falls as g² (ratio 3.88–4.03 per halving of g), and its mobility (t_eff) falls as g² with it.
  - The catch: Q_a stays a fixed fraction (0.28–0.74) of the band width. So "dust" here means "heavy compared with its own band", the lattice's non-relativistic limit. It is not an exact zero at rest.
- **A bound pair of such carriers carries a binding stress of order its binding energy (45C).** The lattice Nordtvedt ratio η = δQ_a/E_b runs from −2.7 to +4.9 across the nine points, and changes sign with the hybridisation. GR's value is 0. So two bodies of equal energy but different binding fractions would fall differently, at order η·E_b/ΔE.
  - The registered clause "|η| > 0.1 at every point" passes only narrowly (minimum 0.104), and at that point the four-band check (C5) moves η to about −0.04. So η does cross zero somewhere in the (m₀, g) plane, but only by tuning: it is not a rule.
- **The exciton of 44's network stays radiation-like at its true minimum (45A, P1 PASS).** Its minimum is away from K = 0 at every coupling, and there its total Q_a/ΔE is 1.40–2.17. The interaction strain exponent that would cancel the whole charge runs 2.25 → 1.33 and 1.69 → 0.88 across couplings: positive and state-dependent, as Codex found at K = 0.

## 45B — one heavy carrier (λ₀ = 0.1, L = 32)

The heavy band is the 4th of 7 (index 3), with its minimum at k = (π, π, π) for every point.

| m₀ | g | e_min (= ΔE) | width W_h | t_eff | admixture | Q_a/ΔE at min | R = Q_a/W_h | moving k = (π/4,0,0): Q_a/ΔE |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 0.3 | 0.025 | 0.2969 | 0.0078 | 1.86e-4 | 0.0026 | **0.0087** | 0.33 | −0.010 |
| 0.3 | 0.05 | 0.2876 | 0.0307 | 7.42e-4 | 0.0101 | 0.0351 | 0.33 | −0.037 |
| 0.3 | 0.1 | 0.2519 | 0.1171 | 2.96e-3 | 0.0371 | 0.147 | 0.32 | −0.105 |
| 0.3 | 0.2 | 0.1254 | 0.4071 | 1.19e-2 | 0.1127 | 0.899 | 0.28 | −0.098 |
| 0.5 | 0.025 | 0.4963 | 0.0075 | 1.56e-4 | 0.0037 | **0.0075** | 0.50 | −0.001 |
| 0.5 | 0.05 | 0.4852 | 0.0296 | 6.26e-4 | 0.0144 | 0.0296 | 0.49 | −0.003 |
| 0.5 | 0.1 | 0.4432 | 0.1136 | 2.51e-3 | 0.0510 | 0.115 | 0.45 | +0.003 |
| 0.5 | 0.2 | 0.3000 | 0.4000 | 1.02e-2 | 0.1429 | 0.476 | 0.36 | +0.139 |
| 0.7 | 0.025 | 0.6953 | 0.0078 | 1.09e-4 | 0.0058 | **0.0083** | 0.74 | +0.007 |
| 0.7 | 0.05 | 0.6817 | 0.0307 | 4.47e-4 | 0.0219 | 0.0321 | 0.71 | +0.029 |
| 0.7 | 0.1 | 0.6310 | 0.1171 | 1.91e-3 | 0.0736 | 0.117 | 0.63 | +0.121 |
| 0.7 | 0.2 | 0.4675 | 0.4071 | 8.44e-3 | 0.1838 | 0.393 | 0.45 | +0.443 |

- At the band minimum, Q_a equals the link admixture 1 − |⟨h|u⟩|² to the printed digits. That is exact at k = (π, π, π): there the link part of the state sits in the uniform link combination, whose hopping expectation is t·|c_link|², so Q_a = t∂e/∂t = admixture at t = 1.
- A moving carrier's Q_a can be negative (m₀ = 0.3), because away from the minimum the hybridised state mixes into link combinations with negative hopping energy.

**Scoring P2: PASS.**
- Q_a/ΔE at the band minimum falls monotonically as g decreases, at every m₀.
- The ratio between g = 0.05 and g = 0.025 is 4.03, 3.96 and 3.88 (m₀ = 0.3, 0.5, 0.7): g² within 3%.
- At g = 0.025 it is 0.0087, 0.0075 and 0.0083, all below 1e-2.
- t_eff rises by 4.00, 4.00 and 4.09 from g = 0.025 to 0.05, so t_eff ∝ g². Mobility and spatial charge vanish together.
- R = Q_a/W_h (recorded, not predicted) is 0.28–0.74. Q_a tracks the band width, not the rest energy.

## 45C — bound pair of heavy carriers (two species, on-site contact, w = 0, L = 32)

U is set so that E_b = 5 / 15 / 30% of W_h. The pair's dispersion minimum (wedge scan at L = 16) is at K = 0 for every point. With both carriers at (π, π, π), that is the free pair's minimum too.

| (m₀, g) | E_b/W_h | E_b | ΔE_pair | E_b/ΔE_pair | Q_a,pair/ΔE_pair | δQ_a | **η = δQ_a/E_b** | ⟨T⟩ | δQ_a/⟨T⟩ | w* = −δQ_a/⟨V⟩ (post hoc) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| (0.3, 0.2) | 0.05 | 0.0204 | 0.2305 | 0.088 | 0.742 | −0.0545 | **−2.68** | 0.0776 | −0.70 | −0.56 |
| | 0.15 | 0.0611 | 0.1897 | 0.322 | 0.645 | −0.1030 | **−1.69** | 0.1115 | −0.92 | −0.60 |
| | 0.30 | 0.1221 | 0.1287 | 0.949 | 0.509 | −0.1599 | **−1.31** | 0.1344 | −1.19 | −0.62 |
| (0.5, 0.1) | 0.05 | 0.0057 | 0.8808 | 0.0064 | 0.118 | +0.0018 | **+0.311** | 0.0191 | +0.093 | +0.07 |
| | 0.15 | 0.0170 | 0.8694 | 0.0196 | 0.115 | −0.0018 | **−0.104** | 0.0273 | −0.065 | −0.04 |
| | 0.30 | 0.0341 | 0.8524 | 0.040 | 0.110 | −0.0079 | **−0.233** | 0.0328 | −0.24 | −0.12 |
| (0.7, 0.05) | 0.05 | 0.0015 | 1.3618 | 0.0011 | 0.0377 | +0.0075 | **+4.87** | 0.0040 | +1.88 | +1.36 |
| | 0.15 | 0.0046 | 1.3587 | 0.0034 | 0.0392 | +0.0095 | **+2.06** | 0.0056 | +1.68 | +0.93 |
| | 0.30 | 0.0092 | 1.3541 | 0.0068 | 0.0398 | +0.0100 | **+1.09** | 0.0068 | +1.48 | +0.63 |

- The single carrier's Q_a/ΔE at these points is 0.899, 0.115 and 0.032.
- The strongly hybridised pair (g = 0.2) has *less* spatial charge than its two free carriers: binding localises it and suppresses hopping. The heavy pair (g = 0.05) has *more*: binding admixes link modes that the free carriers at the band minimum barely use.
- The last column is the strain exponent of the interaction (U ∝ a^−w) that would make the binding stress-free, as in 44. It is computed from recorded δQ_a and ⟨V⟩, after the data were seen. It ranges from −0.62 to +1.36 and is not universal.

**Scoring P3.**
- **|η| > 0.1 at every point: PASS, narrowly.** The smallest is 0.104, at (0.5, 0.1, 15%). That point is also where C5 shows the largest relative projection error (below): with all four upper bands, η there is about −0.04. So this clause is not robust at that one point.
- **Spread of η > 20%: PASS.** η runs from −2.68 to +4.87 and changes sign across (m₀, g). At fixed (m₀, g) it also falls in magnitude as E_b grows.
- **Surprise clause (|η| < 0.01 everywhere): not triggered.** Binding is not stress-free on this lattice.

**P4 (descriptive).** δQ_a/⟨T⟩ runs from −1.19 to +1.88. Binding stress is of the order of the pair's internal kinetic energy, with a sign set by the hybridisation. There is no clean continuum value for a contact pair to compare with.

## 45A — exciton of 44's network at its true dispersion minimum

U is 44's value for each E_b fraction (fixed at K = 0). K* comes from a wedge scan over the full zone at L = 24; the charge is evaluated at K* on L = 48. E_b(K*) is measured from the continuum minimum at that K.

| λ₀ | E_b/gap at K = 0 | K* (units 2π/24) | E₀(K*) | E₀(0) | E_b(K*) | ⟨V⟩ | total Q_a (w = 0) | binding-only Q_a | **Q_a/ΔE** | Q_a,x / Q_a,y / Q_a,z | w cancelling the total |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|
| 0.1 | 0.05 | (12, 12, 12) | 1.2784 | 1.3300 | 0.722 | −0.793 | 1.786 | −0.214 | **1.397** | 0.604 / 0.604 / 0.604 | 2.25 |
| 0.1 | 0.15 | (12, 12, 12) | 1.0772 | 1.1900 | 0.923 | −0.994 | 1.750 | −0.250 | **1.625** | 0.594 / 0.594 / 0.594 | 1.76 |
| 0.1 | 0.30 | (12, 12, 12) | 0.7937 | 0.9800 | 1.206 | −1.278 | 1.699 | −0.301 | **2.141** | 0.580 / 0.580 / 0.580 | 1.33 |
| 0.2 | 0.05 | (1, 1, 1) | 0.7593 | 0.7600 | 0.066 | −0.902 | 1.520 | −0.443 | **2.002** | 0.516 / 0.504 / 0.503 | 1.69 |
| 0.2 | 0.15 | (0, 0, 3) | 0.6779 | 0.6800 | 0.195 | −1.099 | 1.394 | −0.509 | **2.056** | 0.456 / 0.432 / 0.504 | 1.27 |
| 0.2 | 0.30 | (0, 0, 6) | 0.5529 | 0.5600 | 0.500 | −1.358 | 1.199 | −0.543 | **2.168** | 0.343 / 0.323 / 0.532 | 0.88 |

- At λ₀ = 0.1 the minimum is at the zone corner, K = (π, π, π), and sits well below the K = 0 value (by 0.05–0.19). There the strain response is isotropic.
- At λ₀ = 0.2 the minimum is near K = 0 and only 0.0007–0.007 below it. Its directional charges are unequal (z largest along K*), so a moving composite's stress is anisotropic, as a moving body's should be. But the rest-frame value is not isotropic-zero either.
- (The wedge scan orders n_x ≤ n_y ≤ n_z, so (0, 0, 3) stands for its cubic images; the evaluation uses the listed one.)

**Scoring P1: PASS.**
- K* ≠ 0 at all six couplings.
- Q_a/ΔE > 1 at every coupling and both λ₀ (1.40–2.17).
- The total-cancelling w varies by 41% (λ₀ = 0.1) and 48% (λ₀ = 0.2) across couplings, both above 20%.

## Checks

- **C1 (homogeneity): PASS.** Q_N = ΔE for one carrier (0.4432235637 against 0.4432235637, 1.4e-10 relative) and for every pair (≤ 1.3e-11 relative).
- **C2 (g = 0): PASS.** The heavy band is flat at m₀ (spread 0) and Q_a = 0 exactly.
- **C3 (dense vs separable): PASS.** At L = 10, (0.5, 0.1, 15%): 0.8694141941703644 against 0.8694141941703587, difference 6e-15.
- **C4 (Hellmann–Feynman for V): PASS.** U∂E₀/∂U equals ⟨V⟩ to ≤ 8e-9 relative (gate 1e-8; largest at (0.7, 0.05)).
- **C5 (projection error): FAIL at its registered point.** At (0.5, 0.1, 15%), the hopping derivative with both carriers allowed in all four bands above μ differs from the heavy-band value by 57% (L = 4) and 62% (L = 6) of δQ_a. The gate was 5%.
  - The absolute difference is small (0.0011, or 0.065 E_b). It is large relative to δQ_a only because δQ_a is near its sign change there.
  - **Extension (post hoc, `C5_extended.json`).** The same comparison at L = 6 for all nine points. The relative error is 0.5–1.8% at (0.7, 0.05), 3–13% at (0.3, 0.2) and 28–62% at (0.5, 0.1). As a shift in η it is 0.02 (heavy points), 0.08–0.17 (strongly hybridised) and 0.065–0.08 (middle).
  - So the sign and order-one size of η stand at (0.7, 0.05) and (0.3, 0.2). At (0.5, 0.1), η is small and its sign is not fixed by the projected problem.
- **Separate recomputation (`src/exgc45_recompute.py`): PASS.**
  - 45B: an independently built 7-band H(k) (opposite phase convention), eigenvalues only, a 4-point stencil in t, at L = 16. It reproduces every Q_a/ΔE at the minimum to ≤ 8.5e-9 relative.
  - 45C: a dense projected two-body matrix at L = 10, with its own eigenvectors and pair basis, gives η = −1.68713, −0.10420 and +2.06207 at the three 15% points. The driver's separable solver at the same L agrees to ≤ 3e-7. The L = 32 values differ from L = 10 by ≤ 5e-5, so these bound pairs are compact.
  - Shared with the driver: the model definition and the U values.

## Deviations

- **C5 size.** The plan's "65,536 pair states" counted all two-carrier states at L = 4. The check uses total momentum K = 0, which is exact by translation symmetry: 1,024 states at L = 4 and 13,824 at L = 6, dense. L = 6 was added.
- **C5 extended** to all nine points after the registered point failed (post hoc, labelled).
- **45A evaluation size** was L = 48, not the planned L = 32 (finer K grid; the scan at L = 24 maps onto it exactly).
- **w* column in 45C** is post hoc.

## What it means (inference)

- **The rest-energy problem has a lattice answer, with a price.** Putting the rest energy on a site does make a carrier's spatial charge small, as g². But it is small only relative to the rest energy. It stays a fixed fraction of the band width. On this lattice "at rest" still means "in the lowest Bloch state", which hops. So dust is reached only in the limit where the body is heavy compared with its own kinetic scale, the lattice version of the non-relativistic limit. That is acceptable: GR's ρ ≫ p is the same statement.
- **Binding energy does not gravitate as energy.** For a bound pair, δQ_a is of the order of E_b, and its sign depends on how the binding reshapes hopping. In GR the total stress of a stationary bound body is zero (von Laue), so binding energy weighs exactly what it is (Nordtvedt parameter 0; lunar ranging bounds it near 1e-4). Here η is order 1, so the strong equivalence principle fails at order (E_b/ΔE).
- **Why (inference).** A contact interaction with w = 0 carries no stress by construction. So δQ_a is pure kinetic stress, and nothing on the lattice balances it. In the continuum, the balance comes from the interaction's own stress, which is fixed by how the force law scales with distance (the virial theorem). On a lattice with on-site contact, there is no distance for the force to scale with. The exponent w* that would balance it varies from −0.6 to +1.4 between points, so no single strain law for the interaction fixes it.
- **Where this leaves the programme.** 43–45 together say: a fixed lattice with lapse-weighted hopping gets the clock right (Q_N = ΔE) and can make near-dust matter, but cannot make bound bodies stress-free. A von Laue identity needs an interaction whose stress balances the kinetic stress for every stationary state. In the continuum that follows from translation and dilation invariance. A substrate that is to give it must have a dilation-like symmetry in its interaction law: for example, interactions that act across links with strengths set by the same a that sets hopping, so that one scaling acts on both.
