# EX-GC40A — entanglement equilibrium on a pair network

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC40_plan.md). PI closure pending. Cloud container, Python 3.13, numpy 2.5, mpmath. About 0.8 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction (only its axiom, a vacuum made of maximally entangled pairs, is used).

**Counting fixes the geometry side, but the matter side still does not give a constant G.**
- **Pure pairs do not respond to matter at all.** In a vacuum made only of maximally entangled pairs, matter changes no block's entanglement at first order, so it cannot curve this geometry.
- **Area as "cut links" works.** Once pairs share entanglement with their neighbours (λ > 0):
  - every cut link carries the same entropy, wherever it is, to within 3% (all gapped λ);
  - matter now does change entanglement, with the attractive sign.
- **But G is not a constant.** The matter response grows like ℓ^2.4–2.9, not Jacobson's ℓ⁴. So the equilibrium G falls by a factor of 2.3–3 between block sizes 3 and 6.
- Its size, 0.05–0.24, is the same order as 1/(4 s_link) = 0.36, but it is not equal to it.

## Model

- Periodic cubic lattice. Each site carries six fermion modes, one per link end. The pair term joins the two end modes of each link.
- At half filling and λ = 0, every link is one maximally entangled pair. λ couples the six modes of each site.
- The band gap is 2 − 6λ: 2, 1.4, 0.8 and 0.2 for λ = 0, 0.1, 0.2 and 0.3, and it closes at λ = 1/3. **λ = 0.5 is gapless**, so it is reported but kept out of P4.

## Results

**Area (L = 16):**

| λ | s_link = S_B/(cut links), ℓ = 2 → 6 | ℓ² coefficient of a quadratic fit, per link |
|---:|---|---:|
| 0 | 0.693147 at every ℓ (= ln 2 to 9e-14) | 0.693147 |
| 0.1 | 0.6908 → 0.6925 | 0.6933 |
| 0.2 | 0.6851 → 0.6920 | 0.6954 |
| 0.3 | 0.6779 → 0.6940 | 0.7021 |
| 0.5 | 0.6623 → 0.6932 | 0.7094 |

A d = 4 link (two modes per end) gives ln 4 exactly at λ = 0.

**Rewiring at fixed volume (L = 10, dense):** two cut links are rewired into one link inside the block and one outside. The table gives η_count/s_link, with η_count = −δS_B/2.

| λ | min – max over ℓ = 3, 4, 5 and face, edge, corner |
|---:|---|
| 0 | 1.0000 (exact) |
| 0.1 | 0.996 – 1.002 |
| 0.2 | 0.985 – 1.006 |
| 0.3 | 0.971 – 1.012 |
| 0.5 (gapless) | 0.852 – 0.987 |

**Matter (L = 16):**
- The matter is a first-order rotation of 33 valence states into conduction states at k + q, with |k| ≤ 2π·2/16 and q = 2π/16 ẑ.
- The phases are set so that the energy density at the centre is positive.

| λ | k_B, ℓ = 2 / 3 / 4 / 5 / 6 | k_B ∝ ℓ^p, ℓ = 3–6 | G_eq, ℓ = 3 → 6 |
|---:|---|---:|---|
| 0 | δS second order, δρ = 0 at first order | — | — |
| 0.1 | 11.2 / 37.4 / 77.6 / 132 / 200 | 2.42 | 0.159 → 0.053 |
| 0.2 | 14.0 / 53.8 / 120 / 212 / 330 | 2.62 | 0.230 → 0.088 |
| 0.3 | 13.3 / 55.0 / 129 / 238 / 380 | 2.79 | 0.236 → 0.101 |
| 0.5 | 10.4 / 44.8 / 108 / 201 / 326 | 2.87 | 0.195 → 0.087 |

G_eq = 3k_B / (4π s_link ℓ⁴). For comparison, 1/(4 s_link) ≈ 0.36.

## Scoring

- **P1: structural clause PASS, second clause FAIL.**
  - At λ = 0, S_B = ln 2 × (cut links) to 9e-14, and ln 4 × (cut links) for d = 4.
  - "For λ > 0, s_link decreases with λ" holds only for small blocks (ℓ ≤ 4). At ℓ = 6, s_link at λ = 0.3 and 0.5 is at or above ln 2.
  - The fitted area coefficient *rises* with λ, to 0.709 at λ = 0.5. So dressing the pairs adds area-law entanglement; it does not remove it.
- **P2: FAIL as registered, PASS at every gapped λ.**
  - For λ ≤ 0.3, η_count/s_link stays within 3% of 1 at every size and position.
  - At the gapless λ = 0.5 it drops to 0.85 (ℓ = 5, edge). The plan did not exclude gapless λ from P2.
- **P3: PASS.**
  - At λ = 0, δS_B(ε = 1e-3) = 1.2e-6 and δS_B(ε = 1e-4) = 1.5e-8, so δS ∝ ε^1.9. That is second order: |δS|/ε falls about 80× when ε falls 10×.
  - **Deviation, recorded:** the first-order energy density δρ is also exactly zero at λ = 0. The bands are particle–hole symmetric there, and δρ ∝ (e_c + e_v). So k_B is 0/0 at λ = 0, and P3 is scored on δS, its primary clause.
- **P4: sign PASS, size-independence FAIL.**
  - G_eq > 0 (attractive) at every ℓ and λ.
  - Over ℓ = 3–6, G_eq falls by 3.0×, 2.6× and 2.3× at λ = 0.1, 0.2 and 0.3, against a 20% gate.
- **P5 (descriptive):** attractive. G_eq is 0.05–0.24, against 1/(4 s_link) = 0.36. The ratio drifts with ℓ, so no constant can be read off.

## Checks

- **G1: PASS.** FFT and dense correlators agree to 6e-16 at L = 6. At λ = 0, a 3³ block has exactly 54 eigenvalues at 1/2, one per cut link.
- **G2: FAIL at the registered 1%.** |δS(+ε)| and |δS(−ε)| differ by:
  - 39–97% at λ = 0.1;
  - 8–28% at λ = 0.2;
  - 1.6–6% at λ = 0.3;
  - 0.6–3.6% at λ = 0.5.

  The quadratic part is large where the linear part is small (near λ = 0). The central difference used for k_B removes it: central estimates at ε = 1e-3 and 1e-4 agree to 0.26% or better at every row. So k_B is robust, but the gate as registered fails.
- **G3 (first law, separate path): PASS.**
  - δS_B = Tr(K_B δP_B), with K_B from 60-digit (mpmath) Bloch projectors, matches the finite-difference δS_B.
  - Agreement: 7e-7 and 5e-6 (λ = 0.1, ℓ = 2 and 3) and 2e-8 and 1e-8 (λ = 0.3).
  - The block has 10 (ℓ = 2) and 56 (ℓ = 3) exactly pure modes, a structural feature. δP has no weight on them (≤ 3e-18).
  - **Two bugs in the G3 code, recorded:**
    - The first run dropped the on-site λ term on the link matrix element when it rebuilt h(k) in mpmath. It disagreed by 1.2–4×. Record: `g3_first_run_bug.json`.
    - The second run failed on the pure modes (complex logarithm).
    - Neither touched the main results, which come from the double-precision path.
- **G4: PASS at every gapped λ.** The rewired block entropy from a rank-4 Woodbury update of the Green's function (contour integral, no diagonalisation) matches the dense solve to 9e-13 or better. At the gapless λ = 0.5 the Woodbury route has no gap to integrate through and disagrees by O(1), as expected.

## What it means (inference)

- **On the idea that mass comes from pairs sharing entanglement:**
  - Pure pairs are inert to matter: a vacuum of only maximally entangled pairs has no first law.
  - Only once pairs share entanglement with their neighbours does matter register in the entanglement at all. That fits Chain Theory's picture of hubs, and it is what 40C would localise.
- **Counting solves 39A's geometry problem.** Rewiring at a fixed number of sites changes the entropy by exactly s_link per cut link, wherever the change is made. 39A's ambiguity about "fixed volume" is gone.
- **It does not solve the matter problem.** The response grows like ℓ^2.4–2.9, between surface (2) and volume (3), and never reaches Jacobson's ℓ⁴. Two effects contribute:
  - the excitation is a long-wavelength modulation, whose block average falls with ℓ (δρ falls by about 20% from ℓ = 2 to 6);
  - the modular Hamiltonian of a cube on this lattice is not the boost-weighted energy that Jacobson's derivation assumes.
- **Entanglement equilibrium has now failed to give a constant G on three substrates** (39A with metric weights, 40A with counts, and 37A's induced gravity). The "matter sets geometry through entanglement" line looks closed unless the blocks can be made much larger than the cutoff. At ℓ ≤ 6 they are not.

## Files

[Code](../../src/exgc40a.py) · `area.json`, `g1.json`, `rewire.json`, `matter.json`, `g3.json`
