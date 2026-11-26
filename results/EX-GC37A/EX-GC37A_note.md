# EX-GC37A — does gravity's stiffness come from the matter field? (induced gravity on Y1, spin-2 sector)

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC37A_plan.md). PI closure pending. Cloud container, Python 3.13. Under 0.05 core-h. This is our construction, not LaViolette's Model G.

**Partly.** One lattice scalar on Y1 does make spacetime stiff, with the right sign: induced spin-2 gravitons have positive energy. But the stiffness is not GR's:
- the induced gravitons travel at only **0.70–0.78 c** (c_g² = 0.51–0.61 for a light scalar);
- the two polarisations get stiffnesses **2.7× apart**;
- the vacuum produces a graviton mass term about 300–1,800 times larger than the stiffness per unit Q², so it has to be tuned away by hand.

The continuum control, regulated in a Lorentz-invariant way, gives c_g = 1 exactly, as it must. So the cause is Y1's cutoff: a spatial lattice with continuous time. Only stage S1 (the spin-2 sector) was run. Stage S2 (lapse, shift, γ, mode count) was not.

## Method (as planned, S1)

- External transverse-traceless metric h ∝ cos(qz + Ωτ), with q along z.
- Two polarisations:
  - **cross:** h_xy on plaquettes, using 34N's averaged link gradients;
  - **plus:** h_xx = −h_yy at sites.
- Matter couples through √γ γ^{ab} and 1/√γ, expanded to second order. Each h² is placed as an average of squares, so the contact term does not depend on (Ω, q).
- The one-loop Euclidean kernel is contact + bubble, with the frequency integral done analytically.
- Readouts at small Ω, q: A = ∂Π/∂Ω² and B = ∂Π/∂q², giving c_g² = B/A and 1/G_ind = 16πA (EH gives Π = Q²/16πG).

## Results (Brillouin-zone grid 192³; 64³ and 128³ agree to 4 significant figures or better)

| Scalar mass (1/h) | Polarisation | A (kinetic) | c_g² = B/A | G_ind (lattice units) | Mass term / A |
|---:|---|---:|---:|---:|---:|
| 0.1 | cross | 9.49e-4 | **0.607** | 21.0 | 825 |
| 0.1 | plus | 2.59e-3 | **0.512** | 7.7 | 290 |
| 0.5 | cross | 6.93e-4 | 0.531 | 28.7 | 1108 |
| 0.5 | plus | 2.11e-3 | 0.464 | 9.4 | 351 |
| 1.0 | cross | 4.01e-4 | 0.429 | 49.7 | 1809 |
| 1.0 | plus | 1.36e-3 | 0.391 | 14.6 | 515 |

- **P1 (1/G_ind > 0): PASS** for both polarisations at every mass. The induced graviton kinetic term has the healthy sign.
- **P2 (|c_g² − 1| ≤ 2%): FAIL.** The induced gravitons move at 0.70–0.78 c for a light scalar, slower for heavier ones. My recorded expectation was an O(10%) miss. The actual miss is about 40–50% in c_g².
- **The two polarisations differ.** In the continuum they are related by a 45° rotation and are equal. On the lattice their stiffness differs by 2.7×, so the induced spin-2 sector is anisotropic as well as slow.
- **The mass term is enormous.** contact + bubble(0) is about 0.75 in lattice units, against A ≈ 10⁻³. Without a hand-tuned counterterm, the induced graviton would have a mass far above the cutoff. This is the lattice form of the cosmological-constant problem, which the plan subtracted as one recorded tuning.
- **Coefficient fits are stable.** With steps (.02, .04), (.04, .08) and (.1, .2), c_g² changes by less than 1e-4.
- **P3 (γ) and P4 (two modes, λ) are not scored.** S2 was not run. It needs the lapse and shift couplings at second order and their Ward identities, which is the 34S–34Z difficulty. This is a deviation from the plan's two-stage scope.

## Checks

- **G3, exact diagonalisation: PASS.** The static vacuum energy of the lattice Hamiltonian in a cos(qz) background was computed from all normal modes (24³ lattice, m = 0.3, two q, both polarisations). It agrees with the contact + bubble formula on the same grid to 5e-8, the size of the finite-difference error. This also checks the relative normalisation of contact and bubble.
- **G2, continuum Pauli–Villars control: PASS.**
  - Same code with continuum momenta and two Pauli–Villars regulators (M = 3, 6).
  - Π(Ω = s, 0) / Π(0, q = s) = 1.0055, 1.0014 and 1.00035 for momentum cutoffs 60, 120 and 240. That falls as 1/k_max², so the limit is 1.
  - The kinetic term is positive.
  - So the machinery gives Lorentz-invariant induced gravity when the regulator is Lorentz invariant.
- **G1 (Ward identities)** belongs to S2 and was not run. The TT sector needs no gauge fixing at linear order.

## What it means (inference)

- **Sakharov's mechanism is there:** the matter field alone makes space resist bending, with the right sign. That supports the idea that gravity's stiffness need not be a separate ingredient.
- **On Y1 it does not give GR.** The induced coefficient comes from modes near the cutoff. Y1's cutoff is a spatial lattice with continuous time, and its matter modes near the zone edge are far from relativistic. So the gravitons come out slow and anisotropic.
- **That is a design criterion, not a dead end.** If gravity's stiffness is induced by the substrate's matter, the substrate must be locally Lorentz-symmetric right down to its cutoff. Candidates are discrete time treated on the same footing as space, a random or causal-set lattice, or an event network where both are emergent.
- **The giant mass term is the cosmological-constant problem.** No local lattice choice removes it here. Any induced-gravity substrate needs a reason for it to cancel.
- **Comparison with EX-GC37B.** 37B gets the weak-field metric's form from event counting, but it does not explain stiffness, meaning how strongly mass sources N. 37A addresses stiffness and finds it cutoff-dominated. The two are complementary.

## Files

- [Code](../../src/exgc37a.py): S1, G2 and G3 in one file, with separate code paths for G3 (eigenvalues) and G2 (continuum quadrature).
- Results: `s1.json`, `g2.json`, `g3.json`
