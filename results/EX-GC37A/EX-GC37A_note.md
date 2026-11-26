# EX-GC37A — does gravity's stiffness come from the matter field? (induced gravity on Y1, spin-2 sector)

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC37A_plan.md). **PI CLOSED 2026-11-26** (Richard: "close 37A/B"). Cloud container, Python 3.13. Under 0.05 core-h. This is our construction, not LaViolette's Model G.

**Partly.** One lattice scalar on Y1 does make spacetime stiff, with the right sign: induced spin-2 gravitons have positive energy. But the stiffness is not GR's:
- the induced gravitons travel at only **0.78–0.85 c** (c_g² = 0.61–0.72 for a light scalar; corrected 2026-11-26, first run said 0.70–0.78 c);
- the two polarisations get stiffnesses **2.1× apart** (corrected; first run said 2.7×);
- the vacuum produces a graviton mass term about 300–1,800 times larger than the stiffness per unit Q², so it has to be tuned away by hand.

The continuum control, regulated in a Lorentz-invariant way, gives c_g = 1 exactly, as it must. So the cause is Y1's cutoff: a spatial lattice with continuous time. Only stage S1 (the spin-2 sector) was run. Stage S2 (lapse, shift, γ, mode count) was not.

> **Correction 2026-11-26 (after [ChatGPT's review](../../review/EX-GC37AB_chatgpt_review.md)).** The plus-polarisation vertex in the first run was the plain G_a² form, not the 34O action the plan specified. 34O's normal-stress coupling is half plain square and half forward–backward product, which gives K_a² cos²(k_a/2). Rerun with the 34O vertex and a matched contact term (G3 exact diagonalisation agrees to 6e-8):
>
> | Mass | plus (34O vertex): A | c_g² | G_ind | cross (unchanged): c_g² |
> |---:|---:|---:|---:|---:|
> | 0.1 | 4.48e-4 | **0.716** | 44.4 | 0.607 |
> | 0.5 | 2.78e-4 | 0.621 | 71.6 | 0.531 |
> | 1.0 | 1.34e-4 | 0.488 | 148 | 0.429 |
>
> - This reproduces ChatGPT's 0.716.
> - P1 still passes, and P2 still fails: induced gravitons move at 0.78–0.85 c for a light scalar.
> - The polarisation split is now 2.1×, in the reverse order (cross is the stiffer one). The 2.7× in the summary and the "plus" rows below belong to the superseded vertex and are kept for the record.
> - "Mass term" below means the zero-momentum TT response. Calling it a physical graviton mass would need the full constraint and gauge sector (S2), as the review notes.
> - Rows: `s1b.json`; G3 rows `plus34o` in `g3.json`.

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

## Stage S2, static scalar sector (amendment A2, run 2026-11-26)

**P3 FAIL.** In the static scalar sector the induced kernel is not GR's:
- in the gauge h_L = 0 the light-bending parameter comes out **γ = 0.236**;
- the kernel strongly violates the static displacement-gauge identity, so that number depends on the gauge;
- one GR feature does come out right: the induced lapse stays non-dynamical.

**Static kernel, q² coefficients** (m = 0.5; 32³ and 48³ agree to 4 digits):

| | nn | nT | nL | TT | TL | LL |
|---|---:|---:|---:|---:|---:|---:|
| b (q² coefficient) | −9e-6 | 0.0433 | 0.0219 | 0.0918 | 0.0219 | −0.0010 |

- **γ = Π_Tn / (2Π_TT) = 0.236** with h_L = 0, sourced through n. The registered gate was within 2% of 1.
- **The same formula on the 34N Y1 Einstein–Hilbert kernel gives exactly 1** (control), and that kernel's nn, nL, TL and LL entries vanish.
- **Ward residuals** (0 if the induced action were invariant under static displacements):
  - Π_TL/Π_TT = 0.24;
  - Π_nL/Π_nT = 0.51;
  - Π_LL/Π_TT = −0.011.
- **Without a gauge choice the result is not even isotropic.** If h_L is treated as physical, the solve gives h_T = −0.92n and h_L = +1.86n, i.e. "γ_T" 0.46 and "γ_L" −0.93. Not a metric GR would recognise.
- **What does come out right:** Π_nn has essentially no q² term (nn/nT = −2e-4) and is exactly 0 at q = 0. The lapse stays a Lagrange-multiplier-like field, as in GR. That is structural, because the Hamiltonian is linear in N.
- **Separate path:** full real-space 6³ diagonalisation agrees with the block method to 2e-16.
- **Caveat on model dependence.** The linear vertices are 34O's. The second-order (seagull) terms were chosen ultralocal so that they add no q dependence. An exactly displacement-invariant second-order lattice action could change the q² parts. 34S–34Z found no exact cubic gauge identity on Y1, so no such action is known to exist. The S2 numbers are therefore for the stated model, and the Ward violation says that this model, like Y1's cutoff, breaks spatial diffeomorphism invariance at O(1).
- The dynamic part of S2 (λ and the mode count) was not run.

**What S2 adds (inference).** The static, Newtonian-plus-bending sector fails more badly than the spin-2 sector. The induced theory on Y1 is not GR in any sector except the lapse's multiplier structure. That strengthens the 37A lesson: induced gravity needs a substrate whose cutoff respects the spacetime symmetries, not a lattice that only approximates them.

Code: [src/exgc37a_s2.py](../../src/exgc37a_s2.py) · results `s2_run.json`, `s2_check.json`.

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
