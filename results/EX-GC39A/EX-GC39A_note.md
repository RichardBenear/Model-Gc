# EX-GC39A — entanglement equilibrium on the lattice (does matter set the geometry?)

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC39A_plan.md) · [design note](../../spec/EX-GC39_design.md). PI closure pending. Cloud container, Python 3.13, numpy, python-flint. About 0.1 core-h. This is our construction, not LaViolette's Model G.

**Not as Jacobson's argument needs.** On our lattice, the entanglement of a block does respond both to curvature and to matter, and the first law holds exactly. Imposing "entanglement equilibrium" does give a Newton-like G that is roughly the same for every block size. But:
- **Both responses grow like the block's surface (ℓ²), not like Jacobson's ℓ⁴.** They are cutoff effects at the block's boundary, not the bulk curvature and energy balance his derivation uses.
- **Even the sign of G depends on how "fixed volume" is defined on a lattice.** The registered rule gives G_eq ≈ +230 to +300 (attractive). A scale-based rule gives about −9 (repulsive).
- **Neither value matches** 1/(4η) ≈ 11–13 from 38A, or the induced G ≈ 21–29 from 37A.
- **Classical (coherent) matter has no entanglement signature at all,** so it could not gravitate through δS.

The conclusion I draw: on a lattice whose geometry is carried by weights on fixed links, entanglement equilibrium does not single out a geometry. For "matter sets geometry through entanglement" to work, the counting of sites and links would itself have to be the geometry.

## Results (periodic 16³, cube blocks of side ℓ)

| m | ℓ | δS (curvature, per ε) | η_eff (registered) | η_eff (scale-corrected) | k_B (matter) | G_eq (registered) | G_eq (scale-corrected) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.3 | 2 | 1.54e-4 | −0.044 | 1.15 | −657 | 227 | −8.6 |
| 0.3 | 3 | 3.40e-4 | −0.021 | 0.53 | −1540 | 231 | −9.1 |
| 0.3 | 4 | 5.90e-4 | −0.012 | 0.30 | −2767 | 225 | −9.3 |
| 0.3 | 5 | 8.31e-4 | −0.008 | 0.19 | −4270 | 227 | −9.8 |
| 0.5 | 2 | 5.19e-4 | −0.032 | 1.11 | −635 | 299 | −8.6 |
| 0.5 | 3 | 1.16e-3 | −0.016 | 0.51 | −1510 | 291 | −9.3 |
| 0.5 | 4 | 2.02e-3 | −0.010 | 0.29 | −2735 | 267 | −9.5 |
| 0.5 | 5 | 2.85e-3 | −0.007 | 0.18 | −4239 | 250 | −10.1 |

38A's area coefficient: η = 0.0220 (m = 0.3) and 0.0196 (m = 0.5), so 1/(4η) = 11.4 and 12.7.

Definitions:
- **η_eff** is the change in block entropy at fixed volume, per unit change of the block's fixed-volume area under a curvature profile.
- **"Registered"** subtracts the block's response to a uniform ψ equal to the block's mean ψ, as the plan specified.
- **"Scale-corrected"** subtracts ⟨ψ⟩·ℓ dS/dℓ instead: the entropy change from shrinking the block back to its original proper volume. This alternative was computed alongside, and it is descriptive.
- **k_B** is the first-order entropy change per unit local energy density, for a smooth, long-wavelength squeezing of the vacuum.

**Scoring:**
- **P1 (η_eff > 0 and k_B > 0): FAIL.** Both are negative under the registered rule. Their ratio still gives G_eq > 0, i.e. attraction. Under the scale-corrected rule η_eff > 0 but k_B < 0, so G_eq < 0.
- **P2 (η_eff/η = 1 within 10%): FAIL, as expected.** The ratio runs from −2.0 to −0.37 (registered) or from 57 to 8.5 (scale-corrected), and depends strongly on ℓ.
- **P3 (G_eq independent of ℓ within 20%): PASS.** Registered: 225–231 at m = 0.3 (2.5%) and 250–299 at m = 0.5 (18%). Scale-corrected: −8.6 to −10.1 (15%).
  - The reason it passes is not Jacobson's. The curvature response, η_eff × (area change), and the matter response k_B both scale as about ℓ², so their ratio is roughly constant.
  - The area change itself scales as ℓ⁴, as in the continuum. So the entropy response is a boundary (cutoff-scale) effect, not an area-deficit effect.
- **P4 (descriptive):** the registered G_eq is about 20× 1/(4η) and about 8–10× 37A's G_ind. It also differs between m = 0.3 and 0.5, so it is not one universal constant. The scale-corrected |G_eq| ≈ 9–10 is close in size to 1/(4η) ≈ 11–13 but has the wrong sign. I don't read anything into that.
- **P5 (coherent matter has δS = 0): holds, structurally.** A Gaussian state's entropy depends only on its covariances, and a displacement does not change them. So a classical field configuration, which certainly carries energy, leaves every block's entanglement unchanged.

## Checks

- **G1: PASS.** The flat ground state from dense diagonalisation matches the FFT correlators to 9e-16. The slab check against 38A's separated method was not repeated; 38A's G2 already validated the FFT box.
- **G2, linearity:**
  - Matter: PASS, 0.03–0.11%.
  - Curvature at m = 0.5: PASS, 0.05–0.14%.
  - Curvature at m = 0.3: 1.6–2.2%, over the 1% gate, so FAIL. The central difference used for the results cancels the quadratic part.
- **G3, first law, separate path: PASS.** For ℓ = 2 and 3, δ⟨K_B⟩ from an arbitrary-precision entanglement Hamiltonian equals the finite-difference δS_B to 2e-7 and 5e-7.

## What it means (inference)

- **H-A ("matter sets the geometry through entanglement") does not come out cleanly on our lattice.** The ingredients are there: entanglement responds to curvature, and the first law holds. But Jacobson's chain needs two things:
  - the area-law term to measure *proper* area (a cutoff fixed in proper units);
  - the matter term to be a bulk energy integral.

  On a lattice with fixed sites, both responses come from the cut boundary, and the answer depends on how "same volume" is defined. That freedom is exactly what a lattice with weights on fixed links does not fix.
- **This sharpens 37A's lesson.** For geometry to be *derived* from entanglement, the substrate's elementary counts (sites, links, events) have to be the geometry itself. Area would then be "number of cut links" and volume "number of sites", as in causal sets or graph models, not metric weights painted on a fixed grid.
- **The coherent-state result is a general warning.** Classical field energy has no entanglement, so any "gravity from entanglement" mechanism must act through δ⟨K⟩, i.e. energy, not through δS. Then it is energy that gravitates, and entanglement only supplies the coefficient.

## Files

[Code](../../src/exgc39a.py) · `m0.3.json`, `m0.5.json`
