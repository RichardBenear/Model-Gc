# EX-GC38A — is gravity's stiffness the vacuum's entanglement? (area law vs induced G on the Y1 lattice)

2026-11-26 | EXPLORATORY | Claude | **DRAFT: awaiting PI go-ahead. Nothing has been run.** Default cap 2 core-hours; expected under 0.5. This is our construction, not LaViolette's Model G.

## Why

Richard suspects (2026-11-26) that matter gravitates through whatever mechanism produces entanglement. A precise version exists in the literature:
- Sakharov's induced 1/G (computed in 37A) equals the vacuum entanglement entropy per unit area of the same fields, times 4 (Susskind–Uglum 1994; Jacobson 1994).
- Jacobson (1995, 2016) then derives Einstein's equation from entanglement equilibrium.

If the identity S/A = 1/(4G_ind) holds on our lattice, gravity's stiffness *is* the vacuum's entanglement.

## Method

The free lattice scalar is the 37A one: Y1 dispersion, mass m, in its Gaussian ground state.

- **Entanglement across a plane.** The problem separates by transverse momentum into 1D chains. Each chain's half-space entropy comes from the correlation matrices X = ⟨φφ⟩ and P = ⟨ππ⟩ restricted to the region, using symplectic eigenvalues ν of √(XP):

  S = Σ [(ν + ½) ln(ν + ½) − (ν − ½) ln(ν − ½)].

  Use blocks of L sites with two boundaries and divide by 2. Converge in L.
- **Orientations:**
  - plane normal to z (lattice face): on-site term 2 + m² + K_x² + K_y², hopping 1;
  - plane normal to (1,1,0): on-site 4 + K_z² + m², hopping 2 cos(κ/2), area √2 per site.
- **Masses:** m = 0.05, 0.1, 0.2, 0.3, 0.5, 1.0.
- **Induced 1/(4G):** 37A's TT kinetic coefficient, 1/(4G) = 4πA, for both polarisations, at the same masses. Grids are refined until stable at small m.

## Readouts and predictions

The cutoff-dependent parts of both quantities are regulator-specific. The **m² ln m terms are universal**, independent of the cutoff: in the continuum both carry the same coefficient. Expected value for one real scalar: −m² ln(m²)/(24π) per unit area, up to the sign convention of the fit.

- **P1 (scored):** the fitted m² ln m coefficient of S/A (m ≤ 0.3) matches that of 4πA (cross polarisation) within 10%. It also matches the continuum value within 10%, for both plane orientations.
- **P2 (descriptive):** the full ratio (S/A) / (1/4G_ind) for each orientation and polarisation, which tests whether the cutoff-dominated parts also match on this lattice. Expectation: not equal (Y1's anisotropy, 37A), but the same order.
- **P3 (descriptive):** orientation dependence of S/A, compared with the cross/plus anisotropy of 37A.

## Checks

- **G1:** 1D chain entropy against the Calabrese–Cardy massive form, (1/6) ln(ξ) per boundary with ξ = 1/M, at small M.
- **G2, separate path:** brute-force 3D correlation-matrix entropy on a 12³ periodic box cut by two planes. It must match the separated result for the same box.

## Outputs

`src/exgc38a.py`, `results/EX-GC38A/`.
