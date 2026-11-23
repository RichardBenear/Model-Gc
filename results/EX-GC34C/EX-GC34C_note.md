# EX-GC34C — what selects the gravity field's own equations, and hence light bending? (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "go ahead with 34C". Plan `spec/EX-GC34C_plan.md` (52b0c96,
pushed before the analysis). Code `src/exgc34c.py` (symbolic) → `analysis.json`; `src/exgc34c_check.py` (numerical)
→ `check.json`. Static, linear, weak field. Our construction, not LaViolette's Model G.

## Conclusion

**Two rules together select GR's light bending (γ = 1), and neither does it alone.**

1. **"Pure clock":** the lapse is a clock rate for gravity too. Gravity's own energy is multiplied by N, like every
   matter field's, and no energy is stored in gradients of N.
2. **Covariance:** gravity's static energy is one relabelling-invariant functional of the spatial metric. With two
   derivatives the only such functional is √γ R (the spatial curvature).

With both rules the static field equations are exactly those of linearised GR: γ = 1, Newtonian lapse, light bent
twice the lapse-only amount. Drop either rule and γ becomes a free number: 0 for independent fields, −1 for
Nordström, anything in between or beyond for the general pure-clock family.

**Rule 1 is not an extra assumption here.** It is the shared-structure rule from 34A/34B (all energy is static × N
plus kinetic / N) applied to the gravity field itself, which is what "gravitational energy gravitates like any other
energy" asks for. **Rule 2 is the open question for a lattice.** The link weights are not automatically relabelling
invariant. A known lattice construction that has this property in the continuum is Regge calculus (curvature from
link lengths around each hinge). That would be the natural next step.

## Results (registered predictions)

| action | γ (symbolic) | deflection ÷ lapse-only (numerical, b = 8, 12, 16) | prediction |
|---|---:|---:|---|
| A0 independent fields | 0 | 1.000 | P1 ✓ |
| A1 Nordström | −1 | 0.000 | P2 ✓ |
| A3 pure clock + covariance (static Fierz–Pauli = ADM N√γR) | **1** | 2.000 | P3 ✓; lapse identical to A0's Newtonian one at the standard normalisation |
| A2 pure clock, general coefficients | (a c₃ + 2a c₄ − 2b c₁ − 2b c₂ − b c₃) / [4(2c₁² + 2c₁c₂ + 2c₁c₃ + 6c₁c₄ + 4c₂c₄ − c₃²)] | example point: 0.906 = 1 + γ | P4 ✓ (not fixed) |
| A2 with linear gauge invariance (a = −b, c₂ = c₄ − c₁, c₃ = −2c₄) | −b / [4(c₁ + 2c₄)] | example: 0.750 = 1 + γ | still not fixed |

- **The lapse is a pure multiplier in GR's static action.** In Fierz–Pauli the (∂n)² terms cancel exactly
  (coefficient 0 in `analysis.json`), so the pure-clock rule holds in GR. Varying n then gives a constraint that ties
  spatial curvature to energy density. That constraint is the mechanism that makes the spatial metric respond to mass
  even though a body at rest has zero total stress (34B).
- **Linear gauge invariance is not enough.** It removes the dependence on the longitudinal part of h, but the ratio
  between the n–h coupling (from the linear part of the energy) and the h–h stiffness (from its quadratic part) stays
  free. Only a single covariant functional fixes that ratio, because both parts come from expanding the same √γR.
- **Numerical check.** Fields were built on a 128³ grid by inverse FFT, and straight-ray deflection was integrated
  directly. The ratios equal 1 + γ to machine precision at every impact parameter, including the non-isotropic A2
  point. That confirms the longitudinal part B drops out of light bending, so γ = −A/(2n) is the right gauge-free
  readout. This check shares the solved spectra with the symbolic path, so it verifies the readout rather than the
  algebra.

## Limits

Static and linear only. The structure has no shift vector, so moving sources (gravitomagnetism) and gravitational
waves are outside this analysis. Their consistency is a separate question: the full Fierz–Pauli structure needs a
shift-like variable. The 2-derivative restriction excludes higher-curvature terms, which would matter only at short
range.

## Compute

Seconds (symbolic 0.3 s; FFT check < 1 min). Linux cloud, Python 3.13.18, sympy 1.14.0, numpy 2.5.3.

## Suggested next (not approved)

**34D:** give the link weights a Regge-type gravity energy (deficit angles around lattice hinges, multiplied by the
lapse). Check numerically that the static field of a point source and of the knot gives γ ≈ 1 on the lattice, and
measure the bending with the 34A-r machinery. Then ask whether a shift-like variable is forced for moving sources.
