# EX-GC49 — firing-rate changes carried by the same coincidences as light

2026-11-27 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC49_plan.md), with amendment A1 made before any run | About 0.05 core-h on the cloud container (Linux, Python 3.13, numpy 2.x). This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**The literal idea 1 gives the wrong second-order gravity.** Suppose each site's firing rate N is passed to its neighbours by the same two-site coincidences (rate N_xN_y) that carry light. Then the rule derives a field equation for N: N³ is harmonic outside a mass. The field shares the light cone. But the metric's second-order term comes out as PPN β = −1/2, where GR has β = 1.

With 37B's γ = 1, that makes the perihelion advance 1.5 times the GR value, about 64″ per century for Mercury instead of 43″. That is ruled out.

**Only one of the four registered transported quantities gives β = 1: the coordinate slowness 1/N².** That is the refractive index light sees in 37B, where the coordinate light speed is N². If a coincidence passes on 1/N², then ln N is harmonic and the exterior metric is the exponential one, N = e^{−U}. That matches GR at this order.

The lattice reproduces the continuum prediction for each quantity once a lattice artefact that is linear in the source strength is separated out (post hoc, see below). The registered lattice fit, which did not separate that artefact, missed the ±0.05 window for three of the four quantities.

**So the mechanism question becomes sharper.** The coincidence rule can supply both the wave equation for N and the shared light cone. What it does not decide is which quantity gets passed on. Matching GR requires it to be the slowness of light, not the firing rate.

## Scorecard

| Prediction | Result | Score |
|---|---|---|
| P1 lattice β within ±0.05 of −½, 0, ½, 1 (joint fit) | −0.558, −0.055, 0.448, 0.951 | **1 of 4** (slowness passes; the others miss by 0.002–0.008 beyond the window) |
| P2 A = S/4π within 2% | 0.9995 at all 16 points | PASS |
| P3 continuum β within ±0.02, lattice vs continuum N within 1e-3 | Continuum β −0.5004, 0.0001, 0.5001, 1.0001; max lattice–continuum difference 1.25e-4 | PASS |
| P4 light vs δN Shapiro delays within 1% at k₀ = 1 | 0.05% | PASS |
| P4 arrival-difference log-slope in [−2.5, −1.5] | −2.63 | FAIL (falls faster than predicted) |

**Separate path** ([recompute.json](recompute.json)). The exact exterior form N = F⁻¹(−U), expanded with a second finite difference in U and no fitting, gives β = −0.5000, 0.0000, 0.5000 and 1.0000.

## Evidence

**49A, static mass** ([A.json](A.json))
- **Equation.** Σ_y N_xN_y [f(N_y) − f(N_x)] = s_x on a 129³ cube.
- **Source and boundary.** The source is a Gaussian core with σ = 2. The faces carry Dirichlet values from Gauss's law.
- **Strengths.** Q/10 = 0.015 to 0.06.
- **Solver.** The full nonlinear residual was brought below 1e-11, using a DST Poisson solve with Anderson mixing, in 8 to 18 iterations.

| f (what a coincidence passes on) | Harmonic quantity | Predicted β | Continuum fit | Lattice joint fit (scored) | Lattice, artefact separated (post hoc) |
|---|---|---|---|---|---|
| N (firing rate) | N³ | −0.5 | −0.5004 | −0.558 | −0.502 |
| ln N | N² | 0 | 0.0001 | −0.055 | −0.000 |
| −1/N | N | 0.5 | 0.5001 | 0.448 | 0.503 |
| −1/(2N²) (slowness) | ln N | 1 | 1.0001 | 0.951 | 1.006 |

**Why the scored lattice fit is offset** (post hoc diagnosis, not a re-score).
- All four quantities share the same offset, about −0.05. In each, the offset scales as 1/Q: for f = N it runs −0.70, −0.60, −0.57, −0.55 from the weakest to the strongest source. So it is a linear-order effect and has nothing to do with the nonlinearity.
- It comes from the cubic anisotropy of the lattice Green's function. The fit pooled axis and diagonal samples, which have different 1/r² corrections. Taken separately, the axis gives +0.011/Q and the diagonal −0.0074/Q. Pooled, they give −0.058/Q.
- Fitting the x² coefficient against Q (a₂ = c₁Q + c₂Q² + c₃Q³) separates out the linear piece. The Q² part then gives the continuum β to within 0.006, for the pooled, axis-only and diagonal-only samples alike.
- The scored result stays as recorded. The registered window was too tight for a pooled fit at these strengths.

**49B, one cone** ([B.json](B.json))
- **Setup.** A chain at impact parameter 5 across the f = N well at Q/10 = 0.045. Light follows the 37B rule. The firing-rate perturbation follows the linearised f = N rule, including its gradient terms. Both are evolved by RK4.

| k₀ | Shapiro delay, light | Shapiro delay, δN | Relative difference | Arrival difference |
|---|---|---|---|---|
| 0.25 | 7.837 | 8.067 | 2.9% | 0.230 |
| 0.5 | 8.441 | 8.457 | 0.19% | 0.016 |
| 1.0 | 11.329 | 11.335 | 0.05% | 0.006 |

- On a flat chain the two fields arrive at identical times.
- In the well they share the cone, and differ only at long wavelength, where the gradient terms act as a weak potential. The difference falls faster than k₀⁻² between 0.25 and 0.5, and close to k₀⁻¹·⁴ between 0.5 and 1. The fitted slope over all three is −2.63, so the registered window was the wrong shape.
- The one-dimensional chain keeps only the gradient terms along its own line. It tests whether the two cones coincide, not the full 3D linearised operator.

## What this does and does not show

**It shows:**
- A coincidence-carried firing-rate field is a candidate mechanism for the field equation, and it sets its own nonlinear completion.
- The literal form of that field (f = N) fails the perihelion test at order 1.
- The form that passes (slowness 1/N², exponential metric) is singled out by GR's β, not by anything in the counting rule. That is the gap to close if this line continues.

**It does not test:**
- the source side, i.e. what matter sources and whether that source is the energy (44/48);
- moving sources at second order;
- whether a rule that passes 1/N² is causally invariant or stable at large amplitude.

**Suggested next, not approved:** look for a reason inside the event picture for the transported quantity to be the light slowness. One candidate is that a coincidence passes on the delay it caused rather than a rate, since delays add along paths. Then check that rule's large-amplitude behaviour: the exponential metric has no horizon.
