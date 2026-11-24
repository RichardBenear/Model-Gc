# EX-GC34R — leading gravitational feedback on a scalar packet (EXPLORATORY)

2026-11-23 (America/Denver) | ChatGPT | PI go-ahead: “go ahead with GC34R”. [Plan](../../spec/EX-GC34R_plan.md) and registration entry pushed before computation (67d4a89 and cb38836). PI closure pending.

## Registered result

**PASS on R0–R2.** The leading correction to the scalar from its own zeroth-order weak field agrees with the corresponding continuum calculation. The registered observable is the first variation of the proper-volume integral of the squared field, evaluated on the specified synchronous slices. The two contributions are twice the integral of the flat field times its scalar correction and the integral of the half metric trace times the squared flat field.

| Time | Coarse response | Fine response | Continuum response | Fine error | Refinement ratio |
|---|---:|---:|---:|---:|---:|
| 1 | 29.46418 | 29.36886 | 29.33754 | 0.107% | 4.043 |
| 2 | −0.459027 | −0.444332 | −0.438771 | 1.267% | 3.643 |

- **R0:** The initial profiles agree within 1.4e-16 at common sites. The scalar variation matches the GC34O action's numerical variation within 8.3e-11. The largest gravity constraint residual was 1.38e-8.
- **R1:** Both fine responses meet the 5% threshold and improve by more than 3× under refinement. The fine scalar-motion terms alone differ from the continuum by 0.21% and 1.10%. At t=2 the scalar and volume terms partly cancel, which makes the total's relative error more sensitive.
- **R2:** Halving the timestep changed the response by at most 8.7e-6 of the coarse-to-continuum discrepancy. Centered differences at metric strength ±5e-5 match the tangent response to 1.4e-7 relative or better. The maximum metric-component magnitude at the larger finite-difference strength was 9.31e-5.

## Method and limits

The run uses the GC34Q standing preparation with the smaller, band-limited profile registered for this test. The lattice variation retains both the GC34O normal-stress product and the plaquette shear terms. The continuum calculation uses the spectral spatial divergence of the linearized Klein–Gordon action. Both use matched synchronous constraint initialization, and both omit the periodic homogeneous gravity mode.

The scalar correction is evolved on the metric sourced by the original flat scalar. The finite-strength check uses the same fixed metric history and evolves the scalar with the linear-in-metric action. It tests the tangent implementation without updating the gravity source. A reciprocal evolution with the corrected scalar resourcing its metric, including the required next-order conservation terms, remains open. This result applies to one generic standing packet at the two registered times and the specified time slicing.

No registered criterion changed. Run `python src/exgc34r.py`. The run used Python 3.13.2 and NumPy 2.2.4 on Windows and 114.5 CPU seconds (0.032 core-hours, cap 2).

[Full numerical record](results.json) · [Runner](../../src/exgc34r.py)
