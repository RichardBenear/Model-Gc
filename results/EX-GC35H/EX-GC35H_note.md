# EX-GC35H — native staggered Skyrme stress

2026-11-24 | EXPLORATORY | Design, implementation and runs ChatGPT. [Plan](../../spec/EX-GC35H_plan.md) · [Claude review](../../review/EX-GC35EF_claude_review.md). **PI CLOSED 2026-11-24.**

**Full static Ward gate failed.** The new coupling reduces the full-source incompatibility on both tightly relaxed h=.2 knots. It leaves a finite residual and gives an explicit obstruction to extending the free-scalar identity by simply adding components constrained to a unit vector.

| Knot | Degree integral | Projected force density | Inherited residual | Staggered candidate residual | Energy-only residual |
|---|---:|---:|---:|---:|---:|
| B1 | −0.93012 | 9.60e−6 | 1.654% | **0.962%** | 4.36e−16 |
| B2 | −1.87864 | 9.31e−6 | 1.310% | **0.679%** | 4.16e−16 |

The flat energy and the knot relaxation are the GC35D ones. The candidate changes how the action varies with the metric. The sigma-model term uses the forward/backward product for diagonal stress and face-averaged gradients for shear, as in 34O. The quartic term averages pair areas from face corners for diagonal stress and evaluates the shear Gram minors on the face with both transverse gradient orientations. Flat energies matched GC35D to roundoff. Finite differences of the sigma and quartic metric variations agreed to 6.7e−11 and 5.0e−9 relative. The static Y1 residual uses the same nonzero mode cut as GC35D. No source divergence projection was applied.

Coupling just the sigma stress or just the quartic stress with the same energy density left residuals around 10–11%. Combining them reduced the residual below 1% for both knots. Thus the candidate retains substantial cancellation between the stresses, with an unresolved finite-spacing remainder.

## Why the unconstrained sigma identity does not transfer directly

Write uₓ for the four-component unit field and Gₓ,ₐ=(uᵧ−uₓ)/h on the link from x to y=x+eₐ. The normalized staggered sigma stress satisfies

```
(D·T_sigma)ₓ,ₐ = (Δuₓ+Δuᵧ)·Gₓ,ₐ / (12π²).
```

The implementation checked this on random 3D fields to 5.4e−16 relative. In a constrained static sigma solution, Δu=λu. The remaining link term is (λᵧ−λₓ)(1−uₓ·uᵧ)/(h·12π²). The finite-difference chord has radial components at its endpoints, so tangent force balance at each site does not eliminate this reaction term.

An exact periodic rank-one control uses a sequence of collinear unit vectors with one sign flip and is constant along the other axes. Its tangent force and quartic stress vanish identically. The proposed sigma stress has a nonzero divergence with maximum 4.22 in the implemented normalization, exactly equal to the reaction formula. This is a discrete counterexample to a universal identity for this candidate. It gives no general obstruction to a different staggered constrained action.

Both knots meet the relaxation tolerance. The full-source residuals miss the 1e−10 compatibility gate, so the conditional moving-knot runs were not started. The two runs used 136 CPU seconds. Reproduce with `python src/exgc35h.py 1` and `python src/exgc35h.py 2`.

[Source](../../src/exgc35h.py) · [B1 data](b1.json) · [B2 data](b2.json)
