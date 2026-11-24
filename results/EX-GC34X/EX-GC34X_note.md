# EX-GC34X — short-range cubic symmetry coefficient solve (EXPLORATORY)

2026-11-24 | ChatGPT | PI go-ahead: “go ahead with the proposed discrete symmetry idea experiment.” [Plan](../../spec/EX-GC34X_plan.md) and log registered before execution (7417256 and 0b2cb0d). PI closure pending.

**X0 PASS, X1 FAIL.** No coefficient vector in the tested 72-parameter family satisfies even the registered on-shell sample to tolerance. The failure provides a lower bound for this family, beyond reporting the residual of one fitted action. X2 and X3 were retained as diagnostics. No mixed-matter or reciprocal evolution was run.

## Family and feasibility test

Independently translate the first two slots of each of GC34W's three ordered cubic contractions by up to one axis step relative to the third slot. Pair translations with simultaneous sign reversal. There are 24 nonzero offset pairs per contraction. Each coefficient multiplies a form factor of the form cos(q₁·u + q₂·v) − 1, preserving the leading continuum vertex for any fixed finite coefficient vector. Responses are varied in the ten independent metric components.

The real linear system has 1,344 rows and 72 columns. It contains the repeated and nonparallel GC34W controls plus 40 seeded nonparallel mode pairs, with both polarizations. The eight GC34W control contractions reconstruct within 2.3e-16, and the largest free input residual is below 3.5e-16. The normalization uses the unshifted source norm and each gauge-vector norm.

| Quantity | Result |
|---|---:|
| SVD rank at relative cutoffs 1e-10, 1e-12 and 1e-13 | 72 at all three |
| Least-squares coefficient norm | 24.704 |
| Maximum sampled contraction, unshifted | 0.15243 |
| Maximum sampled contraction, least-squares fit | 0.01542 |
| Least-squares residual L2 norm | 0.07642 |
| Lower bound on the maximum complex gauge contraction | 0.002948 |
| Difference between the solve residual and left-space projection | 1.26e-14 |
| Registered maximum contraction tolerance | 1e-10 |

The lower bound is the norm of the component of the target outside the sampled column space divided by the square root of the number of complex conditions. It is a numerical infeasibility check for the specified family and sampled constraints, not a bound on other cubic actions.

## Holdout and continuum behavior

With the same coefficients, the normalized holdout maxima were 1.41 at h=.5 and 0.0214 at h=.125, relative to the unshifted source norm times the gauge-vector norm. X1 already failed, so these are diagnostics. On GC34W's fixed nonparallel controls the worst relative full-source difference from the continuum was 12.05%, 3.08%, and 0.773% at h=.25, .125, and .0625. That supports preservation of the leading continuum interaction while showing the finite-spacing symmetry problem remains.

## What is left open

Nearest-neighbor translations of the three vertex-averaged GC34W contractions are insufficient. A distinct candidate would change where tensor contractions are made or introduce other interaction terms, and test the complete cubic gauge identity including the first nonlinear gauge transformation. Fitting more coefficients of the same form without changing the family cannot repair these sampled equations.

Run `python src/exgc34x.py`. Python 3.13.2, NumPy 2.2.4, Windows, 1.56 CPU seconds for the recorded run. The first launch reached output serialization and failed on NumPy integer metadata. Converting those integers to Python integers fixed output only. The seeded sample, equations, thresholds and scoring were unchanged, and both runs took under five CPU seconds. The reconstruction uses GC34W's action factors. The residual is independently compared with the left-space projection of the matrix assembled from those factors.

[Numerical record](results.json) · [Runner](../../src/exgc34x.py)
