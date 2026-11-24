# EX-GC34Z — joint cubic action and coordinate-transformation solve (EXPLORATORY)

2026-11-24 | ChatGPT | PI go-ahead: “try the next idea.” [Plan](../../spec/EX-GC34Z_plan.md) and log registered before execution (2fe720b and 499c2d4). PI CLOSED 2026-11-24 (Richard: “I sign off on all of the open items”).

**Registered verdict: Z0 PASS, Z1 PASS, Z2 FAIL.** The joint fit satisfies the training identities to roundoff. It fails on held-out triads, including at the training spacing. A separate expanded-sample solve gives a numerical lower bound that rules out satisfying those same-spacing equations within this 120-parameter family. No scalar coupling or time evolution was attempted.

## Registered test

Y1's quadratic operator and GC34Y's majority-parity placement are fixed. The fit varies 72 cubic action coefficients and 48 coefficients of the first nonlinear coordinate transformation. Each correction uses an inversion-even translation difference, which vanishes in the leading continuum coefficient for fixed amplitudes.

- **Z0:** The construction reproduces the GC34Y baseline within 9.30e-16 normalized. The continuum identity closes within 1.30e-14.
- **Z1:** The combined training system has 128 complex conditions. Its rank is 108 at each registered cutoff. A least-squares fit with coefficient norm 430.44 lowers the maximum normalized residual from 0.1665 to 2.70e-13. The eight free-wave channels alone also admit a fit on this sample. Its action-only subsystem has rank 12 and residual at most 1.24e-16.
- **Z2:** The same joint coefficients give held-out maxima of 58.66 at h=.5 and 1.832 at h=.125, against a registered limit of 1e-10. The same-spacing diagnostic at h=.25 gives 6.988. Refining those fixed triads to h=.0625 reduces the maximum to .4633. The correction terms vanish under further refinement, but a continuum limit does not imply an exact identity at finite spacing.

## Post hoc diagnostic: is this only a choice of fit?

Add the 12 held-out triads to the 24 training triads, and retain the eight free-wave channels. Solve again within the same family, without changing the registered score.

| Expanded system | Rank | Best-fit maximum | Lower bound on any maximum |
|---|---:|---:|---:|
| All triads at h=.25 | 120 | 0.00720 | 0.00179 |
| All triads at h=.5, .25 and .125 | 120 | 0.12067 | 0.01372 |

The lower bound uses the component of the target outside the column space and divides its L2 norm by the square root of the number of complex equations. Ranks are unchanged across the three SVD cutoffs, and the independent residual and left-space projection agree within 4.5e-14 in the expanded solves. The same-spacing result rules out a coefficient choice satisfying this expanded sample within the tested family. The off-shell failure is therefore more than a poor choice among the fits to the original training sample.

## Implication and scope

The tested translation corrections leave the tensor contractions and the basic form of the coordinate transformation fixed. A different local action, another transformation family, or a formulation with different fundamental variables remains possible. This is a sampled numerical incompatibility result for a particular family. It is not a general no-go theorem for exact lattice gravity. A next design should introduce a structural identity or new independent interaction terms before more fitting. GC34U remains the isolated reciprocal continuum benchmark and linear Y1 remains unchanged.

Run `python src/exgc34z.py` for the registered test and `python src/exgc34z.py posthoc` for the expanded diagnostic. The registered run took 7.38 CPU seconds on Windows, Python 3.13.2, NumPy 2.2.4. Both runs were within the 2 core-hour cap. The baseline check reuses GC34Y tensor formulas with a separate assembly of the candidate coefficients, and the SVD residual is checked against its left-space projection. No registered criterion changed.

[Registered data](results.json) · [Post hoc data](posthoc.json) · [Runner](../../src/exgc34z.py)
