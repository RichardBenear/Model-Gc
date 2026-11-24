# EX-GC34Y — majority-parity cubic contractions (EXPLORATORY)

2026-11-24 | ChatGPT | PI go-ahead: “go ahead and try the next idea.” [Plan](../../spec/EX-GC34Y_plan.md) and log registered before execution (5da145f and 50bdcc8). PI closure pending.

**Y0 PASS, Y1 FAIL, Y2 FAIL.** Moving products to staggered contraction sites improves the vacuum gauge residuals, but neither the on-shell nor the full cubic condition becomes exact. No matter coupling or reciprocal evolution was attempted.

## Construction

Keep the GC34N quadratic operator. For each primitive product of three factors in the cubic Einstein–Hilbert connection-product action, assign its contraction to the majority parity position component by component. For example, when two factors lie on links in a particular direction and the third lies on a vertex in that direction, average the vertex factor to the link. The unmodified all-vertex prescription is the control. Independently, implement the first nonlinear coordinate variation as the continuum Lie derivative with Y1's derivative symbol and average its two factors to the output component's placement. Test the sum of the linear gauge variation of the cubic action and the nonlinear gauge variation of the quadratic action.

The new assembly reproduces the GC34W vertex-only contractions to (2.85\times10^{-16}) in the registered control. With continuum derivatives both on-shell and full off-shell identities hold, with maxima of (1.07\times10^{-16}) and (1.45\times10^{-14}), respectively. The full off-shell control tests the relative sign and normalization of the cubic and quadratic variations.

## Registered results

| Spacing | Vertex on shell | Majority on shell | Vertex full cubic | Majority full cubic |
|---:|---:|---:|---:|---:|
| 0.5 | 0.04118 | 0.009991 | 9.602 | 0.8356 |
| 0.25 | 0.009958 | 0.002574 | 2.905 | 0.2188 |
| 0.125 | 0.002468 | 0.000648 | 0.7617 | 0.05532 |

Entries are maxima over the registered normalized contractions, with the on-shell and full tests using their own normalizations. Full cubic residuals use the sum of the absolute values of the three individual continuum terms in each case. The median full-cubic majority residual on the finest grid is about 0.006. Finite-spacing majority residuals drop approximately fourfold as the grid spacing halves. Both maxima remain above the registered (10^{-10}) gate.

## Interpretation

Placement is a meaningful lever: it lowers the worst on-shell residual by about fourfold and the worst full-cubic residual by about fourteenfold at the finest spacing. It leaves a nonzero on-shell defect, so changing the first nonlinear coordinate transformation alone cannot repair this particular cubic action. A different cubic action is necessary. If pursuing exact finite-spacing symmetry further, solve the action and transformation coefficients together from the algebraic identity, with the continuum coefficients fixed. This result does not rule out other action families and does not alter the linear Y1 or isolated continuum benchmarks.

Run `python src/exgc34y.py`. The verification reuses GC34W's averaging factors and GC34N's linear Einstein tensor. The vertex control independently compares the expanded component contractions against GC34W's tensor contractions, and the continuum control checks the full identity. Python 3.13.2 and NumPy 2.2.4 on Windows. 3.48 CPU seconds against the 2 core-hour cap. No registered criterion changed.

[Raw results](results.json) · [Runner](../../src/exgc34y.py)
