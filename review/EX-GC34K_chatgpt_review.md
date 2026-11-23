# ChatGPT review of EX-GC34K

2026-11-23 | Independent review of the exploratory results

**Recommendation.** Keep J1 as the linear vacuum-wave benchmark. Before changing the kinetic action or refining the grid, test a moving source that satisfies J1’s exact lattice gauge identities. The current source violates the re-slicing identity at first order in grid spacing, the same order as the reported curvature error.

## What GC34K established

- **Registered vacuum evolution passed.** With constraint-compatible initial data, J1 ran for 40 wave periods in both tested directions. Constraint residuals stayed below 7.8 × 10⁻¹². Permitted zero-frequency seeds changed the tidal signal by less than 7 × 10⁻¹³ relative to the wave. Together with the spectrum scan, this supports removal of the original growing vacuum mode in J1 for the tested conditions. Two low-frequency waves remain, with the extra edge mode at the lattice cutoff.
- **The moving-source benchmark failed.** The worst curvature error was 1.25% at kh = 0.05 against a 1% limit. It fell to 0.75% at kh = 0.03. All 33 sampled kernels showed first-order convergence, but the leading coefficient of J1’s worst curvature component varied substantially across the sample. This corrects GC34J’s claim that the identity fixes the coefficient. The sample does not establish that every allowed kernel must have a first-order error.
- **The centered-gradient obstruction is structural for the tested formulation.** At nonzero wavevectors whose components are 0 or π, the centered gradient vanishes and the lapse–edge coupling does not. No bounded kinetic kernel can satisfy the identity there, regardless of stencil width. The same argument applies to the staggered shift averaged back to vertices in the plan. A genuinely staggered formulation with a different coupling remains open.

## Moving-source compatibility

For the Hermitian lattice equation `Mx = f`, an exact gauge null direction `z` requires `z†f = 0`. A least-squares solver can return a field even when this condition fails.

In `I0g`, the re-slicing vector has lapse component `−iω` and shift components `−g`, where `g` is the forward lattice gradient and μ = −1. In the baseline `e4_detail` calculation, the forcing is proportional to `(1, −v, 0)`, with continuum frequency `ω = k·v`. The solve fixes the spatial displacement gauges and uses least squares. The benchmark does not report the re-slicing overlap or the residual of the full forced equations.

For the axis-parallel case that sets the reported error, the overlap is proportional to `iω + g*·v`. With `q = kh`, dividing its magnitude by the time-derivative contribution gives

`|cos(q) − 1 + i(q − sin(q))| / |q| = |q|/2 + O(|q|³)`.

At kh = 0.05 this is about 2.5%. This is a source-conservation mismatch, measured relative to the time-derivative term. Its contribution to the separate 1.25% curvature error has yet to be measured.

## Recommended next study

Hold J1 and the GC34D static kernel fixed.

1. Audit the existing source. For all nine cases, report its overlap with every gauge direction, the residual of the original forced equations, and the individual curvature components.
2. Derive moving-source coupling from a discrete continuity equation consistent with J1. Verify the source identities through the order being tested. Use a numerical projection as a diagnostic if helpful, with the derived source as the primary test.
3. Repeat the four-spacing curvature comparison with the corrected source. Use a curvature readout invariant under the lattice gauge transformations, or verify its invariance to the tested order. If the linear error disappears, the vacuum action can stay in place. If the error remains, test whether it can be eliminated within the allowed kernel space before changing the staggered geometry or potential.

Refining the grid reduces the current benchmark error. The source audit can distinguish errors in the geometry, matter coupling, and curvature readout. This is directly relevant to Model Gc’s larger question of whether moving matter can source the same geometry that its waves experience.

## Scope and records

I reviewed the registered plan, results, and relevant source, and independently evaluated the source overlap. I did not rerun the forced-field solves.

- [GC34K registered plan](../spec/EX-GC34K_plan.md)
- [GC34K results](../results/EX-GC34K/EX-GC34K_note.md)
- [GC34K source](../src/exgc34k.py)
- [Zone-boundary diagnostic](../results/EX-GC34K/posthoc_zb.json)
