# EX-GC35F — asynchronous local-event order in a linear Y1 host

2026-11-24 | EXPLORATORY | ChatGPT | [Registered plan](../../spec/EX-GC35F_plan.md) · [Claude review](../../review/EX-GC35CD_claude_review.md). **PI CLOSED 2026-11-24.**

**Schedule test passed on the 5³ periodic toy.** All 24 seeded legal schedules at each of two timesteps produced bitwise-identical scalar and metric fields at every recorded depth. A stale-read control gave different results under different schedules. Common-slice constraints started at roundoff and developed a finite timestep defect.

## Event rule

The scalar event at site x and depth n+1 reads the fixed depth-n site and its six neighbors and the depth-(n−1) site. It writes one new scalar value by the centered Klein–Gordon difference. The spatial metric event reads its old and current versions, the depth-n neighborhood for the spatial Y1 stencil and the GC34O stress, and the one same-site scalar depth-(n+1) value needed for the centered scalar velocity. Each metric event writes its six component values at their inherited staggered placements. This is the one-way coupling used in the linear benchmark. Metric events do not change scalar inputs.

The staggered metric stencil was recovered by multiplying the Y1 symbol by the relative input and output placement phases and Fourier transforming it to integer-site offsets. The one-cell support follows the products of the single-step lattice derivatives. Outside the retained patch and in the imaginary part, weights were below 9e−16. Reapplying the recovered real stencil to 128 off-grid wavevectors gave a maximum symbol difference of 5.69e−14. The independent synchronous reference used the original Fourier acceleration and the GC34O stress code.

The scheduler only emits an event when all its versioned inputs are available. No barrier is used. For example, a scalar event at depth two ran after 114 of the 250 events at depth one. This is possible because the rule only depends on a local subset. A finite dependency graph with immutable inputs and single writers gives order independence by induction over its dependencies. The check here establishes that the prescribed graph reproduces this discretization and permits overlapping depths. The graph, integer depth labels, supplied initial slices, and one writer per field version are part of the construction.

## Results

| Δt | Events per schedule | Matching schedules | Schedules with depth overlap | Max metric difference vs synchronous | Metric difference between two stale schedules | Nonzero-mode constraint residual at t=.24 |
|---:|---:|---:|---:|---:|---:|---:|
| .04 | 1,750 | 24/24 | 24/24 | 1.39e−17 | 1.03e−3 | 6.60e−3 |
| .02 | 3,250 | 24/24 | 24/24 | 2.43e−17 | 4.83e−4 | 1.64e−3 |

Scalar values matched the synchronous reference bit for bit. The stale control removes the same-cell future scalar dependency; when that input is absent, it substitutes the current scalar value for the future one. Its different output confirms that order independence relies on versioned dependencies.

Initial nonzero-mode constraint residuals were below 5e−15 and fell by 4.01× at the final common slice when the timestep was halved. The spatial evolution equations agreed with the synchronous acceleration at roundoff. Those are separate checks from schedule equality. The toy evolves the homogeneous spatial metric locally and omits the homogeneous constraint from the comparison. Global compatible initial data were constructed using the Fourier solver. The test does not establish constraint preservation on irregular event cuts, equivalence of arbitrary competing graph rewrites, or a microscopic origin of the shared geometry.

Total CPU for the final runs was under 10 seconds. [Source](../../src/exgc35f.py) · [Results and hashes](results.json)
