# EX-GC03 (EXPLORATORY): which coupling would make free fall universal?

Run by Claude, 2026-11-18, under the PI's go-ahead ("I think the coupling test should be next"). About 0.15 core-hours.
Radial 3D, linear dipole theory, 3 couplings × 7 particle types. Exploratory: this supports no claim on its own.

## Result
- **Active and passive gravitational mass agree.** For both working couplings (consumed G → X, consumed G → Y), the
  push a gradient exerts on a particle, p, is proportional to the particle's own G consumption M within about 5%
  across all seven types (p/M = 0.92–1.07 of the mean). So a particle's "charge" for receiving gravity matches its
  "charge" for making it. That is the Gc analogue of Newton's third law, and it holds without tuning.
- **Inertia does not follow M.** The drag is dominated by the X/Y (Brusselator) pattern that every type shares: the X
  weight is 1.11–1.15 of the total for all types. So the drag is nearly the same for all of them, and the fall speed
  μ = p/γ simply grows with M. The spread across types is 35% (the k = 0.1 particle falls at half the mean speed, the
  k = 0.4 particle at 1.7 times it).
- **"G feeds the core species" gives no gravity at all.** When consumed G becomes W, the particle still makes a 1/r
  well (M = 4–12), but its mobility is exactly zero. The reason is that X and Y no longer feel W or G, so the pattern
  that sets the particle's position is blind to the field. This is a structural zero: a source that does not fall.
- **Slowing G helps only when tuned.** Slowing G adds a drag from the particle's own G cloud ("field drag") that grows
  with M, faster than M itself (the G weight rises about 8-fold for M from 3.3 to 7.8). At τ_G ≈ 600–800 the spread
  falls to about 4–5%, but it rises again on either side. Optimising all four response times gives the same 4.4%
  floor (at the bounds τ_G = 30, τ_X ≈ 0.15). That is tuning, not structure.
- **No common low-drag (inertial) regime.** EX-GC02's low-drag regime would need the drag of every type to vanish at
  the same response times, which needs proportional weight vectors. The W weights differ by a factor of 3.5 across
  types, so it does not exist. (The optimisation in section C of `results.json` minimised a signed drag and drove
  every type into self-propulsion. It is not meaningful and is superseded by this argument.)

## Diagnosis and design implication
Gravity in Gc acts on M, but inertia comes from a structure (the X/Y core) that is common to every particle and does
not scale with M. Universality needs drag (or mass) ∝ M. Two routes:
1. **Field inertia:** inertia carried by the particle's own G field, which scales with M. The slow-G result shows the
   idea works in direction, but the G cloud's weight grows faster than M, so an exact law would need a field whose
   inertia is linear in its source (as electromagnetic field energy is quadratic, which is why that route failed
   historically). Not structural here.
2. **One elementary particle type (most promising).** If all matter is built from identical Gc particles, a composite
   of n particles has M = nM₁, push np₁ and drag nγ₁ at linear order, so its fall speed equals that of a single
   particle automatically. Violations then come only from binding effects, which are small if binding is weak. This
   is the testable next step: does a bound pair of identical Gc particles fall at the same speed as a single one, and
   by how much does binding change it?

## Records
`src/exgc03_coupling.py`, `results.json` (per-type weights, p, M, eigenvalues; summaries),
`paper/figures/EX-GC03_coupling.png`. The τ_G scans are in the figure script, computed from `results.json`.

## Limits
- Linear theory, radial single particles, seven types near one reference. Only types reachable by continuation.
- GW particles were stable (as their eigenvalues show) but carry no gravitational response.
