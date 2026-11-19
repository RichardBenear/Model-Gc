# EX-GC12 (EXPLORATORY): soliton route, Stage 1 – a stable lattice knot

Claude, 2026-11-19, at the PI's request (decision card: "Particle as a field soliton"). Design in
docs/SOLITON_PARTICLE_DESIGN.md. Our construction.

## Answer
**Yes: the lattice rules hold a degree-1 knotted texture, and the degree is exactly conserved.** Under the three-node
lift rule, a hedgehog texture shrinks within about 1000 sweeps to a compact knot of about 150 mismatched edges and
then stays there. Its degree stays at exactly 1 in every check (shown as −1 by our orientation convention). Under the
pairwise rule it unwinds to degree 0, as predicted. The particle has no Derrick-type collapse through the lattice,
because the lift rule never lets a neighbouring pair of lifts get more than 120° apart. The map stays "continuous",
so its degree cannot change.

## Runs (frames only, 25³, outer face fixed, 10 000 sweeps, three random test points for the degree)
| Start | Rule | T_f | Energy (from → late, min) | Degree throughout |
|---|---|---|---|---|
| hedgehog R = 10 | lift | 0.6 | 3528 → 164 (min 148) | 1, constant |
| same, seed 1 | lift | 0.6 | 3528 → 168 (min 162) | 1, constant |
| hedgehog | lift | 0.3 | 3528 → 132 | 1, constant |
| hedgehog | lift | 1.0 | 3528 → 3304 (disordered) | 1, constant |
| uniform (control) | lift | 0.6 | 0 → 30 (thermal) | 0 |
| hedgehog | pairwise (control) | 0.6 | 3528 → 36 | lost: 1 → 0 |
| two textures, +1 and +1, R = 5.5 | lift | 0.6 | 2040 → 436 | 2, constant |
| two textures, +1 and −1 | lift | 0.6 | 2040 → 566 | 0 (they do not annihilate within the run) |

- The degree uses three random points of S³, and all three always agree.
- Neighbouring lifts stay within 120° (minimum dot product 0.5) under the lift rule; the pairwise rule breaks this
  (minimum dot product −1).
- Thermal floor: about 30 mismatched edges in the uniform state at T_f = 0.6. The knot costs about 120–140 edges above
  that, so it is a particle of fixed size, a few lattice spacings across.

## What it means
- **The soliton route has a stable, conserved particle number.** The degree adds (two textures give 2), and a
  texture and an anti-texture carry 0. This is the object the Finkelstein–Rubinstein argument needs.
- **It is built from the same rule as route 3.** It uses the same lattice rotation system and the same three-node rule
  as EX-GC08 and EX-GC10, with no new field. The difference is the shape of the twist: wrapped around all directions
  (degree 1) rather than about one axis (degree 0 plus a fixed core).
- **The exchange sign is not computed yet.** That is Stage 2: build the 2π rotation and the exchange as explicit move
  sequences and compute their loop class.

## Method
`src/exgc12_soliton.py`:
- `Q`: the 48 lifts.
- `hedgehog`: snapped texture, with products for several centres.
- `degree`: Kuhn triangulation, signed covering count.
- `Field.sweep`: checkerboard Metropolis under `exgc08_belt.move_table`. Each move takes the lift a short step away
  from the old one.

Records: `*.json` (history every 500 sweeps). Platform: Linux cloud container, Python 3.13.18, numpy 2.5.3. Compute:
about 5 minutes.

## Limits
- Frames only, with no Gc chemistry, so the knot does not yet consume G (Stage 3 in the design note).
- The +1/−1 pair did not annihilate within 10 000 sweeps. Total degree 0 allows it, but there is a barrier.
- "Degree exactly conserved" is observed and argued, not proven. The argument is that steps under the lift rule are
  short, so each tetrahedron's image changes continuously.
