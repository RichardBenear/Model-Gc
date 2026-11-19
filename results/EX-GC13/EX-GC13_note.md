# EX-GC13 (EXPLORATORY): soliton route, Stage 2 – the exchange sign

Claude, 2026-11-19, at the PI's request ("GO AHEAD"). Design: docs/SOLITON_PARTICLE_DESIGN.md. The theorem being
checked is Finkelstein & Rubinstein (1968). The knots, the lattice rule and the test code are our construction.

## Answer
**Yes: for knot particles (degree-1 textures), swapping two of them is the same kind of loop as turning one by 360°.
Both are "knotted" (class 1), and the loops that should be trivial come out trivial. Each of these loops can be
carried out on our lattice entirely by allowed (no-cut) moves.** This is the classical form of the fermion signature:
a particle that needs 720° to return to itself also picks up the "minus" when two are swapped.

| Loop | Expected class | Computed (3 regular values each) |
|---|---|---|
| one knot turned by 360° | 1 | 1, 1, 1 |
| one knot turned by 720° | 0 | 0, 0, 0 |
| one knot carried around a circle (no turning) | 0 | 0, 0, 0 |
| a degree-2 knot turned by 360° | 0 | 0, 0, 0 |
| **two knots swapped** (each keeps its orientation) | **1** | **1, 1, 1** |
| two knots swapped, one also turned by 360° | 0 | 0, 0, 0 |
| two knots swapped twice | 0 | 0, 0, 0 |
| the pair turned rigidly by 360° | 0 | 0, 0, 0 |

The 24 computations agree with theory in every case, using different regular values. The first four rows check the
method; the last four are the test.

**On the lattice** (spacing 0.2, knot radius 10 lattice units, as in EX-GC12):
- The 360° turn, the swap, and the swap plus turn were sampled at 64 (and 128) times and snapped to the 48 lattice
  lifts.
- Every change, about 9 000 to 81 000 node changes per loop, could be made by an allowed move under the three-node
  rule. No change was stuck or too large.
- Neighbouring lifts stayed within 120°, and the degree stayed constant (1, or 2 for the pair).

So the lattice loops are discrete copies of the continuous ones. A straight-line interpolation between nearby lifts
never passes through zero, so the lattice loops carry the same classes.

## What it means
- **Spin and statistics now come together.** With the particle as a knot of the graph's own rotation system:
  - a 360° turn is a non-trivial loop and a 720° turn is trivial (spin-½-like);
  - swapping two knots is the same non-trivial loop (fermion-like);
  - degree-2 objects (two bound knots) behave as bosons.
- **Contrast with EX-GC10/11.** The tethered Gc particle had the 720° property but no exchange sign. Here both come
  from one object.
- **What is not shown:**
  1. A classical Z₂ of loops becomes an actual sign of an amplitude only once the model has a quantum layer (sums over
     histories with phases). This shows the topology needed for fermions is present; it does not show the particles
     are fermions.
  2. The knot is not yet a Gc particle. It does not consume G, so there is no gravity yet (Stage 3 of the design
     note).
  3. The loops were built by hand from smooth textures, not produced by the lattice dynamics. The lattice part shows
     they are allowed, not that the dynamics would choose them.

## Method
- `src/exgc13_loopclass.py` uses a Pontryagin–Thom construction:
  1. For a loop U_t, form V_t = U_t U_0⁻¹, a loop based at the identity.
  2. Take a regular value p near −1. Its preimage in (0,1) × R³ ≅ R⁴ is a set of closed curves, found from a grid
     search and traced by pseudo-arclength continuation.
  3. Frame each curve by J⁺ (the pull-back of a basis of T_pS³), with the tangent oriented so that det(T, F) > 0.
  4. Track the SO(4) frame (T, F) through its action on self-dual 2-forms, which gives an SO(3) loop, and lift it to
     SU(2) continuously.
  5. Score each component as (its π₁ class + 1) mod 2; the +1 is the bounding framing of a round circle. The loop
     class is the sum mod 2.

  Textures are compactly supported hedgehogs, F = π(1 − r²/R²)³ for r < R, with R = 2. Pairs use a product ansatz
  with separation 2d = 8. Supports are disjoint, so the configuration after a swap is identical to the start.
- `src/exgc13_lattice.py`: snaps each sampled configuration to the 48 lifts and replays every transition as
  single-node moves under EX-GC08's move table, with the new lift a short step from the old. It counts stuck or
  too-large changes, including the closing transition.
- Records: `loop_*.json` (24), `lattice_*.json`. Platform: Linux cloud container, Python 3.13.18, numpy 2.5.3. Compute:
  about 15 minutes.

## Limits and deviations
- The invariant is computed on smooth maps, not directly on lattice data. The lattice result rests on the
  straight-line homotopy argument, which is reasoned rather than proven.
- The first lattice attempt (spacing 0.5, then 0.25 for the swap) had steep textures: neighbouring lifts 90° apart and
  some stuck changes. The finer spacing 0.2, the same knot size as EX-GC12, is clean. Both are recorded here; only the
  clean runs are kept as files.
- The pseudo-arclength tracer needed a larger offset of p from −1 (0.5 rather than 0.25) for the 360° turn, where
  exactly −1 has a degenerate preimage. That change was made before any class was read off.
