# EX-GC08 (EXPLORATORY): route 3, stage 1 – a discrete belt trick

Claude, 2026-11-19, at the PI's request ("yes, go ahead", Stage 1 of docs/ROUTE3_RIBBON_DESIGN.md). Exploratory
track; our construction.

## Answer
**Yes: a purely local graph rule reproduces the 720° topology.** A core tethered to a fixed outer boundary through a
3D lattice graph behaves like this:
- twisted by **one turn**, it can never be untangled (it settles at about 530 mismatched edges in 3 runs);
- twisted by **two turns**, it untangles completely (0 mismatched edges in 2 of 2 runs).

But the rule has to be of a specific kind. A rule that only looks at **pairs** of neighbours does not work: it either
lets the one-turn twist slip away or freezes everything. The rule must also look at **triples**: a node may change
its edge ordering only if the change does not swing any connection through the half-turn ambiguity. That is the
discrete version of "ribbons may not be cut".

## Set-up
On a cubic graph a node's rotation system (which incident edge counts as its +x, +y, +z) is one of the 24 rotations
of the cube. So the ribbon structure is a frame per node, carried by the graph's own edge ordering. A local move
changes one node's frame by one allowed step. The step sets considered were 90°, 120° and 180° rotations and their
unions.

## Part A: which local rules have the right topology? (algebra, exact)
The space of configurations reachable by local moves is a complex built on the 24 cube rotations. Its loops tell us
which twists can be undone.
- **Pairwise rules** (any step set S; all 15 unions of the 4 rotation classes were checked). Either
  1. the one-turn loop is a boundary, so the twist can be undone and the 720° distinction is lost (every S that
     contains 120° steps or a 180° class); or
  2. there are no or few moves that remove loops, so even two turns are stuck (90° only: H₁ has rank 49; 90° plus
     180° edge: rank 121).

  None gives the 3D-rotation topology.
- **Three-node ("no half-turn ambiguity") rule** with S = 90° + 120°. Every rotation smaller than 180° has a unique
  "short" lift to SU(2), the quaternion with positive scalar part. A move is allowed only if, for every neighbour,
  the short lifts compose correctly. Result: 384 of the 416 triangles survive. The homology is H₁ = Z₂ over GF(2) and
  0 over Q, which is the signature of the rotation group's Z₂. The one-turn loop is **not** a boundary. The product
  of short lifts around one turn is −1 and around two turns is +1.
- **The index is conserved by construction.** The sign of the product of short lifts along any path from the core to
  the boundary cannot change under allowed moves; it follows from the triangle condition in two lines. One turn has
  index −1, while the untwisted state (and any zero-energy state) has +1. So one turn can never reach zero energy.

## Part B: does the dynamics find the untangling? (annealing)
Cubic shell of half-size 10 (9261 nodes). The core (|x|∞ ≤ 1) and the outer face are held at the identity frame.
The initial twist is spread radially as steps of 90° about z. Metropolis annealing uses single-node moves under the
rule, and the energy is the number of edges whose two frames differ.

| Rule | Turns | Runs | Final energy (from) | Index |
|---|---|---|---|---|
| three-node, S = 90°+120° | 0 | 1 | 0 (0) | +1 constant |
| three-node | 1 | 3 (4000, 4000, 8000 sweeps) | 527, 528, 546 (3384) | −1 constant |
| three-node | 2 | 2 (4000, 8000 sweeps) | 0, 0 (7248) | +1 constant |
| pairwise only, S = 90°+120° | 1 | 1 | **0** (3384): the twist slipped away | changed |
| pairwise only, S = 90°+120° | 2 | 1 | 1865 (glassy, not topological) | changed |
| pairwise only, S = 90° | 2 | 1 | 3792 (frozen) | — |

The second two-turn run reached 368 after 4000 sweeps and 0 when rerun with 8000 sweeps (same seed). The 8000-sweep
record replaced the 4000-sweep one.

## What it means for route 3
- The 720° distinction can be built into a hypergraph rule set with **local** moves and **no extra field**. The
  ordering of each node's edges is the only extra structure.
- The physics is in the three-node rule. It forbids a connection from flipping through a half-turn in one step, which
  is exactly the discrete form of "ribbons cannot be cut". Without it the twist leaks away; with only nearest-
  neighbour constraints it freezes.
- The one-turn state carries an energy cost that does not anneal away (about 530 edge mismatches here). If this
  carried over to Gc, a particle turned once relative to its tether would hold a twist that only a second turn
  removes.

## Limits
- This is a test of the rule set's topology, not yet of Gc. There is no chemistry and no G flow; the tether is the
  whole shell. Stage 2 (Gc tokens on the ribbon graph) and Stage 3 (the tether formed by a particle's G supply) are
  still to come.
- Only the twist about z was tested, and the frames are restricted to the 24 cube rotations (steps of 90° and 120°).
  Finer frame sets are not needed for the topology, but they would smooth the energy.
- Annealing is a heuristic. The one-turn floor is guaranteed by the conserved index, but the exact value of the
  floor is not optimised.
- H₁ = Z₂ is consistent with the full loop group being Z₂, but it does not prove it. The two-turn runs reaching zero
  are the direct evidence that two turns can be undone.

## Method and records
`src/exgc08_belt.py` (Part A: `part_a`, `part_a2`; Part B: `part_b2`, with `z2_lift` as the index), records
`partA.json`, `partA2.json`, `B2_*.json`, `summary.json`; summary and figure `src/exgc08_summary.py` →
`paper/figures/EX-GC08_belt_trick.png`. Platform: Linux cloud container, Python 3.13.18, numpy 2.5.3. Compute: about
0.4 core-hours.

## Deviation from the log entry
The log predicted that a pairwise step set without 180° steps would keep the Z₂. It does not (Part A). The working
rule needs the three-node check, which was added before Part B. Part B was then run with both rules as planned, with
the pairwise rule as the control.
