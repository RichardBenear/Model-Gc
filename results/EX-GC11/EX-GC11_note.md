# EX-GC11 (EXPLORATORY): route 3, exchange sign – Step 0, what the frame-field tether can do

Claude, 2026-11-19, at the PI's request ("let's do the exchange sign test now"). Exploratory track; our construction.

## Answer
**In the tether we built for EX-GC10, an exchange cannot produce a sign. The test as first planned therefore cannot
work there, and a different kind of tether is needed.**
- Exchanging two untwisted cores leaves nothing behind. The tether returns to its thermal level (72 → 78 mismatched
  edges) and both Z₂ indices stay +1.
- A core carrying a one-turn twist cannot move at all. After the twist settles against the core surface (as in
  EX-GC10), no step of the core is allowed without cutting, even after 5000 relaxation sweeps.

## Why (argument written before the runs)
- **Classical spin–statistics needs tethers that cannot pass through each other.** The classical form (Dirac belts;
  Bolker 1973, "The spinor spanner"; Feynman's 1986 Dirac lecture) says that exchanging two objects tethered to a
  wall is equivalent to turning one of them by 360°. That relies on the belts having identity: a belt cannot pass
  through another belt or through an object.
- **Our tether is a frame field, so it has no such identity.** EX-GC10's tether is one frame per lattice node, and a
  "belt" is just any path through that field. Paths can be redrawn on the other side of a core at no cost, so tethers
  pass through each other freely. Moving a core by pure translation through an untwisted field changes nothing.
- **Hence no exchange sign.** Each core's index (the product of short lifts along a path from the core to the wall)
  does not depend on the path, because a loop around a core can be slid off it in 3D. Allowed moves conserve it, so no
  translation, and therefore no exchange, can change it.
- **The field version of spin–statistics needs solitons.** That version (Finkelstein & Rubinstein 1968) needs the
  particle itself to be a topological soliton of the field, a texture of degree 1 in π₃. Our twisted particle is
  "object plus twisted belt". Its index is a π₁ class tied to the fixed core, not a soliton winding.

## Runs (frames only, 25³ lattice, lift rule, T_f = 0.6, cores 3×3×3, half-circle radius 4)
| Case | Start: E, indices | End: E, indices | Note |
|---|---|---|---|
| exchange, both untwisted | 72, (+1, +1) | 78, (+1, +1) | no step ever blocked |
| A twisted by 1 turn, no exchange (control) | 1472, (−1, +1) | 542, (−1, +1) | the twist settles against A's surface |
| A twisted, then exchange | — | — | core A cannot take a single step (stopped after 500 retries × 10 sweeps) |
| A and B twisted, then exchange | — | — | same |

Code: `src/exgc11_exchange.py` (`Pair`, `try_step`: a core may enter a node only by an allowed move to the identity;
`lift_z`: index along +z). Compute: about 1 minute. Platform: Linux cloud container, Python 3.13.18, numpy 2.5.3.

## What it means
1. **Spin and statistics are decoupled in the frame-field tether.** It gives the 720° property (EX-GC08, EX-GC10) but
   no exchange sign.
2. **A twisted particle cannot move under the no-cut rule.** This is the conflict between "no cutting" and motion
   flagged in docs/ROUTE3_RIBBON_DESIGN.md §6(c). A twisted particle would need a move that carries its twist along
   with it.
3. **An exchange sign needs tethers with identity:** supply ribbons as individual paths that cannot cross each other
   or pass through a core. That is closer to the original ribbon-graph idea, but a new and larger model. The other
   option is particles that are themselves solitons of the frame field (route 2's kind of object).

## Deviation
The log entry was written right after the four quick runs (about 1 minute of compute), not before them.
