# Route 3: a 720° (spinor-like) particle from the structure of a Gc hypergraph (design note, 2026-11-19)

Author: Claude, at the PI's request (Richard Benear, 2026-11-19: "let's do route 3"). Design and analysis only; no
runs. Each point is marked **established** (standard mathematics or physics, or our earlier records), **design**
(a choice we would be making), or **speculative**. The model is our construction, not LaViolette's.

## Summary
A 720° property needs the particle to be **tethered** to its surroundings by connections that remember twisting.
In a hypergraph, the tether can be the graph's own connectivity. What the graph must also carry is an **ordering of
each node's connections** (a "rotation system" or ribbon structure), so that "twisted" has a meaning. Gc already
supplies a natural tether: every particle is fed by a steady inflow of G from the far reservoir, and in a token graph
that inflow runs along definite edges. This answers the "what creates the twist in reality" objection to route 2:
nothing extra is created. A rotation of the particle twists the supply lines it already has.

I propose three stages, each a separate exploration:
1. **Discrete belt trick.** Do our local rewrite moves on a ribbon graph reproduce the 720° topology (one turn cannot
   be undone, two turns can)?
2. **Gc on a ribbon graph.** Tokens and conserved G on the same graph: does a particle form, and is the far field
   still 1/r?
3. **Tethered Gc particle.** Rotate a Gc particle and test whether its G-supply tether shows the 720° behaviour, and
   whether exchanging two particles gives the matching sign.

## 1. Why a tether is the essential ingredient (established)
- The rotations of 3D space form a group whose loops come in two classes: a 360° turn is a loop that cannot be
  shrunk away, while a 720° turn can (π₁(SO(3)) = Z₂). An isolated object never notices this: after 360° it is
  simply back where it started.
- It becomes physical when the object is attached to fixed surroundings by flexible connections. A 360° turn leaves
  the connections tangled in a way no rearrangement can undo. After 720° they can be untangled without moving the
  object. This is Dirac's belt or plate trick, and it works for any number of connections, even a continuous
  spread of them in every direction.
- In quantum theory, spin-½ wave functions change sign under 360°. The tethered object is the classical picture of
  that sign: the tangle after one turn is the "minus".

## 2. What a graph needs to represent twisting (design)
A bare hypergraph has no notion of rotation: it only records which nodes are linked. Two standard additions supply
one:
- **Rotation system (ribbon graph).** Each node stores the cyclic order of its edges, and each edge carries a twist
  bit (0 or 1, a half-twist). This is how graph theory describes graphs drawn on surfaces, and it is the discrete
  form of a "framing". Edges become thin ribbons, and twisting a node's star of edges is a well-defined operation.
- **Local frames.** Alternatively, each node carries a discrete orientation (one of the 24 rotations of a cube), and
  each edge carries the relative rotation between its ends. This is closer to lattice gauge theory and easier to
  check numerically, but it is nearer to route 2's "field", so the ribbon version is preferred.

Either way, the twist lives in the *connections*, not in an extra substance. For a hypergraph model this is natural,
because connectivity is all there is.

## 3. The tether in Gc (design, motivated by our results)
- In Gc every particle consumes G at a steady rate M, and the G it consumes arrives from the far reservoir (EX008,
  MG195; Gauss's law is G conservation, see docs/PHYSICS_MAPPING.md). In a token graph that supply is a set of token
  paths from the boundary into the core.
- If the edges those tokens travel along are ribbons, rotating the core rotates the inner ends of the supply ribbons
  while the outer ends stay fixed at the reservoir. That is exactly the tethered-object geometry of §1.
- **The key design question is reconnection.** In the continuum, G diffuses, so flux lines have no identity and any
  twist relaxes instantly. In a graph, a twist survives only if the rewrite rules **cannot cut and re-glue ribbons
  freely**. They may slide ribbons and pass twists along, but they may change the total twist only by a full turn at
  a time. With that restriction the 360° tangle is protected; without it, it is not. This is a rule-design choice we
  would make deliberately and then test.

## 4. Stage 1: a discrete belt trick (proposed exploration, cheap)
- **Set-up.** A core node-cluster joined to a fixed outer boundary by N ribbons (N = 3, 6, 12, and a dense spherical
  shell), embedded in a 3D lattice for bookkeeping only. The allowed moves are local ribbon moves: slide a ribbon over
  a neighbour, pass a half-twist along a ribbon, and swap two adjacent ribbons with a compensating twist. No cutting.
- **Test.** Apply a 360° turn of the core (a cyclic shift of its rotation system plus the matching twists). Search for
  a sequence of allowed moves that returns the original graph while the core is held fixed: breadth-first search for
  small N, and simulated annealing on a "tangle energy" (total twist plus crossings) for large N. Repeat for 720°.
- **Prediction.** 360°: no sequence exists (energy stays above a floor). 720°: a sequence exists (energy returns to
  zero). Also predicted: a Z₂ invariant (total twist mod 2, suitably defined) that flips at each turn and that no
  allowed move changes.
- **Why it matters.** This checks that our local rules encode the topology of 3D rotations. If they don't, nothing
  built on them can give spin ½. It is mathematics more than physics, but it fixes the rule set for Stages 2 and 3.

## 5. Stage 2: Gc on a ribbon graph (proposed, larger)
- Port Gc1 to tokens on a 3D-like ribbon graph. Species are counts on nodes; hops go along edges; G is conserved
  exactly, with a reservoir boundary. The parent ModelG tick-based token engine (MG104+) and the StochasticVacuum
  token reservoirs are the starting code (docs/RELATED_REPOS.md).
- Questions: the minimum token number Ω at which a particle survives (open since the Gc set-up), and whether the far
  G field is 1/r on the graph. There is no 1D Gc particle (EX008), so this must be 3D-like from the start.

## 6. Stage 3: the tethered Gc particle (proposed, the actual test)
- Track supply-token paths as ribbons, rotate a particle (or let a rotating one form), and apply Stage 1's test to
  the supply tether. Then exchange two particles: in 3D, exchanging two tethered objects is linked to a 360° turn of
  one of them (the spin–statistics connection, which has a classical tether version). This would test whether a
  "minus sign" under exchange comes with the 720° property, which would be the fermion signature.
- **Speculative caveats.** (a) A steady Gc particle does not rotate by itself (EX-GC05–07), so Stage 3 would rotate
  it by hand, or ask whether rotation costs energy that a twist stores. (b) Whether a classical tangle ever becomes a
  quantum sign depends on the programme's longer-term link to quantum mechanics, so Stage 3 can at most show the
  geometry is there. (c) If the rules forbid cutting ribbons, G supply paths can never be rerouted. That may conflict
  with particles moving, which needs paths to change, so Stage 2 must check that the restriction still allows motion
  and 1/r.

## 7. Recommendation
Start with Stage 1 as EX-GC08. It is cheap (graph search and annealing on small graphs, well under 1 core-hour) and
decides the rule set before any chemistry is ported. Stages 2 and 3 are larger and would each get their own entry.

References: Dirac belt and plate trick (e.g. Penrose & Rindler, Spinors and Space-Time, vol. 1, §1.5); ribbon graphs
and rotation systems (Gross & Tucker, Topological Graph Theory, 1987); framed links and ribbon categories (Kauffman,
Knots and Physics, 1991); Bilson-Thompson, "A topological model of composite preons" (arXiv:hep-ph/0503213, 2005);
Finkelstein & Rubinstein, J. Math. Phys. 9, 1762 (1968).
