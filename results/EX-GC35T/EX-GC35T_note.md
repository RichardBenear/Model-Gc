# EX-GC35T — persistent-vertex absorption on a hypergraph

2026-11-24 | EXPLORATORY | ChatGPT | [Design and test plan](../../spec/EX-GC35T_plan.md). **PI CLOSED 2026-11-25.**

**Persistent-vertex absorption passes the anchored graph test when each carrier owns its consumed relation.** Adding a shared mass-counter relation restores destructive overlaps and graph failures.

Keep the particle as a vertex P. Each carrier owns an individually identified directed relation along its route, which terminates at P. A hop consumes that carrier's current relation and emits its next route relation. Absorption consumes the final carrier relation and appends an individually deposited mass relation on P. P is shared but never consumed. Count one mass relation per absorbed carrier, so free carriers plus deposits are conserved. Initial carrier relations and final deposits are fixed boundary ports. The geometry and routes are prescribed.

The finite screen used 1–4 carriers and every tuple of 1–3 hop counts, plus the four zero-hop direct controls. The positive screen exhaustively traversed the reachable progress states, testing every enabled event pair for a commuting diamond. Each complete graph is assembled from the resulting fixed per-carrier event identifiers and edges. Dynamic programming counted the serial histories represented by this state exploration. The positive family had no shared consumed relation and all 124 configurations passed both graph comparisons.

The negative control additionally consumes and rewrites one common counter relation at P on every absorption. Its state includes the absorption order, allowing the exhaustively represented terminal DAGs to be compared. Direct two-carrier absorption still gives an unlabelled chain, but its anchored input and deposit dependencies differ. With unequal route lengths (1,2), the plain graph already depends on which carrier absorbs first. In the routed screen 108 of 120 configurations fail the plain comparison and 117 fail the anchored comparison. The control therefore reproduces the serialization mechanism with a shared relation.

This demonstrates a parallel absorber design with a persistent shared vertex and independent consumptions. The passing cases do not resolve a destructive critical pair. No spatial inflow, clock calibration at an absorber, or attraction was run. [Code](../../src/exgc35t.py) · [State counts and graph checks](results.json)
