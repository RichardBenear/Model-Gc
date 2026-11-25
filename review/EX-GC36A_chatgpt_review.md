# Independent review of EX-GC36A — Rule 55744

ChatGPT | [Results](../results/EX-GC36A/EX-GC36A_note.md) · [Registered plan](../spec/EX-GC36A_plan.md)

**GC36A reproduces the particle transition under the published rule. Causal invariance remains open.** The strongest finding is that changing vertex numbering while preserving the underlying unlabelled phase-B hypergraph can change its legal successors.

## Independent one-step check

I reconstructed the ten relations of isolated phase B from the recorded consumed and read relations. I enumerated the five distinct vertex roles under each of the 120 possible vertex renumberings, applied version L including its `u < v` gate and duplicate suppression, and counted distinct successor relation sets.

| Distinct one-step successors | Vertex renumberings |
|---:|---:|
| 1 | 40 |
| 2 | 40 |
| 3 | 40 |

The published numbering gives two successors, each with 14 relations. Mapping vertices (0,1,2,3,4) to (0,1,3,4,2) admits an additional successor with 15 relations. Its relation count distinguishes it even after erasing vertex labels. If vertex order is physical state, the published transition remains valid for that ordered state. This check does not settle causal invariance of that ordered system. [Implementation](../src/exgc36a.py) · [Phase-B match data](../results/EX-GC36A/results.json)

## Interpretation

- **Particle persistence and branch reconvergence are separate.** Both version L successors contain recognized particle phases. The independent CRR-simulator replay found no isomorphic states through four further events at matching event counts. Unequal-length joins and longer continuations remain open. Version C changes the transition relation by admitting both orientations. It tests that alternative rule.
- **The vacuum results are descriptive.** The positive control joined only one of three asymmetric forks within the registered horizon. The symmetric star checks immediate isomorphism and supplies little evidence about delayed joins. The unfinished vacuum stages further limit the verdict. [Replay code](../src/exgc36a_recompute.py)

## Priorities

1. Calibrate delayed joins against a known witness. Then compare state frontiers across different event counts before spending compute on causal-history isomorphism.
2. Replay the four vacuum forks classified as non-overlapping and trace both orders of the two events. Check whether a match disappears, the fresh numbering changes eligibility, or output suppression changes the live relations or their producers. Duplicate-output suppression is a testable hypothesis. [Join and fork implementation](../src/exgc36a.py)
3. Treat particle persistence in an ordered rule as the current Model Gc result. Test separately whether persistence survives without vertex ordering determining which branches exist.
