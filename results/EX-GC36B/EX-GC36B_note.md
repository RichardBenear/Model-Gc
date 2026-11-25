# EX-GC36B — order-free particles, odd vacuum forks, and a calibrated delayed-join test

2026-11-25 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC36B_plan.md). PI closure pending. Cloud container, Python 3.13. About 0.69 core-h of the 2 core-h cap.

**Only 1 of CRR's 80 closed-particle rules (66427) keeps its particle once the vertex-order condition is removed. Rule 55744 does not.** The repaired test now finds known delayed joins, including ones of unequal length. Under it, the published Wolfram rule's forks all rejoin. But the branches of CRR particles still do not reconverge. For 55744's particle that holds through 6 events per side. For 66427's order-free particle it holds through 5 per side, where the state cap was reached. In the 55744 vacuum, only 1 of 28 competing forks rejoined within 6.

## Step 1 — order-free closure (80 rules)

| Quantity | Rules |
|---|---:|
| Particle closed under version L with the published numbering (sanity check) | 80 |
| Successor set the same under all 120 renumberings (version L) | **0** |
| Every phase closed under version C (no vertex order) | **1** (rule 66427) |

- **55744:** under C, phases A and B each have 4 moves, of which 2 stay particles. Phase C has 3 moves, of which 1 stays. **Q1 holds.**
- **ChatGPT's count is reproduced:** under L, phase B gives 1, 2 or 3 successors, 40 renumberings each. Across renumberings there are 8 different successor sets.
- **66427:** phases with 2 and 3 order-free moves, all of which land in the particle class. Under L, by contrast, some renumberings leave it with *no* move.
- **Q2** (my expectation of fewer than 10) holds.
- **Separate path** ([exgc36b_recompute.py](../../src/exgc36b_recompute.py)): CRR's own `x6_generic` helpers with an independent orientation-free gate search. It also gives exactly [66427]. The first version of this check wrongly used `normalize_subset`, which sorts the support and so ignores the permutation. It returned no rules. I fixed it before recording this result.

## Step 2 — the four odd 36A vacuum forks

All 4 are class (ii), so **Q3 holds**. Both orders are legal. In each case one event's output already exists, so it is suppressed. In the other order, the first event consumes that relation, so the output is added. The final live sets differ by one relation. Duplicate-output suppression, a rule detail in `x6_generic.py`, makes otherwise independent events fail to commute.

## Step 3 — state-frontier comparator F (D = 6)

**Gate G3 holds:**

| Witness | F | 36A equal-length history test |
|---|---|---|
| chains m = n = 1…6 | join at (m, m) | join at m |
| (2, 4), (3, 5) | join at (2, 4), (3, 5) | no join (as expected) |
| negative (different terminals) | terminal, distinct | terminal, distinct |

**Applications:**
- **A1, published `xy,zy→xw,yw,zw` on the asymmetric seed:** all 3 forks state-join, at (0,0), (1,1) and (1,1). **Q4 holds.** This is Wolfram's branch-pair criterion, which compares states. The 36A failure of this control therefore came from its stricter history comparison: live-relation producers in equal-length histories. *Inference:* 36A's "unresolved" verdicts were too strict a reading.
- **A2, 55744-L particle fork at B:** unresolved through 6 events per side.
- **A3, 55744-L vacuum, 28 overlapping forks from the 36A start-state multiway:** 1 joins, at (6, 6), from depth 4. The other 27 are unresolved through 6 per side.

## Step 4 — rule 66427 (order-free) particle forks, version C

| Phase | Moves | Fork | F |
|---|---:|---|---|
| 0 | 2 | overlapping | no join through 5 per side, then the cap of 20,000 states |
| 1 | 3 | overlapping | no join through 5 / 4 per side, then cap |
| 1 | 3 | overlapping | no join through 5 per side, then cap |
| 1 | 3 | non-overlapping | join at (1, 1), as expected |

## What it means (inference)

- **Vertex order is load-bearing in CRR's particles.** 79 of 80 lose closure without it, so the particle exists only on ordered hypergraphs. 66427 is the single exception, and the natural candidate for any further causal-invariance work.
- **With a calibrated test, the published rule resolves and CRR particles do not, within these bounds.** The particle's internal branching keeps it a particle but sends the world down histories that do not meet. Bounded results are not proofs. But the contrast with the published control, which joins within 1 step, is now a like-for-like comparison.
- For the programme, this is the same pattern as 35L–35X: persistence comes with branches that don't rejoin. In the rule family CRR found, we still have no rule where a particle and causal invariance coexist.

## Files

[Code](../../src/exgc36b.py) · [Recompute](../../src/exgc36b_recompute.py) · [results.json](results.json) · [recompute.json](recompute.json)
