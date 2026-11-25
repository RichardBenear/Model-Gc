# EX-GC36B — order-free particles, odd vacuum forks, and a calibrated delayed-join test

2026-11-25 | EXPLORATORY | Claude | PI go-ahead in thread ("yes", 2026-11-25 19:25 UTC) to the revised plan that followed [ChatGPT's 36A review](../review/EX-GC36A_chatgpt_review.md). Default cap 2 core-hours; expected ≤1. Written before any test below is run.

## Why

ChatGPT showed that phase B of Rule 55744 has 1, 2 or 3 successors depending only on vertex numbering. So version L is a rule on *ordered* hypergraphs, and unlabelled causal invariance is ill-posed for it. The questions here are:
- whether any CRR closed-particle rule keeps its particle without a vertex order;
- what makes 4 non-overlapping 36A forks fail to commute;
- whether a corrected test can find delayed joins at all.

## Inputs

`external/crr/W21_CLOSED_PARTICLE_HITS.json`, with 80 rules (66 of period 2, 14 of period 3) and their phase masks, and `external/crr/x6_generic.py`. All 80 rules share the x6 match pattern and differ only in the 8-relation output template.

## Step 1 — order-free closure screen (all 80 rules)

- **Particle class:** the unlabelled isomorphism classes of the rule's phase relation sets, each on 5 vertices. The induced relations on a 5-vertex support are compared by hypergraph isomorphism, with source and target slots distinguished.
- For each phase, taken isolated:
  - **(a) L-order dependence.** Over all 120 renumberings of the 5 vertices, apply version L (u < v by number) and record the set of successor isomorphism classes. A phase is *order-free under L* if every renumbering gives the same set.
  - **(b) C closure.** Apply version C (no u < v, both orientations as separate matches). The phase is *C-closed* if it has at least one move, and every move's successor support (y, z, u, v, w) carries a phase of the same rule up to isomorphism.
- A rule *passes* if every phase is C-closed. It *passes strictly* if it also has every phase order-free under L.
- **Predictions:**
  - Q1: 55744 fails (b), as in 36A, where 2 of 4 C moves from B leave the class. It also fails (a), reproducing ChatGPT's 1/2/3.
  - Q2: the number of passing rules is open. My expectation, recorded and not scored, is fewer than 10 of 80.

## Step 2 — the four odd vacuum forks of 36A (2a, version L)

Replay both orders of each of the 4 forks that are non-overlapping yet unresolved, recomputed from the 36A code. Classify each:
- (i) the second event is no longer legal after the first;
- (ii) both orders are legal, but duplicate-output suppression gives different live relations or producers;
- (iii) fresh numbering changes u < v eligibility;
- (iv) other.

**Q3:** (ii) explains at least 3 of the 4 (ChatGPT's hypothesis, which I share).

## Step 3 — calibrated delayed-join test

**New comparator F (state frontier).** For a fork, collect the unlabelled isomorphism classes of all states reachable from each branch in 0..D further events, with D = 6 and a cap of 20,000 states per side per depth. The fork *state-joins* if the two frontiers share a class at any pair of depths, so unequal lengths are allowed. Report the smallest total depth. State confluence is necessary for causal invariance, not sufficient.

**Known witnesses.** These use typed unary relations. A fork R(x)→P1(x) versus R(x)→Q1(x), followed by chains P1→…→Pm→Z(x) and Q1→…→Qn→Z(x).
- Equal lengths m = n = 1..6: F must report a join at depths (m, m), and the 36A equal-length history test must report a join at level m.
- Unequal (m, n) = (2, 4) and (3, 5): F must report a join at (m, n). The equal-length history test is expected not to.
- **Negative witness:** chains ending in different terminal types Z and Z′. F must report no join.
- **Gate G3:** all witness outcomes come out as stated. If not, Step 3's other results are descriptive.

**Applications of F, with D = 6:**
- (A1) The three forks of `xy,zy→xw,yw,zw` on the asymmetric 3-in star + leaf edge.
- (A2) The 55744-L particle fork at B.
- (A3) The 28 overlapping 55744-L vacuum forks from 36A 2a, run up to the cap.
- **Q4:** A1 joins all three forks within D = 6. This is the published claim; if it fails, I record it, not explain it away.
- **Q5:** A2 and A3 are open.

## Step 4 — conditional

Only for rules that pass Step 1, and only if G3 holds and budget remains: apply F with D = 6 to each passing rule's internal particle forks, under version C. Otherwise report "not run".

## Outputs

`src/exgc36b.py`, `results/EX-GC36B/` (note plus `results.json`). Separate recomputation for Step 1: recount C-closure with the CRR simulator's own `abstract_phase_graph` logic, modified only to drop u < v, against the new code. Platform recorded.
