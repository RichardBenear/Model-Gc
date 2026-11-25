# EX-GC36 design note: is CRR Rule 55744 causally invariant, and can it carry a Gc particle?

2026-11-25 | Claude | DESIGN ONLY, nothing registered | Requested by the PI after the [35U–X review](../review/EX-GC35UVWX_claude_review.md).

## Why this study

- 35L–35X found no small interacting rule whose competing histories rejoin. Rules pass only when nothing competes, and fail once particles move and interact on a fixed lattice (35W).
- The CRR investigation (RichardBenear/CRR-Investigation) is already the "space is the graph" route. It has one rule, and a particle is a recurring loop whose identity lives in relations, not vertices. Its Rule 55744 has a **persistent gate** relation, the same device as our 35T persistent vertex. Its particle class also **branches internally and stays itself**: A→B, C→B, B→A or B→C.
- CRR never tested causal invariance. It updated with one global random choice per event. Its Y165 two-particle "interaction" was traced to that shared event sequence, which is exactly an order dependence.

So 55744 is the best candidate we have for real branching plus a persistent particle. The 35-series tools are the missing check for it. Model Gc is our own construction, and so is any Gc reading of 55744. CRR's rule is Richard's earlier construction, not LaViolette's.

## The rule (from `external/crr/x6_generic.py`, CRR Stage X6, unchanged)

- **Relation:** a ternary relation T(a; {b, c}), with a source a and an unordered target pair.
- **Match:** five distinct vertices x, y, z, u, v with u < v, and four relations:
  - T(x;{y,z}), T(y;{z,x}) and T(x;{u,v}), which are consumed;
  - T(v;{y,z}), the gate, which is read but not consumed.
- **Output:** one fresh vertex w and eight relations, from the template [[1,4,5],[2,1,3],[2,4,5],[3,1,4],[3,2,5],[3,4,5],[4,2,5],[5,2,3]], with roles 1–5 = y, z, u, v, w.
- **Start state:** T(0;{1,2}), T(1;{2,0}), T(0;{3,4}), T(4;{1,2}).
- **Update:** one legal match chosen uniformly at random per event.

**A load-bearing detail.** The tie-break u < v uses vertex numbers, and vertex numbers are birth order under the *global* event counter. Two vertices born in causally unrelated events have no order of their own. Under a different schedule their numbers can swap, which changes which output roles they get. *Inference:* this alone could break causal invariance, independently of any overlap. The study must test two versions:
- (L) the rule exactly as written;
- (C) the tie-break replaced by a schedule-independent one: causal order where it exists, and a symmetric treatment of {u, v} (both role assignments are separate matches) where it does not.

## Pre-design check (run, cloud, ~4 min; [script](../external/crr/anchor_check.py))

| Quantity | W37 row 0 (CRR record) | This check (Python 3.13, cloud) |
|---|---:|---:|
| Nodes at 5,000 events | 5,005 | 5,005 |
| Relations | 21,131 | 21,127 |
| Max degree | 402 | 638 (2-section) / 1,327 (incidence) |

**The reproduction anchor does not match.** The same simulator and seed give a different history here. The likely causes are:
- (a) `rng.choice(tuple(gates_set))` depends on Python's set ordering, which can change between Python versions;
- (b) W37 came from the Stage W simulator rather than `x6_generic.py`;
- (c) a different degree definition.

This needs settling before any result is compared with CRR's records. It does *not* block a causal-invariance test, which compares all schedules and does not depend on any one random stream.

The same check counted competition. Over the first 2,000 events there were on average about 600 legal gates per state (minimum 1, maximum 1,183). Of 79.1 million gate pairs:
- 105,707 (0.13%) share a consumed relation, which is a genuine overlap;
- 124,645 (0.16%) have one gate's persistent relation consumed by the other, a read/write overlap like 35W.

So 55744 has real competition, unlike every rule that passed in the 35 series.

## Questions and stages

**Stage 0. Anchor** (no compute here; Richard's PC or the original environment). Reproduce W37 row 0 with the CRR Python version, or identify which simulator produced it. Record the platform.

**Stage 1. Calibrate the test** (≤0.2 core-h). This is step 1 of the 35U–X review. Score published causally invariant rules on asymmetric seeds under two conventions:
- (i) unlabelled states and causal graphs, the Wolfram definition;
- (ii) the anchored local certificate.

Stages 2–3 use whichever convention passes the published positive controls, and report the other alongside.

**Stage 2. Vacuum multiway test** (≤0.2 core-h). From the start state, enumerate all schedules to depth 6–8. Merge states by hypergraph isomorphism. Then:
- Do all branches reach isomorphic states at equal depth?
- Are the causal graphs to a common depth isomorphic?

Run L and C. Also take the ~0.3% of overlapping gate pairs from a grown vacuum, around 2,000 events, and resolve each critical pair locally to 3–4 further events.

**Stage 3. Particle test** (≤0.1 core-h). Use the phase masks from `W21_CLOSED_PARTICLE_HITS.json`, which is needed from the CRR files.
- For the particle's internal fork at B (B→A versus B→C), do the two branches rejoin, with matching causal graphs, when completed through one period?
- Does the same hold for a particle–vacuum contact fork?

The particle class closing (CRR's result) is a weaker statement than this. It says both branches stay particles, not that their histories match.

**Predictions (to be frozen before running).** *Inference:* I expect L to fail at Stage 2 through the birth-order tie-break. I expect C's outcome to be open. I expect Stage 3 to be the most likely place for a pass, because CRR already found the branching closes on the particle class.

**Stage 4. Only if Stage 2 or 3 passes under C: a Gc reading** (1–2 core-h). In 55744 every event consumes 3 relations and makes 8, so vacuum growth is a flux of fresh relations. Does a particle, which recycles its own relations, act as a sink for that growth? Measure growth density or ball-growth against graph distance from the particle, with a no-particle control, and look for a deficit that falls off with distance. Caveat from StochasticVacuum's Rule 55744 plan, which tabulates CRR Stage W: the dimension estimate was still rising (about 3.6 to 4.7 from 5k to 80k events) and hub degree grows, so distances may be small-world. Any falloff has to be read against that.

## What each outcome would mean

- **Stage 2 or 3 passes under C:** a rule with real branching, a persistent particle and causal invariance. This is the combination the programme needs for both the gravity and QM links. It would lead straight to Stage 4.
- **Both fail under C:** CRR-style particles depend on event order. That supports path A (the deterministic Y1 event graph) for gravity work, and is a result the CRR paper should state.
- **L fails but C passes:** the birth-order tie-break was the only obstruction. That is worth reporting back to CRR on its own.

## Inputs committed

`external/crr/` holds CRR files as uploaded by Richard on 2026-11-25, unchanged:

| File | SHA-256 |
|---|---|
| `x6_generic.py` | 827cf676fe04f4d26ca39213173055baf2cc8ce8a639bf5b2d2cfc9eef523f77 |
| `W37_RULE55744_MULTI_SEED_5K.csv` | 7d523fb0b874b4f1abf2d1e9a3f99e5f080e15b0af7e8f0c07fcecd0c67ec8ea |

Still needed from the CRR files: `W21_CLOSED_PARTICLE_HITS.json` (particle phase masks) and, for Stage 0, the CRR Python version or the Stage W simulator.

Total for Stages 1–3: about 0.5 core-hours, cloud. Stage 4 is 1–2 core-hours, conditional.
