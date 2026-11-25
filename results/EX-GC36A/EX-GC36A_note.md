# EX-GC36A — causal-invariance audit of CRR Rule 55744

2026-11-25 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC36A_plan.md) · [Design](../../docs/EX-GC36_rule55744_causal_invariance_design.md). PI closure pending. Cloud container, Python 3.13, networkx 3.7. About 1.6 core-h of the 2 core-h cap.

**No competing pair of events in Rule 55744 rejoined within three further steps, in the vacuum or in the particle. But the calibration gate was not met, so this is a bounded, descriptive result, not a finding that 55744 is non-invariant.** The particle's two branches stay particles, as CRR found. Their histories do not meet, though: no isomorphic states appear for up to four further steps on each side. CRR's particle closure also depends on the rule's birth-order tie-break. Remove that tie-break and two new events appear in phase B, and both leave the particle class.

## Stage 1 — calibration (gate P1 not met)

| Rule and seed | Forks | Join under U (unlabelled) | Join under A (anchored) |
|---|---:|---:|---:|
| `xy→xy,yz`, two loops | 1 (no overlap) | 1 | 1 |
| `xy,zy→xw,yw,zw`, 3-in star + leaf edge | 3 (overlapping) | **1** | 0 |
| `xy,zy→xw,yw,zw`, 4-in star | 15 (12 overlapping) | 15, all at step 0 by symmetry | 3 (only the non-overlapping) |
| `xy,xz→xy`, 2-out star + leaf edge (negative) | 1 | 0 (terminal, distinct) | 0 |

The negative control fails as it must. But the published overlapping rule joins only 1 of its 3 forks on the asymmetric seed within k = 3, so **P1 fails**. Per the plan, stages 2–3 are descriptive only. The separate path through ChatGPT's 35V engine gives the same 1-of-3 ([recompute.json](recompute.json)). Two readings remain open, and this run cannot tell them apart:
- k = 3 is too short for this rule's resolutions;
- the published claim holds for the rule's own initial conditions, not for arbitrary seeds.

## Stage 2a — vacuum, start-state multiway to depth 4, version L

| Fork type | Forks | Joined (U) | Unresolved at k = 3 (U) | Anchored (A) |
|---|---:|---:|---:|---|
| Overlapping (shared consumed or read/write) | 28 | 0 | 28 | same as U |
| Non-overlapping | 36 | 32 (all at step 1) | 4 | same as U |

**P2 holds** in the bounded sense: under L, no overlapping fork joins. I did not diagnose the 4 non-overlapping forks that fail to commute. A likely cause is set semantics, where an output already present is not added again, so one event can change what the other produces. That is *inference*.

Version C (2a) and both versions of 2b (the grown-vacuum critical pairs) were **not completed**:
- 2a-C reached the 25-minute command limit with no output (deviation D1);
- I stopped 2b-L at about 0.17 core-h to stay inside the cap;
- 2b-C was not started.

So P3 is unscored.

## Stage 3 — particle fork at phase B

| Version | Events in isolated B | Successor phases | B→C vs B→A fork, U | Same fork, A |
|---|---:|---|---|---|
| L (as written) | 2 | C and A (closure reproduced, **P4a holds**) | unresolved at 3 | unresolved at 3 (**P4b holds**) |
| C (no u < v) | 4 | C, A, and two events whose successor is **no particle phase** | unresolved at 3 | unresolved at 3 |

In version C, the extra non-overlapping fork joins at step 1, as it should. The other five forks are unresolved.

**P4c** (the main question): unresolved. The separate path used the CRR simulator for the dynamics and plain hypergraph isomorphism of states. It found no isomorphic left and right *states* at any depth up to 4 further events (left 65, right 9 states at depth 4). State confluence is weaker than causal-graph equivalence, so the histories do not rejoin in this window.

## What it means (inference)

- 55744 has genuine competition, and in every case tested here the competing branches went separate ways for at least 3–4 steps. Nothing so far supports causal invariance.
- The particle's "branching that stays itself" holds at the level of the particle's state class, not of its history. Different branches leave different leftover structure behind.
- The closure itself uses u < v by global vertex number, which is a non-local, schedule-dependent ingredient. That is worth stating in the CRR record.
- We still cannot certify a known positive control with this bounded test. Until that is fixed, "unresolved" here cannot be read as "not causally invariant".

## Files

[Engine and stages](../../src/exgc36a.py) · [Recomputation](../../src/exgc36a_recompute.py) · [results.json](results.json) · [recompute.json](recompute.json). Inputs in `external/crr/` (SHA-256s in the design note). The engine reproduced the CRR simulator's legal events and outputs for 40 steps of seed 0 before use.
