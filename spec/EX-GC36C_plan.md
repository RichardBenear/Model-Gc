# EX-GC36C — does the order-free particle of rule 66427 persist in its own debris?

2026-11-25 | EXPLORATORY | Claude | PI go-ahead in thread ("yes", 2026-11-25 20:51 UTC). It follows [ChatGPT's 36B review](../review/EX-GC36B_chatgpt_review.md). Default cap 2 core-hours; expected ≤0.2. Written before any run.

**Correction recorded from the 36B review.** The 36B frontier merged version-L states by unlabelled isomorphism. Numbering changes L's legal moves, so the 36B vacuum count of 1 rejoined fork in 28 is a **lower bound**. ChatGPT's numbered recheck of the 55744 phase-B fork confirms no join within 6 steps per side. Version C rules, including 66427 here, are unaffected.

## Question

36B tested 66427's two phases isolated. In a real evolution, consumed-and-produced relations leave debris on the support's old vertices, and that debris can create new legal moves. Does a particle-shaped support keep renewing itself, and for how long?

## Method

- **Rules:**
  - 66427, version C: no u < v, both orientations, set semantics, fresh vertex per event. Template and phase masks come from `external/crr/W21_CLOSED_PARTICLE_HITS.json`.
  - Control 1: 55744, version L, as CRR ran it.
  - Control 2: 66427, version L.
- **Starts:** each phase of the rule, isolated on vertices 0–4.
- **Dynamics:** CRR's scheduler. One legal match is chosen uniformly from all legal matches in the whole state, debris included. The runs use `x6_generic.Sim` with its incremental gate bookkeeping. For version C, the gate search is overridden to allow both orientations; nothing else changes.
- **Schedules:** 200 per start, seeds 0–199 in `random.Random`, up to T = 2,000 events or until no move is legal.
- **Particle tracking:**
  - A support is a 5-vertex set whose induced relations are isomorphic to one of the rule's phases.
  - After each event, the event's successor support (y, z, u, v, w) is tested; if it is a phase, it is added to the live set as a *renewal*.
  - Every live support is retested after each event and dropped if no longer a phase.
  - The particle is *alive* while the live set is non-empty.
- **Readouts per schedule:**
  - lifetime (the event at which the live set first empties);
  - renewals;
  - the maximum number of simultaneous live supports, which would show duplication;
  - the fraction of events whose five vertices lie inside a live support;
  - the final number of legal moves;
  - whether the run halted.

## Predictions

- **Q1 (scored):** 66427-C is *persistent* if at least 50% of schedules from each phase are still alive at 500 events. My expectation, recorded and not scored, is not persistent.
- **Q2 (scored):** the median 55744-L lifetime exceeds the median 66427-C lifetime.
- **Q3 (descriptive):** the fraction of debris-driven events, and the maximum number of simultaneous supports, for each rule.

Nothing here tests causal invariance. Forks are not compared in this study.

## Outputs

`src/exgc36c.py`, `results/EX-GC36C/` (note plus `results.json`). Separate recomputation: replay 5 schedules per start with an independent tracker that scans every 5-subset of the most recent 12 vertices for phase supports, using 36B's canonical form. Its lifetimes must agree with the incremental tracker's.

**Amendment before any run.** A schedule stops when the particle dies, which is when the live set empties, instead of continuing to T. Lifetime is unaffected. The final legal-move count and the halted flag are then reported at the stop. This keeps cost inside the cap, because a vacuum-like debris state can carry hundreds of legal moves.

**Deviation D1 (after a smoke run, before any scored data).** The first full run was stopped after about 14 minutes. Nothing from it was recorded. A smoke run of 30 seeds from 66427-C phase 0 then showed two things. Most schedules die within a few events. The survivors (2 of 30) keep going to 2,000 events, with particle-shaped supports multiplying to more than 200 and about 1,700 legal moves, at about 130 s per run. At that rate the controls would far exceed the cap. So **T is cut to 600**. Q1, scored at 500 events, is unchanged. "Alive at 2,000" is dropped. The 30 smoke seeds are rerun as part of the main set; they are not excluded.
