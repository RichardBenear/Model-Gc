# Claude review of EX-GC35P, with next steps

2026-11-25 | Reviewer Claude | Study designed and run by ChatGPT | EXPLORATORY | PI CLOSED 2026-11-24.

Scope: an opinion on what 35P means and what to do next. I read the plan, note, HANDOFF and LOG on `main` at b13ee94 and recomputed nothing. The critical-pair examples below are worked by hand. Statements marked *inference* are my reading, not something a run showed.

## Summary

All five registered predictions held. Only the no-overlap control (A) passes both graph comparisons. Every rule with a real choice fails, even with site labels dropped: B, which picks a subset of tokens; C and D, which pick a batch; and E, which decides whether to capture. The reciprocal clock loop runs cleanly with the passing rule. That loop is valid but expected, because rule A has nothing to compete over.

*Inference:* the common cause is that in B–E the choice changes **how many events happen, or which inputs an event has**. D can send two tokens in one event or in two. E can capture, or let a token leave. C can batch an arrival now or later. Graphs with different event counts can never be isomorphic. So the whole family of "site holds a pool of tokens" rules is the wrong shape for causal invariance. Better scoring won't change that.

## What a passing interacting rule probably looks like

Wolfram-model rules that are causally invariant either have no overlapping matches or have overlaps that resolve into branches with the same causal structure. Gorard's paper, already cited in 35P, frames it this way. The 35P rules carry token *multisets* at sites. A different representation avoids the subset and batch choices: **each site is one element with one state**, such as vacant, R, L or particle. An event rewrites a fixed pattern of adjacent elements and consumes whole elements. Interactions then come from adjacency, not from shared pools.

Two hand-worked critical pairs show this can work:

- **Two-species exclusion with pass-through.** The rules are R0→0R, 0L→L0 and RL→LR. In R0L, the moves R0→0R and 0L→L0 compete for the same vacancy.
  - Branch 1 gives 0RL → 0LR → L0R.
  - Branch 2 gives RL0 → LR0 → L0R.
  - Both branches reach the same final state in three events. Both causal graphs are the same triangle: first→second, first→third, second→third.
  - They are isomorphic once labels are dropped, and related by the R↔L mirror when labels are kept.
- **A Model Gc absorber.** A particle element P absorbs G tokens from either side: GP→P and PG→P. In GPG the two absorptions share P. Either order gives a two-event chain, each event fed by one initial G. Unlabelled, the graphs are isomorphic. Labelled, they differ only in which G was absorbed first.

These are single critical pairs, not proofs. But they show interacting rules with real competition whose branches rejoin with the same causal structure. 35L–P found no such case. The absorber example also connects directly to the Gc mechanism (core-consumed carrier).

## Next steps (in priority order, rough costs)

1. **Site-element rule catalogue** (≤0.1 core-h; minutes). Reuse the 35L/35P scorer, with every site a single element. Candidates, with predictions to register before running:
   - (a) **One-species exclusion, 10→01:** no self-overlap. Predict pass.
   - (b) **Two-species exclusion with pass-through,** as above. Predict pass unlabelled, and labelled only modulo the mirror symmetry.
   - (c) **Two species without pass-through (RL blocked):** predict fail, ledger included. R0L jams as either 0RL or RL0, so the end states differ. This control shows that the exchange step is what lets (b) reconverge.
   - (d) **Absorber:** GP→P and PG→P plus G hopping into vacancies. Predict pass unlabelled.
   - (e) **Absorber with recoil:** absorbing from the left moves P right, as a discrete EX-GC16. Prediction open; this is the interesting one.
   - (f) **35L paired rule** as the failing control.

   Also add an **exhaustive critical-pair check** for each rule: all minimal overlaps, each branch completed and compared. For terminating rules this is the local test that suggests global causal invariance. Whether it is sufficient here should be stated as a check, not assumed.

2. **Emergent clocks instead of imposed deadlines** (design note, then ≤0.1 core-h). The 35I/K/O hosts impose the lapse through local deadlines, which removes branching altogether. In a causally invariant rule there are no deadlines. A site's "proper time" is the depth of its events in the causal graph. *Inference:* in an exclusion or absorber rule, crowded regions allow fewer events per causal step, so a lapse-like slowdown could emerge from occupancy. That would put gravity-like time dilation inside the causal graph, as the programme goal requires. Measure the events per unit causal depth at each site against the local density in rules (a), (b) and (d), with randomized schedules.

3. **Only after 1 and 2:** check whether the absorber rule (d or e) gives a steady G inflow toward P with a 1/r-like profile in 2D or 3D. That is the discrete analogue of MG195/EX-GC16. Then check whether two absorbers drift together. About 0.5 core-h.

4. **Keep paused:** knot-on-Y1 continuum work and the lens line, as recorded in HANDOFF. White paper v0.5 is frozen through 35O. 35P and the steps above would start the v0.6 record.

Total for steps 1–2: under 0.2 core-hours. Step 3 is about 0.5 core-hours.
