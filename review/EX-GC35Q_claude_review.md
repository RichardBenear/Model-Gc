# Claude review of EX-GC35Q and the step-2 entry gate, with next steps

2026-11-25 | Reviewer Claude | Study designed and run by ChatGPT | EXPLORATORY | PI closure of 35Q pending.

Scope: what 35Q means, what the step-2 stop means, and what to do next. I read the plan, note and HANDOFF on `main` at 1260606. I recomputed nothing. The bordered critical-pair argument below is worked by hand. Statements marked *inference* are my reading, not something a run showed.

## Summary

- **35Q was step 1.** It ran the catalogue I proposed, with predictions written down before the run. My two hand-worked examples, pass-through (B) and the absorber (D), resolve on their minimal overlaps, as I predicted. But both fail plain causal-graph equivalence in slightly larger words: 0R0L and G0PG. **My review's prediction that B and D would pass was wrong.** The single-pair check I relied on is too weak, and 35Q showed it directly.
- **Step 2 was not a failure.** It stopped at the entry gate the design asked for, and no clock run took place. One correction to the HANDOFF reading, though. Rule A (one-species exclusion, 10→01) *is* interacting: tokens block each other and jam. It is causally invariant because its matches never overlap, not because nothing interacts. So step 2 can run on A now. See next step 2.

## What 35Q means

**Why the local checks misled (inference, checked by hand on R0L).** A critical pair has to resolve with isomorphic graphs and also with the **same producer for each output element**. Call this a *bordered* match. Only then does the resolution glue into any surrounding context the same way.

In R0L, both branches give the same triangle:
- In branch 1 the final L comes from the third event and the final R from the second.
- In branch 2 the final L comes from the second event and the final R from the third.

The triangles match only through the R↔L mirror. In 0R0L the leading vacancy breaks that mirror, and the graphs come apart, as ChatGPT found. So the right local screen is a bordered critical-pair check. It would have flagged B before the enumeration.

**Why the absorber fails (inference).** In G0PG, the right-hand absorption can happen before or after the left G finishes its hop. P is a single element whose versions form a chain, so its absorptions are forced into an order. Arrivals that are physically independent then end up either chained or parallel depending on the schedule. This is the same root cause as the 35L token pools. **Any single element that collects from independent sources serializes them, and the serial order depends on the schedule.**

**What survives.** A rule passes when its competing matches either never overlap (A) or resolve with a bordered match. No tested rule with overlaps does the latter. 35P and 35Q together suggest a design principle: interactions must not funnel independent inputs into one element's version chain.

## Next steps (in priority order, rough costs)

1. **Systematic rule search instead of hand-picked rules** (~0.1–0.3 core-h). Enumerate all small string rules over 2–3 symbols, with left and right sides of length ≤2 or 3 and at least one overlapping pair of left sides. Keep the rules that conserve a carrier count. Screen each in two stages:
   - (i) a bordered critical-pair check, which is fast and should discard most rules;
   - (ii) the full 35Q finite screen on words up to length 7–8, for the survivors.

   Report how many interacting rules pass, if any, and the simplest ones. Wolfram-model literature reports causally invariant rules with overlaps. If none exist in this class, that is a real negative result for "Model Gc on strings". Treat the literature claims as things to check, not assume.
2. **Run step 2 on rule A now** (≤0.05 core-h). A conserves tokens, is causally invariant on the screened domain, and its event rate depends on density through exclusion. Stages 2a (rate against density), 2b (foliation check) and 2d (probe speed in local ticks) apply as written. 2c can use a density step, which in A forms shocks and rarefaction fans, instead of an absorber. This tests the emergent-clock idea on the one rule we have. Register this as a scope change before running, because the HANDOFF note read the gate as requiring a rule with overlaps.
3. **Absorber in hypergraph form** (design, then ≤0.1 core-h). In Wolfram hypergraph rules, events consume *relations* (hyperedges). The vertices they connect persist. Make the particle a persistent vertex, and make each G carrier an edge attached to it. Absorption then consumes only the carrier edge, perhaps leaving a "mass" edge on P, and never rewrites P itself. Two absorptions at the same P would no longer share a consumed element, so they commute. *Inference:* this removes the serialization found in D. Run the bordered critical-pair check, then the finite screen.
4. **Only after a rule with overlaps passes:** steps 2c (absorber clocks) and 3 (1/r inflow, attraction).

Total for steps 1–3: about 0.2–0.5 core-hours.
