# Claude review of EX-GC35R, 35S and 35T, with next steps

2026-11-25 | Reviewer Claude | Studies designed and run by ChatGPT | EXPLORATORY | PI closures of 35R–T pending.

Scope: what the three results mean and what to do next. I read the plans, notes, HANDOFF and LOG on `main` at 3b145f4. I did not recheck the enumerations. One check was recomputed (35S rate law, below). Statements marked *inference* are my reading, not something a run showed.

## Summary

- **35R:** in the short-string slice, no rule with competing matches has a local certificate that glues into any surrounding history (zero strict-border survivors). Four pass the finite plain-graph screen. I think the strict screen is slightly too strict, and that the ternary survivor is the 35Q absorber failure in disguise (see below).
- **35S:** the density-rate law is real and reproducible, but it is the known rule-184 law, and the density step does not behave like a lapse. Exclusion gives crowding, not time dilation.
- **35T:** a persistent particle vertex with per-carrier consumed relations absorbs in parallel and passes. A shared counter fails. This confirms the design principle from the 35Q review. It passes because nothing competes, not because a competition resolves.
- **Pattern across 35L–35T:** everything that passes has no overlapping matches (independent hops, exclusion, persistent-vertex absorption, the 35F Y1 event graph). Everything with a real competition fails. That is now the main fact to plan around.

## 35R — systematic string search

The search did what the 35Q review asked, with a cleaner two-border design than I proposed. The headline negative is solid: 0 of 734 sets pass with both input and output ports anchored.

**The strict border may over-reject (inference, worked by hand).** Take the ternary survivor 01→1, 12→1 on 012. In one branch the input 0 feeds the first event; in the other it feeds the second. Anchoring that port is the right test only if a 0 can be *produced* by an earlier event in some larger context. Here neither rule ever produces a 0 or a 2, so those ports are always initial elements at depth zero. The mismatch can then never show in any reachable history, which is consistent with the plain-graph pass at length 7. A correct local certificate should anchor only ports whose symbol some rule can produce. The binary survivor (001→01 produces a 0) does not get this excuse.

**The ternary survivor is the 35Q absorber (inference).** 01→1 and 12→1 is a carrier 1 eating 0s from the left and 2s from the right. That is GP→P, PG→P with P = 1. It passes only because its food is never delivered by transport. Add any hop rule that produces 0s or 2s and the anchored mismatch becomes reachable, which is the G0PG failure in 35Q.

So I read 35R as: on strings, every competing-match rule in this slice either fails or passes for a reason that disappears once carriers move. I would stop string enumeration after one cheap fix (step 1).

## 35S — exclusion as an emergent clock

**The rate law is rule 184 (recomputed, separate path).** For a confluent rule, the causal depth of an event equals its time under the greedy parallel update. The greedy parallel update of 10→01 is the cellular automaton rule 184, whose flux per site per step is min(ρ, 1−ρ). Each move counts at two sites, so the site rate is 2 min(ρ, 1−ρ). A 30-line rule-184 simulation on a 64-ring reproduces the 35S table exactly at all five densities (0.25, 0.5, 1.0, 0.5, 0.25). So the ρ(1−ρ) prediction was the random-sequential law, and depth layering selects the parallel one. That is a useful lesson in itself: **causal-depth time is maximal-parallel time.**

**Why it is not a lapse.** A lapse must slow every process at a place by the same factor. Here (a) the forward and reverse depth rates in the step region differ by about 21%, and (b) the probe's calibrated speed differs by 42% between plateaus. The probe moves on its own events, so its own clock is degenerate (one site per tick, as ChatGPT noted). A density step in rule 184 is also a moving shock or fan, not a static region, so forward/reverse asymmetry is expected there. My conclusion is to drop exclusion as the clock mechanism. The 35I result (prescribed lapse = event rate) still stands; what is missing is a rule that produces it.

## 35T — persistent-vertex absorber

This works as I hoped, and the counter control fails as predicted. Two limits matter for what comes next:

- Routes are prescribed and carriers do not interact. The positive family is causally invariant for the same reason as 35P rule A: no two events ever compete.
- An absorber fed by straight inflowing carriers is the EX-GC16 / Le Sage picture. EX-GC17 and EX-GC25 already showed that gravity from absorbed carriers defines a rest frame and gives drag. So 35T should not be developed into an inflow-and-attraction model on its own. *Inference:* its value is as a building block for a particle that consumes without serializing.

## The strategic point

The program needs causal invariance for the gravity link and branching for the QM link (multiway / branchial structure). Rules with no overlapping matches give causal invariance for free, but their multiway system has no real branching, so they cannot supply the QM side. The hard target is a rule with overlaps whose branches rejoin with the same causal graph. Strings have not produced one. Hypergraph rules are the remaining place to look, and the Wolfram-model literature reports causally invariant rules with overlaps there. Those reports should be checked, not assumed.

## Next steps (priority order, rough costs)

1. **Close the string line properly** (≤0.05 core-h). Rerun 35R with a "producible-port" border that anchors only ports some rule can produce. Report whether the ternary and binary survivors become certified. Then add one transport rule (e.g. 0 and 2 hopping toward the 1) and confirm the ternary survivor fails, as predicted above. Register these predictions first. After that, stop enumerating strings.
2. **Small hypergraph critical-pair search** (0.3–1 core-h). Wolfram-style rules on binary relations with 1–2 relations on the left side (signatures up to 2₂→3₂ or so), one relation type conserved as the carrier. Use the 35R two-stage screen (bordered critical pairs, then finite multiway screen), with the producible-port fix. Include the literature's reported causally invariant rules as positive controls. This is the decisive test of "overlaps that resolve" in the representation where 35T showed persistence helps.
3. **Recoiling absorber in persistent-vertex form** (≤0.1 core-h). Let each absorbed carrier deposit a direction relation on P, and let P step one vertex when the net deposit reaches a threshold. Moving P changes which carriers it can absorb, so this is the first case where the physics forces an overlap. Prediction open. This is the discrete EX-GC16 recoil with identities kept explicit.
4. **Clock rule, only if Richard wants it now** (≤0.05 core-h). A probe whose hops need a background event at its site, in the dilute regime where rate ∝ ρ. That ties the probe to the background clock, which 35S lacked. I rank this below 2 and 3.
5. **Keep paused:** 1/r inflow and two-absorber attraction (drag already known from EX-GC17/25), and continuum Y1/knot work.

Total for steps 1–3: about 0.5–1.2 core-hours, cloud is fine. Step 4 adds about 0.05.
