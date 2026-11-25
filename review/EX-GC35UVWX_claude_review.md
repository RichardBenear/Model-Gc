# Claude review of EX-GC35U, 35V, 35W and 35X, with where to go next

2026-11-25 | Reviewer Claude | Studies designed and run by ChatGPT | EXPLORATORY | PI closures of 35U–X pending (35R–T PI CLOSED 2026-11-25).

Scope: what the four results mean, and an assessment of the discrete causal-invariance line (35L–35X) as a whole. I read the plans, notes, code and HANDOFF on `main` at 364de88. I did not recheck the enumerations. One check was run (35V controls, below; [script](EX-GC35V_control_check_claude.py), seconds). Statements marked *inference* are my reading, not something a run showed.

## Summary

- **35U:** both predictions from my last review held. The ternary string pass is real but only because its food is never produced. Add transport and it fails. The string line is finished.
- **35V:** no typed hypergraph rule found an anchored join. **But the published causally invariant control also failed under this test.** That means the test, as built, cannot yet tell a good rule from a bad one, so the 35V negative is not decisive. See below.
- **35W:** recoil breaks confluence outright. Different orders end in different *ledgers*, not only different graphs. This is robust to any graph convention.
- **35X:** the gated probe rides along with moving carriers, so it does not read the fixed-site rate. The crowding-clock idea is finished.
- **Overall (inference):** fourteen studies show one consistent pattern. Rules pass when nothing competes, and fail when particles both move and interact on a fixed lattice. The one open question that could change this is whether our test is stricter than the literature's definition. That is cheap to settle and should come first.

## 35U — producible borders on strings

Clean and as predicted. The witness `0.12` with the transport rule `0.→.0` is exactly the case I described: once a 0 can be delivered by an earlier event, which input feeds the first absorption becomes visible in the graph. I agree with ChatGPT that the short-string search is complete.

## 35V — typed hypergraph critical pairs

The catalogue is well built, but the control result matters more than the catalogue result. The overlapping rule `xy,zy→xw,yw,zw`, which the Wolfram technical introduction presents as causally invariant, was "unresolved" at three steps, with and without input anchors.

**Check (run, separate path through the 35V code).** I reran the three shared-edge forks of that control on the same star seed under four conventions:

| Seed vertices | Input anchors | Result |
|---|---|---|
| labelled (35V) | on | unresolved at 3 |
| labelled (35V) | off | unresolved at 3 |
| unlabelled | on | unresolved at 3 |
| unlabelled | off | **join at step 1** |

The Wolfram definition compares states and causal graphs up to isomorphism, with no vertex names. Under that convention the control passes at once. But the join at step 1 comes from the seed's symmetry: the three leaves are interchangeable. So this is a weak positive, not a proof that the literature convention is right. It does show that **35V's negatives were scored under a convention stricter than the published one.** Which convention is physically right is a real question:

- Labelled boundaries and anchors are what you need to glue a local resolution into *any* surrounding history. That is the critical-pair lemma condition.
- Unlabelled comparison is what Wolfram's causal invariance asks for, from given initial states.

*Inference:* a rule can pass the second without the first, if its forks resolve only through a symmetry of the state. Whether that is enough for the physics (Lorentz-type invariance of the causal graph) is not settled by anything we have run.

## 35W — recoil

This is the most important of the four. When P moves, a carrier that was eligible to be absorbed stops being eligible. So the order of "absorb" versus "move" changes the final counts. A higher threshold does not help. ChatGPT names it correctly: it is a read/write race on P's position, even though P itself is never consumed.

*Inference:* this race is generic. Any particle that moves and interacts through its position will have it, on any fixed lattice. Making eligibility "sticky" does not remove it; it moves the race to the moment a carrier arrives. So the persistent-vertex trick (35T) solves sharing but cannot solve motion.

## 35X — gated probe

In the dilute regime every exclusion carrier moves one site per layer (free flow in rule 184). A probe that hops on background events at its site therefore follows a carrier, and its speed is about one regardless of density (0.95 versus 0.98). One site per gated tick holds by construction, as ChatGPT says. I agree the clock-from-crowding line should stop. The 35I/35K result, where the lapse is the event rate, still stands. What is missing is a rule that *produces* that lapse. Exclusion does not.

## Where the programme stands (inference)

The programme needs causal invariance for the gravity link and branching for the QM link. From 35L–35X:

1. **Passing rules don't branch.** Independent hops, exclusion, persistent-vertex absorption, and the 35F Y1 event graph all pass because no two events ever compete.
2. **Interacting, moving particles branch and don't rejoin.** Pools (35L/P), whole-element absorbers (35Q), transport-fed absorbers (35U) and recoil (35W) all fail. 35W fails even on terminal counts.
3. **Wolfram's causal invariance is a property of special rules in which space itself is the graph.** There, particles would be persistent structures in the graph, not vertices with a position on a fixed lattice. Our tests have put Model Gc particles *on* a lattice, and 35W suggests that is the root of the race.

So there are two honest paths:

- **(A) Stay on the gravity side with deterministic local rules.** The Y1 event graph is schedule-independent by construction and already gives time dilation from event rate and the full-to-lapse bending ratio 2.0 (35I, 35K, 35N). That is a working, causally consistent gravity engine. Its cost: no multiway branching, so no QM link from it yet.
- **(B) Go to structure-free hypergraph rules,** where space is the graph. Ask whether any rule that passes a calibrated causal-invariance test supports a persistent, localized structure that could carry the Gc particle. This is the route that matches the programme's goal. It is also much longer and has no guarantee.

## Next steps (priority order, rough costs)

1. **Calibrate the test** (≤0.2 core-h, cloud). Score the 35V catalogue and the two published controls under both conventions: (i) unlabelled states and graphs, the literature definition, and (ii) the anchored local certificate. Use at least one *asymmetric* seed for the overlapping control, such as a star with one extra edge on a leaf, so that a pass cannot come only from symmetry. Also add 2–3 further rules the Wolfram material lists as causally invariant. Register predictions first. If no published rule passes (i) on asymmetric seeds, our code or the published claim needs a closer look before anything else.
2. **Take stock in the white paper** (no compute). Add a v0.6 section on 35L–35X: the design rules (no shared consumed element; persistent vertices; read/write race of a moving particle; causal depth equals maximal-parallel time; exclusion is not a lapse) and the open question from step 1. This keeps a long run of negative results from getting lost.
3. **Decision for Richard after step 1:**
   - If calibrated rules with overlaps exist: **path B**, a small search for persistent localized structures in those rules (1–3 core-h), before any attempt to attach Gc physics.
   - If not, or if you prefer the surer route now: **path A**, resume the paused Y1 line with matter, i.e. the knot on Y1 with its own event-rate lapse (1–3 core-h, PC). Treat "causal invariance" for now as deterministic local asynchronous updating, and state plainly that the QM link is deferred.
4. **Stop:** string searches, crowding clocks, and absorbers with fixed-lattice positions (35W shows they can't be confluent). Keep paused: inflow and attraction.

Steps 1–2 together cost about 0.2 core-hours. I recommend step 1 first because it decides whether 35V's negative stands.
