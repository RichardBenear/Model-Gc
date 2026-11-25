# Claude review of EX-GC35L, 35M, 35N and 35O, with next steps

2026-11-25 | Reviewer Claude | Studies designed and run by ChatGPT | EXPLORATORY | PI closures of 35L–35O pending (35A–35K PI CLOSED 2026-11-24).

Scope: an opinion on what the four results mean and what to do next. I read the plans, notes, raw-result keys, the literature note `docs/EX-GC35_local_matter_conservation.md` and HANDOFF on `main` at 5b84ea5. I re-scored 35L with two coarser graph comparisons as a post hoc check (below, and `src/exgc35l_posthoc.py`, `review/EX-GC35L_posthoc_claude.json`, 34 CPU-s). Otherwise I recomputed nothing. Statements marked *inference* are my reading, not something a run showed.

## Summary

The four results split cleanly into two sides.

- **Gravity side: works in the event host.** Local clocks set by a matter-sourced lapse give GR's factor of two in bending (35K, 35N: 2.009 and 2.007). In a 1D toy, two-way feedback runs as local events with schedule-independent output (35O). The leftover 4.4% absolute gap does not shrink with spacing, so it is a comparator effect, not the event mechanism.
- **Matter side: the real obstacle is now clear.** Nonlinear field matter can't easily get exact conservation on a fixed lattice (35J, then 35H/35M). Token matter gets exact conservation for free, but the simple carrier rule is **not causally invariant**. Ledgers agree, causal graphs don't (35L, 35O). My post hoc check shows this is not an artefact of labelling: it survives dropping site labels altogether.

*Inference:* the programme's critical path is now a carrier rule that is causally invariant at the graph level. Continuum knot work on Y1 can pause.

## EX-GC35L: competing consumptions

**What it means.** Both registered predictions held. Ledgers and per-site firing counts agree in all 494 configurations, as the Abelian (chip-firing) property guarantees. The token-producer graph fails in 294 configurations, starting from (2,2). My suggested "site-history" graph, with interchangeable FIFO arrivals, also fails, in 176 configurations. So treating tokens of one kind as identical does not rescue causal invariance for this rule.

**Post hoc check (mine, not a registered readout).** The (0,5,0) witness is left–right symmetric, and 35L compared graphs with fixed site labels. A Wolfram-style comparison usually treats events as unlabelled, so I re-scored the same enumeration with the same rule code, imported unchanged:

| Graph | Labelled (35L) | Labelled up to mirror | Unlabelled (node and edge kinds only) |
|---|---:|---:|---:|
| Token producer | 294 | 294 | 294 |
| FIFO site history | 176 | 173 | 167 |

The labelled counts reproduce 35L exactly. Allowing the mirror symmetry or dropping site labels rescues only 3 and 9 configurations. The first unlabelled FIFO failure is (0,4,2), which is not symmetric. So the failure is genuine causal-structure dependence on order, not a choice of labels.

**Why it happens (inference).** Each event draws from a shared site pool. Whether an arrival lands before or after a local firing changes which earlier event the firing depends on. The final counts converge because the rule is Abelian, but the dependency structure does not. Confluence of states is weaker than causal invariance, and this rule has only the first. In Wolfram's picture that is the expected split: branching with matching final states is the "multiway" part, but relativity-like behaviour needs the causal graphs themselves to agree.

## EX-GC35M: local constraint damping

**What it means.** This particular damping (trace correction from the Hamiltonian residual plus an adjoint correction from the momentum residual, α=β=.12) does not do what the non-local projector did. The fine residual stays at 0.045 and the tidal error reaches 27.6% at h=.25, against 8.94% projected at the same spacing.

The held-static knot audit is the more important half. With the full (O(h²)-incompatible) stress the residual grows to 0.15–0.18. With energy-only sourcing it stays at roundoff. *Inference:* this is expected in principle. If a source is not exactly conserved, no metric satisfies all the constraints, so damping can only spread the error around, not remove it. Damping protects against numerical drift. It cannot repair a source that breaks the Bianchi identity.

So 35M tested one ad hoc form. A principled Z4-type formulation (extra constraint variables with their own damped evolution) was not tested. Even if it worked, it could not fix the knot's source incompatibility. I would not spend more on this now.

## EX-GC35N: lens at h=.08

**What it means.** The full/lapse ratio is stable: 2.009 at h=.12 and 2.007 at h=.08. The gap left after the leading dispersion correction also stays put: 4.37% then 4.38% for the full metric, and 2.44% then 2.53% for lapse only. A gap that does not change with spacing is not a lattice error. 35K had already shown that the varied-clock run matches the common-clock run to 0.004°, so it is not the event mechanism either. *Inference:* the gap lies between a finite-width wave beam and the straight-ray, waist-averaged comparator. A continuum wave solve of the same beam would settle it, but the factor-two result doesn't need it. I'd close the lens line here.

## EX-GC35O: reciprocal local clocks

**What it means.** In a 1D reduction, packet energy sources a retarded lapse wave, and the lapse sets each site's next deadline. The loop closes with every legal schedule byte-identical (8/8 at both steps). The coupling is a real effect: turning feedback off changes the packet by 0.447, against a step-halving change of 0.062. That's a working local, two-way gravity–matter clock loop, a first for the programme.

The caveats are as stated in the note: 1D, a lapse equation imported rather than derived from Y1, first-order integration, and no convergence gate. Schedule independence again comes from the deadline protocol, which works like conservative parallel discrete-event simulation. It is structural, not a physics result.

The token part reproduces 35L inside the clock loop. Counts, fields and hop times agree, but the carrier graphs differ. The lapse does feed into token timing (last hop at 1.368 with feedback, 1.360 without), so the hookup itself works.

## Next steps (in priority order, rough costs)

1. **Search for a causally invariant carrier rule** (≤0.1 core-h; minutes). Use the 35L domain and code. Score each candidate on ledgers, the labelled graph and the unlabelled graph, with a prediction for each written down before running. Candidates:
   - (a) **No-overlap control:** single tokens hop in a fixed direction carried by the token. It should pass trivially and checks the scorer.
   - (b) **35L paired rule:** should fail. This is the negative control.
   - (c) **Whole-site firing:** a site with ≥2 tokens consumes all of them and sends a fixed split, so there is no subset choice.
   - (d) **Per-link FIFO queues consumed whole** (an Abelian-network processor).
   - (e) **Collision rule:** two opposite-moving tokens meeting at a site.

   Report which pass. Also report, for any that pass, whether they still conserve tokens and allow a lapse-set rate. A short literature check (Gorard's papers on causal invariance and confluence in Wolfram models, Bond–Levine on Abelian networks) can go in the same note. If nothing in the catalogue passes, write down the structural reason. That would be informative too.
2. **Pause the knot-on-Y1 continuum line** (0 core-h). Record its status in HANDOFF: energy-only sourcing is exact; full stress is O(h²) incompatible; local damping cannot fix an unconserved source; exact nonlinear identities are blocked for the tested couplings (35H, 35J, literature note). Reopen it only if step 1 points to token matter that needs a continuum check.
3. **Close the lens line** (0 core-h). Record that the factor 2.007–2.009 is stable and that the 4.4% gap is spacing-independent, so it is a comparator effect. Optional and low priority: a continuum wave solve of the same beam (~0.2 core-h).
4. **If step 1 finds a passing rule:** put it into the 35O clock loop in place of the paired rule (~0.1 core-h). Check that the causal graphs stay equivalent when the lapse sets the rates. That is the first combined test of gravity-set clocks with causally invariant matter. Only after that would 2D, or deriving the lapse wave from Y1, be worth the cost.
5. **Housekeeping:** PI closures for 35L–35O. I'd suggest freezing white paper v0.5 at this point. The gravity-in-event-host results (35F, 35I, 35K, 35N, 35O) and the matter obstacle (35H, 35J, 35L, 35M) make a natural chapter break.

Total for steps 1, 3 (optional) and 4: about 0.2–0.4 core-hours.
