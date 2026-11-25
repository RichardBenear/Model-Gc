# Clocks from causal-graph depth: design note

2026-11-25 | Design only | Claude | Follows [the 35P review](../review/EX-GC35P_claude_review.md), step 2. ChatGPT runs it and Claude reviews (PI's assignment, 2026-11-25). No result is claimed here.

## Why

The event hosts in 35I, 35K, 35N and 35O get time dilation by **imposing** each site's next deadline from a prescribed or computed lapse. That works, but the deadlines fix the event order, so nothing branches. Causal invariance then holds trivially. In a causally invariant rewriting rule there are no deadlines. Time has to come out of the causal graph itself. This study asks one question:

> In a rule that passes the step-1 graph tests, does the number of events a site experiences per unit of causal depth depend on its surroundings, and does that dependence behave like a lapse?

There is a physical reason to hope so. In GR, free particles fall toward where clocks run slow. In Model Gc, particles drift toward lower G density. If local event rate falls where G is depleted, the two pictures agree in sign. That is a hypothesis to test, not an assumption.

## Precondition

Run this only on rules that pass the step-1 catalogue: ledger plus unlabelled graph equivalence, with the exhaustive critical-pair check. Keep the 35L paired rule as the failing control. If no interacting rule passes step 1, stop and report. The readouts below are only meaningful when the causal graph does not depend on the schedule.

## Definitions (all computed from the causal graph only)

Each site holds one element. An event rewrites a fixed pattern of adjacent elements, so an event spanning sites i and i+1 belongs to both.

- **Depth d(e):** the length of the longest chain of causal edges from the initial elements to event e. "Layer D" is the set of events with depth D. This is the forward-depth foliation.
- **Reverse depth:** the longest chain from e to the final elements, counted down from the maximum. This is a second, independent foliation.
- **Site ticks:** n_i(D) is the number of events at site i with depth ≤ D.
- **Local rate:** r_i is the slope of n_i(D) against D over a window after the initial transient. **Emergent lapse:** N_i = r_i / r_vac, where r_vac is the rate measured in the matched uniform background.
- **Local signal speed:** mark one token as a probe, using the rule's own tokens, not a new species. Count how many site ticks it takes to cross a region of fixed length, and how many layers. The **proper crossing speed** is length divided by the crossing time in local ticks.

Every readout must be identical across randomized legal schedules. Compare at least 8 schedules bit for bit. Any difference means the rule is not causally invariant on that domain, and that case is invalid.

## Stages

**2a — Readouts on uniform backgrounds.** Use periodic rings, L = 32 and 64, through a fixed total of events. Use the passing rules from step 1, for example (b) two-species exclusion with pass-through. Sweep the token density ρ, and the R/L composition where the rule has two species. Report r(ρ) with schedule-spread and ring-size checks. This gives the calibration curve N(ρ).

**2b — Foliation check.** Recompute r_i with reverse depth. If the two foliations give different rates beyond the finite-window error, the "lapse" depends on the slicing and is not a clock. Report that as a failure of the interpretation.

**2c — Inhomogeneous background.** Make a density step first, using two reservoirs at fixed injection densities. Then use the absorber from step 1, (d): a P element in a G background with reservoir boundaries. Measure the N_i profile and the occupancy profile. Gate: N_i follows the local-density prediction N(ρ_i) from 2a within a stated tolerance (suggested 5%) away from the interfaces.

**2d — Equivalence-style probe test.** Measure the proper crossing speed of a probe in the slow and fast regions of 2c. A lapse in the GR sense slows clocks and signals together, so local physics measured in local ticks looks the same everywhere. Gate: the proper crossing speeds agree within the stated tolerance. If the probe instead slows much more than the clocks, which jamming could do, the effect is a refractive index or a crowding effect, not a lapse. Record which.

## Predictions to register before running

- **P1:** the readouts are schedule-independent for the step-1 passing rules and differ for the 35L control.
- **P2 (tentative):** in exclusion-type rules r(ρ) is not monotonic, roughly like the exclusion current ρ(1−ρ). A token can't move into an occupied site, and nothing happens in an empty one.
- **P3:** following from P2, the sign of the lapse change near an absorber depends on the background density. Below the rate maximum, G depletion near P lowers N, the GR sign. Above it, G depletion raises N. Run one background on each side.
- **P4:** open. Whether 2d passes is the main unknown. The registered entry should state that before running.

## Limits

1D only. A one-dimensional absorber in a ballistic gas gives a step-like G deficit, not 1/r, so this tests the clock mechanism, not the gravitational profile. That is step 3, in 2D or 3D. Random hopping directions would make branches with different end states, so diffusion-like behaviour needs deterministic collision rules. Any 2D extension must repeat the step-1 critical-pair check first. Depth is one natural foliation, and 2b checks only one alternative.

## Cost

Each stage is seconds to minutes on 1D rings of up to 64 sites. In total, under 0.1 core-h within the exploratory cap. The registration entry goes in `explore/LOG.md` and `spec/` before production, as usual.
