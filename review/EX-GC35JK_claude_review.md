# Claude review of EX-GC35J, 35K and the competing-event design, with next steps

2026-11-25 | Reviewer Claude | Studies designed and run by ChatGPT | EXPLORATORY | PI closures of 35A–35K pending.

Scope: an opinion on what the results mean and what to do next. Richard asked about "35G/J/K". 35G was reviewed in [the 35G–I review](EX-GC35GHI_claude_review.md). Since then ChatGPT has only added the follow-up recipe lines (C≈1.44, the 0.56% extrapolation, and the non-local projector note). Those match that review, and I have nothing new to add on 35G. I did not recompute anything. I read the plans, notes, design doc, HANDOFF and LOG on `main` at 5639e76. Statements marked *inference* are my reading, not something a run showed.

## Summary

- **35J**: the soft-constraint fix I proposed fails the exact momentum identity, and the reason is general. The 34O scalar is exact because it is linear. A nonlinear site potential does not obey the discrete chain rule across a link. So "exact energy and exact local momentum" is probably the wrong gate for nonlinear matter on a fixed lattice. I'd stop chasing it.
- **35K**: a good result. The full Y1 metric (lapse plus spatial links), sourced by matter and run as local proper-time events, bends a packet 2.009× as much as the lapse alone. That is GR's factor of two, now in the event host. The absolute angle is still 10.6% high. That excess is shrinking with spacing and is mostly stencil dispersion.
- **Competing-event design**: sound. It is the right next run, and it is cheap. I'd add predictions made in advance, plus a third graph definition.

## EX-GC35J: soft-constraint Ward gate

**What it means.** The momentum identity needs the link term (∇Vₓ+∇Vᵧ)·(uᵧ−uₓ)/2 to equal Vᵧ−Vₓ. That is the trapezoid rule, and it is exact only for quadratic V. The 34O scalar passed because its mass term is quadratic. ChatGPT's two-point example (remainder 0.75μ/h) shows the failure plainly. The relative defect, 0.713 for every μ, says it is structural, not a tuning problem. The quartic Skyrme candidate also fails its static test (1.003). So moving to a soft constraint removes the 35H reaction term, but it brings in an equivalent chain-rule term.

*Inference:* this is likely a general property rather than a defect of this construction. On a fixed lattice, nonlinear fields keep exact energy, because the time direction stays continuous. Exact local momentum generally needs a continuous translation symmetry, which the lattice does not have. Known "momentum-conserving discretizations" of nonlinear Klein–Gordon (for example Kevrekidis, Physica D 2003, and Dmitriev and co-workers) use discrete-gradient forms like (F(uᵧ)−F(uₓ))/(uᵧ−uₓ). My recollection is that they usually give up the standard energy structure in exchange. That recollection should be checked before anyone relies on it. If it holds, the 35B→35J sequence has been aiming at a target that may not exist for any local nonlinear matter.

**What follows.** There are two sensible paths:
- *Accept an O(h²) source and control the constraints locally.* Energy-only sourcing is already exact (35B/35H). The full-stress mismatch is O(h²) (35D). 35C showed that an O(h²) source error drives a secular constraint drift. 35E/G removed that drift with a non-local projector. Adding local constraint-damping terms would serve the same purpose and could run as a local rule.
- *Get exact conservation from counting instead of calculus.* A rule that moves conserved tokens conserves them exactly by construction, nonlinear or not. That is also the setting of the competing-event design. *Inference:* this is the more natural road to the programme's goal. The PDE matter fields would then serve as continuum checks, not as the fundamental layer.

## EX-GC35K: matter-sourced lapse and spatial links

**What it means.** The metric comes from a static Y1 solve sourced by a 3D Gaussian mass (kernel residual 1e-16). Each site fires at its own proper-time rate. On the finest grid the full bend is 3.135° and the lapse-only bend is 1.560°, a ratio of 2.009. Lattice errors largely cancel in that ratio, which is why it is the robust number. The absolute values run high:
- Full: 10.6% above the matched finite-path estimate of 2.834°, down from 16.0% at h=.16.
- Lapse-only: 8.6% high.

The centred-stencil dispersion correction sec²(kh/2) accounts for most of the excess: it gives 3.00° against the measured 3.135°. That leaves about 4.5% for finite-width, ray-path and periodic-source effects. The spacing ratio of 1.33 should cut a pure h² error by 1.78×. The excess fell only 1.51×, consistent with a mix of h² error and something that does not depend on spacing. So the absolute angle has not been shown to converge yet. It is not contradicted either.

**The most important lesson is the clock interface.** With the 35I linear predictor, the field norm grew 2.93× and the angle came out wrong. With quadratic reconstruction the varied-clock run matched the common-clock control to 0.004°, and the norm ratio went from 0.974 to 0.997 when the step was halved. *Inference:* "lapse = event rate" works only if neighbours that tick at different rates read each other through a good enough reconstruction. In the event graph that reconstruction is part of the rule, and it is approximate. A truly discrete version would need either rates in small rational ratios (as in 35I, 2:3) or an interface rule that is exact by design.

**Limits** (as the note says): the metric is static and external, the packet doesn't source anything, and the equal-time reversal check is bookkeeping only.

## Competing-event design

The design is sound. It is careful to separate graph equivalence from ledger agreement, and that is the key distinction. Two additions:

1. **Register predictions.** The paired-redistribution rule is a chip-firing (Abelian sandpile) system at the level of site counts. For such systems the final configuration and the number of firings at each site are known to be the same for every legal order (the Abelian property). So the ledger comparison is expected to pass by a known theorem. That makes it a check on the code, not a discovery. *Prediction (inference):* the token-identity dependency graph will fail at same-site forks, because the two candidate pairs come from different producers by construction. Writing both predictions down before running makes the outcome informative either way.
2. **Add a third graph definition between ledgers and token identity.** In this "site-history" graph, each event depends on the previous event at its own site and on the hop events that brought it the tokens it needed, with tokens at a site treated as interchangeable (the event uses the earliest arrivals, or arrivals as a multiset). If that graph is order independent while the token-identity graph is not, the result is a positive one. It would mean causal invariance holds once particles of one type are truly indistinguishable, and that connects directly to the QM side of the goal.

## Next steps (in priority order, rough costs)

1. **Run the competing-event enumeration** (well under 0.1 core-h; seconds to minutes). Use short marked lines and small token counts. Readouts are ledgers, firing counts, the site-history graph and the token-identity graph. Register the two predictions above first. Report the smallest inequivalent fork if one exists. This is the crux of the programme's goal and the cheapest step on the list.
2. **Stop the exact-Ward search for knots. Test local constraint damping instead** (~1–2 core-h, could use the PC). Use energy-exact sourcing plus the O(h²) stress, and add damping terms on the Y1 stencil. Do this first on the 35C scalar run, where projected 35E/G is the reference, then on a static knot. Gate: pre-reset-equivalent residual bounded at ~h², and tidal error within a stated factor of 35G. Only if that passes, run the moving knot. *Optional, paper only:* a short literature check on whether exact energy and exact momentum can coexist for local nonlinear lattice fields. If they cannot, record that 35H/35J hit a general barrier.
3. **35K absolute convergence** (optional, ~0.3 core-h). One h=.08 full plus lapse-only pair, to see whether the excess over the dispersion-corrected estimate keeps falling. The factor of two does not need it, so this is low priority.
4. **Later, once 1 and 2 report:** a reciprocal event host, in which the packet (or token matter) sources the Y1 metric and the event rates update from it. Combined with step 1, that becomes "token matter with lapse-set rates", which is the first version where gravity, discreteness and a real causal-invariance test all sit in one rule.

Housekeeping: PI closures for 35A–35K are pending. The v0.5 draft is at revision 5 through 35I, and 35J/K can be folded in when Richard wants.

Total for steps 1–3: about 1.5–2.5 core-hours, all within the exploratory caps.
