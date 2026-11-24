# Claude review of EX-GC35G, 35H and 35I, with next steps

2026-11-24 | Reviewer Claude | Studies designed and run by ChatGPT | EXPLORATORY | PI closures of 35A–35I pending.

Scope: an opinion on what the three results mean and what to do next, as Richard asked. I did not recompute the numbers. I read the plans, notes, HANDOFF and LOG entries on `main` at d44281f, and checked two points in the code (the 34O mass term and the 35G recipe's projector). Statements marked *inference* are my reading, not something a run showed.

## Summary

- **35G**: the projected-Y1 scalar benchmark is second-order convergent. The 2% gate missed only because h=.125 is not quite fine enough. I would close it as a working benchmark and stop refining it. Its real limitation is that the projection step is non-local.
- **35H**: the most informative result of the three. It finds why the knot's stress never matches the gravity constraints exactly: a reaction force from the hard unit-vector constraint. That points to a specific fix (a soft constraint) rather than more stencil tuning.
- **35I**: a clean demonstration that "lapse = local event rate" gives redshift and refraction in an event graph. But the lapse is prescribed and nothing competes for inputs, so the order independence holds by construction. It shows the mechanism works. It is not yet a causal-invariance test.

## EX-GC35G: projected benchmark at h=.125

**What it means.** The reconstructed tidal error fell by 3.98–4.00× at every checkpoint when h was halved, and it levels off at about 2.25% by t=16–20. That is clean h² behaviour with a bounded late error. So the error from h=.25 (≈9%) to h=.125 (≈2.25%) predicts about 0.56% at h=.0625. That figure is an extrapolation, not a run. The 2% gate was set before we knew the size of the h² coefficient. Missing it by 0.25 points does not change the conclusion. The projected Y1 scheme reproduces linearised GR tides to second order over the window tested.

**Caveats.**
1. *The projection does a lot of work.* The free control reaches a 95.5% tidal error and a 9.1e-2 constraint residual. Cumulative velocity corrections sum to 0.84, and the per-step impulse is about 10% of the tidal norm late in the run. The projected solution is a constrained evolution, which is standard in numerical relativity. But the projected solution is not what the Y1 equations produce on their own.
2. *The projector is non-local.* It is a per-Fourier-mode pseudoinverse (recipe §2). A Fourier projection cannot run as a local event rule. So the benchmark is a valid numerical reference, but it cannot be carried into the event-graph route as written. *Inference:* any version that has to be causal and local needs the constraint to be preserved, or damped, by local terms.
3. The late-time levelling is seen over t=16–20 only. That matches the plan's wording ("no asymptotic claim").

**Recommendation.** Close 35G as a second-order benchmark with the recipe as written. Do not run h=.0625.

## EX-GC35H: native staggered Skyrme stress

**What it means.** The 34O-style staggered coupling roughly halves the knot's full-source residual: B1 1.65% → 0.96%, B2 1.31% → 0.68%. The energy-only source stays at roundoff. The useful part is the diagnosis. For a unit-vector field the sigma-stress divergence leaves a link term

(λᵧ − λₓ)(1 − uₓ·uᵧ) / (h·12π²),

where λ is the constraint multiplier. The finite-difference chord between two points on the sphere is not tangent to the sphere. So the constraint force does not drop out of the momentum balance the way it does in the continuum. The rank-one periodic counterexample shows this term alone, exactly. Since 1 − u·u ≈ h²|∇u|²/2, the term is O(h²). *Inference:* this is very likely the h^2.3 residual seen in 35D. If so, the "knot stress incompatibility" thread from 35B onwards has one identified cause, and the cause is the hard constraint rather than the gravity side.

The moving-knot run was correctly not started, as the plan gated it.

**What it suggests.** There are two ways to remove the reaction term:

- **(a) Soft constraint (recommended).** Drop |u|=1 and make u an unconstrained four-component field with a local potential, V = (μ/4)(u·u − 1)². Each component is then a plain scalar. 34O already gives exact energy and momentum identities for a scalar with a site-local mass term (`src/exgc34o.py` lines 26–28, 84). A site-local potential should enter the same way, because it couples through the lapse and volume factors only. The sigma-model limit is μ → ∞. Open points: (i) whether the quartic Skyrme term has an exact staggered identity once u is unconstrained, since 35H only established its first variation; (ii) the degree is protected only while u ≠ 0 everywhere. With a finite potential barrier the knot can unwind through a zero at a site. So μ must be large enough that the barrier is well above the knot energy, and the degree must be monitored.
- **(b) Geodesic links.** Build link gradients from the sphere's arc length (arccos uₓ·uᵧ) so the link difference is tangent. This keeps the topology exact, but the metric variation gets messier, and I don't see an obvious exact identity.

## EX-GC35I: lapse sets the local event rate

**What it means.** This is the step-3 idea from the 35E/F review working in a minimal setting. With equal proper-time steps per event, the N=1.5 side fires 1.5× as often per coordinate time. Static observers see the frequency ratio 0.6671 (target 2/3, 0.07% off at h=.1). The transmitted wave normal turns to 32.85° (Snell 32.79°). All 18 randomized ready-order schedules matched the reference bit for bit. Time dilation now sits in the counting of events, not in a coefficient. That is the right shape for the programme's goal.

**What it does not show yet.**
1. *The lapse is prescribed.* Nothing in the graph produces N(x). Matter has to set it before this can count as gravity.
2. *The order independence is guaranteed by construction.* Each version has one writer and fixed, earlier dependencies. So any legal schedule computes the same DAG, and that would hold for any update formula. This confirms the bookkeeping is right. It does not test causal invariance in the Wolfram sense, where rewrites compete for shared inputs and different choices could give different causal graphs.
3. *The lapse covers only half of light bending.* 34A showed that the full GR deflection needs the spatial metric as well (GR/lapse-only = 2.005). An event-graph light-bending test needs the spatial link weights too.
4. The 0.86% frequency-fit offset in the uniform control is as big as the effect a weak-field test would look for. Any later redshift test at small potential needs a better frequency estimator, or a ratio taken against a matched uniform run.

## Next steps (in priority order, rough costs)

1. **Close out 35G and record it** (0 core-h). Close it as a second-order benchmark: error ≈ C h², C ≈ 1.45 at t=20, with 0.56% at h=.0625 given as an extrapolation only. Add one line to the recipe saying that the projector is a non-local reference tool and not part of a local rule. *Optional, ~0.5 core-h:* replace the projection with local constraint-damping terms (a Z4/Gundlach-style damping of the constraint on the same stencil). The test is whether the pre-reset residual stays bounded at ~h² without any global step. This is the version a local event rule could use. I'd do it only after step 2.

2. **Soft-constraint Skyrme on Y1** (~0.5–2 core-h, gated in three parts):
   - 2a: exact Ward identities for the unconstrained four-component field with potential, sigma and quartic terms, on random fields and the 35H rank-one control (seconds). If the quartic term breaks the identity, stop and report.
   - 2b: relax B1 and B2 at h=.2 for two or three values of μ. Gates: full-source residual at roundoff, degree within the 35D tolerance, no zero of |u|. Report how the knot energy approaches the hard-constraint value as μ grows.
   - 2c: only if 2b passes, run the gated moving-knot coupled test from 35H.

3. **Dynamical lapse in the event graph** (~0.2–1 core-h). Take N from the Y1 lapse equation sourced by a static mass, and keep the 35I firing rule. Add spatial link weights, then measure the deflection of a packet passing the mass against the lapse-only value and 4GM/bc² (as in 34A/34C). That moves 35I from "prescribed geometry" to "geometry set by matter". The lapse can stay static for this step.

4. **Design note on a real causal-invariance test** (paper only, 0 core-h; run later). Propose a rule in which events take conserved inputs from a shared pool, such as G tokens hopping between neighbouring sites, so that two ready events can both claim the same input. Then state what has to be order independent: the final causal graph up to relabelling, or only coarse observables. This is the crux of the programme's goal, and it is still undesigned. It should be settled on paper before anyone writes code.

5. **Housekeeping**: PI closures for 35A–35I are pending. Fold 35G–35I into the v0.5 draft when Richard wants.

Total for steps 1–3: about 1–3.5 core-hours, all within the exploratory caps.
