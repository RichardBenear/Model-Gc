# EX-GC41 — event-counting optics on a random network: bending from the front, and shared site clocks

2026-11-26 | EXPLORATORY | Claude | Plan written on PI request ("yes, go ahead", 2026-11-26 16:38 UTC), after the [ChatGPT/Codex review of 40A/B](../review/EX-GC40AB_codex_review.md) (its tests 1 and 2). Written before any run. **Running needs a separate go-ahead.** Proposed cap 4 core-h. This is our construction, not LaViolette's Model G.

## Why

40B showed two things:
- **What worked.** On a random pair network with no lengths, the counting rule gives v ∝ N² and a Shapiro ratio of 2.00.
- **What did not.** Bending was unresolved, for two reasons:
  - single shortest paths were too noisy;
  - the registered target (a ratio of 2) was wrong for fixed endpoints.

  The stochastic version also drew each link's crossing time independently. Real sites share one firing history across all their links.

The two parts:
- **41A** measures bending from the arrival-time front, against a continuum benchmark computed for the same geometry.
- **41B** replaces independent link draws with shared site clocks. It asks whether a reproducible front survives the microscopic rule.

## Common setup

- **Rules.** Unchanged from 40B and 37B.
  - Two-site rule: a signal crosses a link at a coincidence, with τ = 0.01.
  - One-site control: a signal crosses when the sender fires.
- **Geometry.** Every link costs the same in flat space, whatever its length.
- **Lapse.** N = 1 − 2GM/√(r² + r_s²), with r_s = 3.
  - Main lens: GM = 2, with sites at r < 5 removed (as in 40B).
  - Weak lens: GM = 0.5. Here N ≥ 1/3 everywhere and nothing is removed.
- **Plane-wave source.** Every site with x < x_src starts at t = 0 (a multi-source Dijkstra).
- **Receivers.** All sites in a thin strip x ∈ [x_rec − 2, x_rec + 2].
- **Pairing.** Every lens run is paired with a flat run on the same graph.
  - Delays are expressed as equivalent length: Δ(y) = (T_lens − T_flat) × v_flat, where v_flat is the flat speed on that graph.
- **Front estimators** (fixed before any run):
  - **Delay** D(b): the mean of Δ over receivers with |y − b| ≤ 5.
  - **Deflection** α(b): minus the slope of a least-squares line Δ(y) over |y − b| ≤ 10.
  - The same estimators are applied to the continuum benchmark, at the same source and receiver lines.
- **Uncertainty.** Confidence intervals are graph-level bootstrap intervals (95%, 2,000 resamples).

## 41A — bending from the arrival-time front (deterministic link times)

**Domain.** A square of side 300 (Poisson density 1, r_c = 1.5, about 90,000 sites). Source line x_src = −140, receiver line x_rec = +120. Impact parameters b = 20, 40, 60.

**Continuum benchmark.**
- A plane-wave ray family. About 2,000 parallel rays start at x = −140 with spacing 0.05 in y.
- Each ray is integrated through the index n = 1/v(N(x)), with v the deterministic per-hop speed (37B rate sum for two-site, N for one-site).
- Recorded for each ray: its arrival time and y at x = +120. Interpolating gives T_c(y) and Δ_c(y).
- **Check:** the ray map y₀ → y must stay monotone over the scored windows (no caustics). This is checked first, and any window that fails is reported and not scored.

**Pilot (not scored).**
- Run on 5 graphs (seeds 100–104).
- From its graph-to-graph variance, choose the scored ensemble size N_g. Target: a 95% half-width ≤ 2% of α_c at b = 40 for the two-site rule. Cap: N_g ≤ 60.
- The pilot's variance and N_g are recorded in the LOG before the scored run.

**Scored run.** N_g new graphs (seeds 200 onward).

**Predictions.** Each is scored for the two-site and the one-site rule separately, before any ratio is taken.
- **P1 (null).** The flat front's slope, the same estimator applied to T_flat × v_flat over the same windows, has a 95% interval that contains 0, with half-width ≤ 1% of α_c(40).
- **P2 (delay).** D(b) is within 3% of D_c(b) at b = 40 and 60, for both lenses.
- **P3 (deflection).** α(b) is within 5% of α_c(b) at b = 40 and 60, for both lenses.
  - The weak lens is the clean test. The main lens checks that the agreement holds at u = 2GM/b up to 0.1.
- **P4 (ratio).** The two-site ÷ one-site deflection ratio is within 5% of the continuum ratio for the same geometry.
  - The continuum ratios are computed and written into the LOG before the scored run. 40B's fixed-endpoint values were 1.77 and 1.87; a plane wave differs.
  - The ratio tends to 2 only as GM/b → 0. The weak-lens ratio at b = 60 is predicted to be within 3% of 2.
- **P5 (isotropy, re-test).** Flat speed over 16 sectors from a point source at the centre, on the scored graphs. The max−min spread is < 3.5 × the bootstrap standard error, and the bootstrap standard error must be below 0.3%. This replaces 40B's fixed 1% gate, which was below its own noise floor.

**What a result would mean.**
- If P2–P4 pass, an irregular network of firing sites carries GR's light delay and bending for a prescribed lapse, with no metric put in.
- If the delay passes but the bending fails, the network carries the time delay but not the transverse wavefront shape. That would point to a defect in how fronts coarse-grain on a random graph.

## 41B — shared site clocks

**Domain.** Smaller, because every firing is simulated.
- A square of side 100 (about 10,000 sites).
- x_src = −45, x_rec = +40.
- Lens: GM = 1, r_s = 3 (N ≥ 1/3, no removal), with b = 10 and 20.

**Firing histories.**
- Each site i has one stream of unit-rate Poisson arrival times, from `runner.seed_for`-style SeedSequence streams.
- Its firing times are those arrivals divided by N_i. So a flat run and a lens run share the same random draws (common random numbers), as the review asks.

**Earliest arrival with shared clocks.**
- The crossing time of link (x, y) for a signal at x at time t is the first coincidence after t: one end fires while the other fired within τ before (two-site), or the first firing of x after t (one-site).
- Both are non-decreasing in t, so a time-dependent Dijkstra gives the exact earliest arrival.
- Firing arrays are generated lazily up to the needed horizon.

**Compared on the same graphs.** All four use the same firing histories where they apply.
1. Deterministic mean link times (40B and 41A).
2. Independent-edge first passage (40B's stochastic rule).
3. Shared-clock first passage, one history.
4. A collective front: the median arrival over M = 21 independent firing histories at each receiver. This is the collective-signal readout the review suggested, fixed here before the run.

**Ensemble.** 10 graphs × 21 histories per graph, paired flat/lens.

**Predictions.**
- **P6 (single history).** With shared clocks, the two-site ÷ one-site delay ratio from one history scatters by more than 10% from graph to graph. Expectation: like 40B's independent draws, a single first passage is not clean optics. Recorded as expected.
- **P7 (collective front).** The median-front delay D(b), in equivalent length, is within 5% of the deterministic D(b) on the same graphs, for both rules at b = 10 and 20. The median-front ratio is within 5% of the deterministic ratio.
  - Expectation: uncertain. First passage favours fast routes, and the median of first passages is still a first passage. If P7 fails, the "collective mode" has to be defined some other way, such as a wave field (37D).
- **P8 (does clock sharing matter).** On the same graphs, flat speed v_flat for shared-clock vs independent-edge first passage, single history. Descriptive: the sign and size of the difference. A difference larger than 3 bootstrap SE means correlations among a site's links matter.

## Checks

- **H1.** The deterministic Dijkstra on a chain reproduces the 37B rate sum (as G5 in 40B).
- **H2.** The ray-family benchmark at GM → 0 gives α_c = 0 and D_c = 0. At GM = 0.5 and b = 60, its α_c agrees with the thin-lens Born integral ∫∂_y ln n dx to 2%.
- **H3.** On a 1D chain with shared clocks, the mean single-hop crossing time reproduces 40B's G6 sampler to 1% (same rule, stationary start).
- **H4.** The time-dependent Dijkstra equals brute-force event simulation (propagating every token over every link) on a 200-site graph, for 5 histories, exactly.

## Outputs

`src/exgc41a.py`, `src/exgc41b.py`, `results/EX-GC41/`.
