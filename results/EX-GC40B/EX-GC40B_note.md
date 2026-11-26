# EX-GC40B — event counting on an irregular pair network

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC40_plan.md). PI CLOSED 2026-11-26. Cloud container, Python 3.13, numpy 2.5, scipy, numba. About 0.1 core-h. This is our construction, not LaViolette's Model G.

**37B's counting rule works on a random network, with no lattice and no lengths.**
- Every link costs the same in flat space, whatever its length, so the geometry is pure counting.
- On a 2D Poisson network, signals still move at a speed ∝ N² (exponent 1.996).
- Their Shapiro delay is 2.00 × the one-site (clock-only) delay at b = 60, as GR needs.

**The bending test was not decisive.** Single shortest paths on a random graph are too jagged, and the registered gate was set below the noise floor.

## Results

Network: 70,000 sites in a disk of radius 150 at density 1, links within r_c = 1.5 (mean degree ≈ 7). The lapse is N = 1 − 2GM/√(r² + r_s²), with GM = 2 and r_s = 3. There are 20 graph realisations.

| Test | Deterministic (mean link times) | Stochastic (sampled firings) | Gate | Score |
|---|---|---|---|---|
| P6 isotropy: spread of speed over 16 sectors at d = 100 | 1.42% | 3.8% (two-site), 4.9% (one-site) | < 1% | **FAIL** |
| P7 speed exponent, two-site / one-site | 1.996 / 1.000 | 1.992 / 1.029 | 2 ± 0.05 / 1 ± 0.05 | **PASS** |
| P8 Shapiro ratio, two-site ÷ one-site, at b = 20 / 40 / 60 | 1.898 / 1.984 / **2.004** | 1.69 / 1.82 / 1.37 | 2 within 3% at b = 60 (deterministic) | **PASS** |
| P9 bending ratio at b = 40 / 60 | 2.90 / 1.33 | 1.70 / −2.3 | 2 within 5% | **FAIL** |
| P9 absolute two-site bending ÷ continuum ray, at b = 40 / 60 | 1.30 / 1.01 | 0.77 / 0.88 | within 5% | **FAIL** |

Delays are in flat travel-time units (delay × flat speed), as in 37B.

**P6:**
- A descriptive add-on, run after scoring ([script](../../src/exgc40b_isonoise.py), `isotropy_noise_check.json`), repeated the test on 19 new graphs (one graph was skipped because its source site was isolated). The spread was 1.50%.
- The two sector patterns correlate only weakly (r = 0.42). The expected spread from sampling noise alone is about 1.8%.
- So the 1.4% is consistent with noise. The 1% gate was below what 20 realisations can resolve, which is a design flaw. I find no sign of a preferred direction, but the registered result is a FAIL.

**P8:** the two halves of the realisations give 2.019 and 1.990, so G7 holds for P8. The ratio falls at smaller b, as the continuum does in a strong field (u = 2GM/b = 0.2 at b = 20).

**P9 failed for two reasons, both in the design:**
1. **Noise.**
   - Each turning angle is a line fit to a shortest path, averaged over 7 receivers and 20 graphs.
   - The standard error is 0.045–0.06 rad, against signals of 0.13–0.38 rad.
   - The halves disagree badly (at b = 60 the ratio is 1.00 in one half and 2.23 in the other), so G7 fails for P9.
2. **The target was wrong for this geometry.**
   - With the source and receiver fixed at x = ∓120, the more strongly bent ray bows further from the mass.
   - So even the exact continuum rays give a ratio of 1.62 / 1.77 / 1.87 at b = 20 / 40 / 60, not 2.
   - A ratio of 2 needs a weaker field or a larger box.
   - A wavefront-tilt readout (the gradient of the delay across many receivers) would average far better. That would be a new readout, so it was not used here.

**Stochastic version (descriptive):**
- First passage picks the fastest of many random routes, so travel times are about 3.3× shorter than the mean link times give. The ratios are noisier.
- The Shapiro ratio falls below 2 (1.37–1.82), and its halves disagree.
- Taken literally, a single firing history does not give clean optics. The same thing happened in 37D: light has to be a collective mode, not one token's first arrival.

## Checks

- **G5: PASS.** On a chain, Dijkstra reproduces 37B's rate sum to 7e-16.
- **G6:** the mean sampled crossing time is 0.72%, 0.87% and 1.32% (± 0.22%) above 1/rate at N = 0.5, 0.8 and 1.0. This is the known start-up offset of order τN (37B measured 0.84% at its operating point). Here both clocks start from stationary backward recurrences.
- **G7:** holds for P8, fails for P9 (above).

**Deviation:** with the registered GM = 2 and r_s = 3, N ≤ 0 for r < 2.65. Sites with r < 5 were removed before any run. The shortest paths never came closer than r = 12 (stochastic) or r = 21 (deterministic), so the removal did not touch them.

**Code fixes during the run, recorded:**
- The first deterministic lens run crashed in its summary step: the continuum ray shooting bracket reached into the core, where N < 0. The bracket was restricted to the primary image (initial angle 0 to 0.35), and the run was repeated.
- At the same time, the delay ratio was put in flat travel-time units, as 37B and the plan's "ratio of delays" intend. Raw delays differ by the 50× ratio of flat speeds.

## What it means (inference)

- **Event counting does not need a lattice.** On an irregular network with no lengths at all, "a clock is one firing, a signal is a coincidence" gives v ∝ N² and GR's factor-2 Shapiro delay. That is 37B's result moved onto the kind of substrate Chain Theory assumes, and the strongest version of it so far.
- **Bending is unresolved, not refuted.** The deterministic absolute bending at b = 60 is 1.005 × the continuum ray, but the noise makes that a coincidence until a wavefront-based readout confirms it.
- **The stochastic optics are poor.** As in 37D, the signal has to be an average over many firings, a wave, not a single first passage.

## Files

[Code](../../src/exgc40b.py) · `checks.json`, `isotropy_{det,stoch}.json`, `exponent_{det,stoch}.json`, `lens_{det,stoch}.json`, `lens_rows_{det,stoch}.json`, `isotropy_noise_check.json`
