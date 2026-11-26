# EX-GC41 — event-counting optics on a random network: front bending and shared site clocks

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC41_plan.md), following the [Codex review](../../review/EX-GC40AB_codex_review.md). PI closure pending. Cloud container, Python 3.13, numpy 2.5, scipy, numba. About 0.3 core-h. This is our construction, not LaViolette's Model G.

**The time delay is robust. The bending is not measurable this way, and sharing site clocks changes the two-site delay by about 16%.**
- **41A, delay: PASS.** On a random network with no lengths, the arrival-time delay behind the lens matches the continuum ray calculation to within 2.3%, at b = 40 and 60. This holds for both rules and both lens strengths. This is GR's Shapiro delay (for the two-site rule) recovered absolutely, not only as a ratio.
- **41A, bending: FAIL.** The slope of the arrival front is dominated by roughness of the first-passage front. That roughness does not cancel between paired flat and lens runs. Even 60 graphs leave ±20% errors.
- **41B, shared clocks.** One firing history gives no usable optics, as expected. The median over 21 histories reproduces the deterministic delay for the one-site rule (within 1.4%). For the two-site rule it reads 16% high, so the two-site ÷ one-site ratio comes out 2.2–2.5 rather than 1.9–2.0.

## 41A — front readout (deterministic link times)

Domain: a 300 × 300 square, about 90,000 sites. Plane-wave source at x < −140, receiver strip x = 120 ± 2. The main lens is GM = 2, with sites at r < 5 removed; the weak lens is GM = 0.5.

**Continuum benchmark** (ray family through the same index n = 1/v(N), same lines):
- The ray map is monotone once rays entering r < 5 are dropped (the graph has no sites there).
- Halving the step changes α by 1e-14.

**Pilot (5 graphs):**
- The graph-to-graph sd of α (two-site, main lens, b = 40) was 0.239, against α_c = 0.265.
- The 2% target would need thousands of graphs, so N_g = 60 (the cap).
- This was recorded in the LOG before the scored run, together with the continuum ratios.

**Scored run, 60 graphs.** The table gives lattice ÷ continuum, with the 95% graph-bootstrap interval for α in brackets.

| | b | Delay D/D_c | Deflection α/α_c [95% CI ÷ α_c] |
|---|---:|---:|---|
| two-site, main | 20 / 40 / 60 | 1.026 / **0.999 / 1.009** | 1.18 / 1.12 [0.91, 1.35] / 0.96 [0.78, 1.16] |
| two-site, weak | 20 / 40 / 60 | 1.036 / **1.023 / 1.003** | 0.99 / 1.43 [1.12, 1.77] / 0.92 [0.67, 1.21] |
| one-site, main | 20 / 40 / 60 | 1.041 / **1.006 / 1.005** | 1.06 / 1.34 [1.06, 1.64] / 0.88 [0.67, 1.11] |
| one-site, weak | 20 / 40 / 60 | 1.034 / **1.023 / 1.002** | 0.88 / 1.41 [1.08, 1.78] / 0.97 [0.68, 1.31] |

**Scoring:**
- **P1 (null): FAIL on precision.**
  - The flat front's slope has a 95% interval containing 0 in every window, so there is no bias.
  - But its half-width is 17–22% of α_c(40), against the 1% gate. That is the same front roughness that sinks P3.
- **P2 (delay within 3% at b = 40 and 60): PASS.** All eight values lie between 0.9985 and 1.023.
  - At b = 20 the graph reads 2.6–4.1% high. That is outside the scored range: the field is strongest there, and the front passes closest to the removed core.
- **P3 (deflection within 5%): FAIL.** The values run from 0.88 to 1.43.
  - Most intervals contain 1. The b = 40 windows read high under all four conditions, because they share graphs and front structure.
- **P4 (ratio within 5% of the continuum): FAIL.** Lattice ÷ continuum ratios run from 0.84 to 1.13, with intervals of ±12–15%.
  - The weak-lens clause ("within 3% of 2 at b = 60") already failed at the benchmark, which gives 1.924. That was recorded before the scored run.
- **P5 (isotropy): FAIL, by a hair, on the precision clause.**
  - The spread over 16 sectors is 1.03%, which is 3.2 standard errors (< 3.5, consistent with isotropy).
  - But the standard error is 0.319%, against the required < 0.3%. No anisotropy is detected at the 1% level.

**Checks:**
- **H1: not rerun.** 41A uses the same link weights and Dijkstra that passed 40B's G5.
- **H2: FAIL as registered.** The ray deflection at GM = 0.5, b = 60 is 0.930 (two-site) and 0.964 (one-site) of the straight-line Born integral, against a 2% gate.
  - The difference is ray bowing: the bent ray arrives from larger y, where the field is weaker. It is not an integration error.
  - At GM = 0.01, ray ÷ Born = 1.005 for both rules. This descriptive check was added after H2 was seen.

## 41B — shared site clocks (10 graphs × 21 firing histories)

Domain: 100 × 100 (about 9,900 sites), GM = 1, b = 10 and 20. Flat and lens runs share each site's random stream (common random numbers).

| | b | Deterministic | Median of 21 shared-clock histories | median ÷ det | One history (mean) | Independent edges (40B rule) |
|---|---:|---:|---:|---:|---:|---:|
| two-site delay | 10 / 20 | 14.77 / 11.87 | 17.08 / 13.72 | **1.156 / 1.156** | 17.38 / 14.14 | 14.72 / 12.36 |
| one-site delay | 10 / 20 | 7.87 / 5.90 | 7.98 / 5.87 | **1.014 / 0.995** | 8.15 / 5.89 | 9.57 / 7.71 |
| two ÷ one ratio | 10 / 20 | 1.89 / 2.01 | 2.23 / 2.45 | 1.18 / 1.22 | — | — |

Delays are in equivalent length (delay × flat speed).

- **P6 (one history is not clean optics): confirmed as expected.** The single-history ratio has a relative sd of 3.5–12.5 across graphs and histories. Some histories give a one-site delay near zero or negative.
- **P7 (median front within 5% of deterministic, delays and ratio): FAIL.**
  - It passes for the one-site rule (1.4% and 0.5%), but the two-site median delay is 15.6% high at both b, so the ratio is 18–22% high.
  - So a "collective mode" defined as the median of first passages does not reproduce the mean-field optics for coincidence propagation.
- **P8 (does sharing clocks matter): yes for the one-site rule, marginal for the two-site rule.**
  - Flat speed, shared vs independent edges:
    - two-site: 0.0731 vs 0.0716 (about 1.9 standard errors);
    - one-site: 2.50 vs 3.63 (about 30 standard errors).
  - With one shared clock, every link of a site crosses at that site's single firing. That removes the route diversity that independent draws had invented.

**Checks: H3 PASS, H4 PASS.**
- H3: the mean single-hop time from a random start matches 40B's sampler to 0.4% (N = 0.5 and 1; standard error 0.7%).
- H4: the time-dependent Dijkstra equals brute-force relaxation on a 200-site graph, exactly, for 5 histories and both rules.

## What it means (inference)

- **Delay is the robust signature.** A random network of firing sites, with every link the same, carries the continuum delay of a prescribed lapse to about 2%, for clock-only and coincidence propagation alike. With coincidences that delay is GR's full Shapiro delay. No metric and no lengths were put in.
- **Bending needs a smoother signal than a first-passage front.** On a random graph, the earliest-arrival front is rough on scales of tens of links, and that roughness swamps transverse slopes of 0.05–0.3 rad. A wave or diffusion-smoothed field on the network, or much wider apertures, would be needed.
- **Microscopic clocks matter for coincidences.** With real shared histories, the two-site delay in a lapse well is 16% larger than the mean-field value. The rule then over-weights the lapse's effect on light, ratio about 2.3 instead of 2.
  - Whether this survives a better definition of the collective front (a wave mode, as 37D suggested) is the open question.
  - As it stands, GR's factor 2 is a mean-field property of the rule. The first arrival of real events does not preserve it.

## Files

[41A code](../../src/exgc41a.py) · [41B code](../../src/exgc41b.py) · `continuum.json`, `pilot.json`, `scored.json`, `41b_checks.json`, `41b_run.json`
