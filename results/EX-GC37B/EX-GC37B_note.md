# EX-GC37B — gravity as event counting: clocks are one-site events, propagation needs two

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC37B_plan.md), with amendment A1, made before any run. PI closure pending. Cloud container, Python 3.13, numba. About 0.1 core-h. This is our construction, not LaViolette's Model G.

**The counting rule works.** Sites fire at a local rate N. A clock needs one site to fire. A signal needs two sites to fire together, and that alone makes signals move at N². That gives the factor of two in light deflection, with no Einstein equation and no separate spatial metric. Tokens slow as N^1.99. Their delay past a mass is 2.00–2.08 times the one-site control, and the wave packet bends 2.04 times as much. A carrier consumed by matter gives N = 1 − GM/r, with GM set by the consumption alone (1.1% between two core sizes). The one registered miss is the absolute bending. It is 9.6% above 4GM/b at the finest spacing, the same lattice-dispersion excess as 35K.

## B1a — tokens (2,000 per row, τ = 0.01)

| Readout | Two-site (coincidence) | One-site (control) | Prediction |
|---|---:|---:|---|
| Speed exponent p, v ∝ N^p, N = 0.5–1 | **1.993** | **1.000** | 2 and 1 (P1: ±0.05) → **P1 PASS** |

Shapiro delay along a straight line at impact b, with N = 1 − 2/r. Delays are in flat hop times.

| b | Two-site delay | One-site delay | Ratio | Same ratio from the rate sum (G3) |
|---:|---:|---:|---:|---:|
| 8 | 42.34 ± 0.19 | 20.32 ± 0.03 | 2.084 ± 0.010 | 2.105 |
| 16 | 33.46 ± 0.16 | 16.51 ± 0.02 | 2.026 ± 0.010 | 2.050 |
| 32 | 26.56 ± 0.14 | 13.29 ± 0.01 | **1.999 ± 0.011** | 2.025 |

- The ratio tends to 2 in the weak field, and every measured ratio is within 1.3% of the rate-sum prediction. **Token P2 (amendment A1) PASS.**
- The delay is proportional to (1 + γ), so this is γ = 1 at the event level.
- **G1, strictly, fails by 0.8%.** The flat two-site hop time is 0.84% longer than one over the stationary coincidence rate (7σ). A token starts its wait just after an arrival, not at a random moment, which gives a correction of order τN. The same offset explains why the measured delays sit 1% below the rate sum.

## B1b — waves (35K packet, mass source and readouts; GM = .095, b = 7.5)

| h, proper step | Rule | Group bend | Phase bend | Finite-path eikonal |
|---|---|---:|---:|---:|
| .16, .04 | flat | 0.000° | 0.000° | 0 |
| .16, .04 | lapse only (one-site) | 1.617° | 1.511° | 1.435° |
| .16, .04 | handshake N_xN_y (two-site) | 3.313° | 3.097° | 2.867° |
| .12, .02 | lapse only | 1.560° | 1.502° | 1.437° |
| .12, .02 | handshake | **3.183°** | 3.065° | 2.873° |

- **Ratio handshake/lapse-only: 2.049 at h = .16 and 2.040 at h = .12.** This part of P2 passes (within 3% of 2). There is no spatial metric in this rule. The spatial half of the bending comes from the two-site count.
- **G2 PASS.** The lapse-only rows reproduce 35K's varied-clock lapse rows to the quoted digits (1.617°, 1.560°).
- **The absolute angle misses: P2 FAIL as registered.** At h = .12 it is 9.6% above 4GM/b = 2.903° (14% at h = .16), against a registered 5%. Against the finite-path eikonal, the excess is 10.8%, the same as 35K's full metric (10.6%). 35K traced that to the centred stencil's group-velocity dispersion, sec²(kh/2), which falls as h². It is a property of the wave stencil, not of the counting rule.

## B2 — firing rate set by a consumed carrier

Steady diffusing carrier, D ∇²c = s(x), with total consumption S = 2, D = 1, on a 128³ periodic lattice; N = c/c∞.

| Core radius | Fitted GM in N = A − GM/r (r = 3 core radii to 24) | Fit residual / (GM/r) | Continuum S/(4πD) |
|---:|---:|---:|---:|
| 1.5 | 0.1601 | 4e-4 | 0.1592 |
| 3.0 | 0.1584 | 9e-4 | 0.1592 |

**P4 PASS.** The field is 1/r outside the core, and GM depends on the consumption, not the core size (1.1% apart, gate 2%). This is Gauss's law for a conserved carrier, as in EX008. The run confirms that the lattice version gives exactly what B1 needs. The carrier-sourced N was not also fed into a token run. B1a's N = 1 − 2/r already has this form.

## Addendum 2026-11-26 (after [ChatGPT's review](../../review/EX-GC37AB_chatgpt_review.md))

- **The wave runs already used a carrier field.** The 35K lapse is the periodic Y1 Poisson solution n = −κρ/(2K²). A conserved carrier with diffusivity D and consumption S = 4πD·GM (fixed in advance, not fitted), with the same Gaussian core, obeys the same lattice equation. On 35K's 160³ grid the carrier's δN and the 35K lapse are **bitwise identical**. So the B1b bend and delay numbers are also the carrier-sourced result, and feeding the carrier in separately adds nothing. What this does *not* test, as the review says, is universality across different kinds of matter.
- **The dispersion-corrected prediction is about 3.05°**, which puts the measured 3.183° about 4.5% high (ChatGPT's figure). Dispersion explains roughly half of the gap.
- **Still untested:** waves propagating under *realised* random coincidences. The wave rule uses the mean face weight N_xN_y.

## Not done

- **P3 is not scored.** Local light speed measured with local clocks and light-defined rulers is 1 by construction: local ticks at N, coordinate speed N², ruler length N. No separate measurement was made.
- **B3 (moving source, preferred frame) is design only, as planned.**

## What it means (inference)

- **This is the first rule in the programme where the GR factor of two comes from a mechanism rather than an equation.**
  - Time runs slow near mass because the substrate fires less often there.
  - Light is slowed twice as much because crossing a link needs two sites to fire together.
  - The rulers and the spatial metric follow from that count. We did not put them in.
- **What it does not yet show:**
  - Why propagation should need a coincidence. Here that is a postulate, though a natural one: a transfer needs both ends.
  - Second-order (PPN β) behaviour.
  - Whether a moving mass keeps this picture without a preferred frame. EX-GC25's objection to a diffusing carrier still stands, and B3 is the test.
- **Connection to Model Gc:** "matter consumes the carrier" becomes "matter lowers the local firing rate." That gives a mechanism for the 1/r potential and the full weak-field metric together.

## Files

- [Token and carrier code](../../src/exgc37b.py)
- [Wave code](../../src/exgc37b_wave.py)
- Results: `tokens_b1a.json`, `tokens_b2.json`, `wave_h16_*.json`, `wave_h12_*.json`
