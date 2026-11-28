# EX-GC61 — the shift as a travelling field: local events, arrival, settling and the detector without refitting

2026-11-28 | EXPLORATORY | Claude | PI closure pending | [Plan](../../spec/EX-GC61_plan.md) (pushed before any run; amendment A1 before any run) | About 1.25 core-h scored, plus about 0.4 core-h of recompute and post-hoc runs, on the cloud container (Linux, Python 3.13.18, numpy 2.5.3; 4 workers). Scoring: `src/exgc61_summary.py` → `summary.json`; recomputation: `src/exgc61_recompute.py` → `recompute.json`; post-hoc: `src/exgc61_posthoc.py` → `posthoc.json`. This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**The shift no longer needs an instantaneous solve.** A local, order-free event rule carries it causally and settles to the free-space field. Before any wall reflection gets back, a probe in the travelling field reads free-space GR's precession with nothing refitted: the J rate is 1.007 of GR and Σ is 1.001. 59's 2.3% box deficit was an artefact of the static solve.

Six of the ten predictions and all checks passed. The four failures are:
- **P2 and P4:** I timed the arrival as if the source were a point.
- **P8a–d:** my plan left out transforming the probe's preparation. With it put in post hoc, the "readout convention" result holds to 10⁻⁸.
- **P10 (one row):** J's vector response at an off-symmetry site in the source field is 41% off.

λ = 4 and the source law are still put in by hand.

- **Events = the law (P1 PASS).**
  - The event leapfrog matches 53's RK4 wave equation to 8.4 × 10⁻⁵ of max|β| at t = 20.
  - The difference falls by 4.00 at dt/2.
  - Shuffling the events within each class changes nothing beyond 2 × 10⁻¹⁵.
- **Settling (P3 PASS).** At t = 50 the travelling field differs from a near-free-space static solve (n = 193) by 0.22% rms over 8 ≤ |r| ≤ 20. It differs from 59's own zero-box static field by 1.8%.
- **Causality and boxes (P4 part).**
  - Up to t = 70, the n = 97 and n = 129 fields at the probe agree to 3 × 10⁻⁵.
  - Over the probe window [40, 50] the probe rates in the two boxes are identical to 10⁻¹³.
  - So nothing reaches the probe faster than the lattice allows.
- **The probe without refitting (P5, P6 PASS).**

  | | Σ rate | J rate |
  |---|---|---|
  | Axis, d = 16: / free-space GR | 1.0012 | 1.0072 |
  | Axis: / 59's static zero-box rate | 1.025 | 1.031 |
  | Diagonal: / frozen two-coefficient prediction (60's coefficients, travelling field's gradients) | 1.0002 | 0.9947 |
  | Diagonal: / rotation alone | 1.272 | 1.0004 |

  J is the clean rotation readout: rotation alone predicts it to 0.04%. Σ needs 60's shear term.
- **Reflection (P7 PASS).** At t = 90–100, after the reflected switch-on pulse reaches the probe, the on-axis rate drops by 80%.
  - This is the radiation from switching on a spin, reflected by the zero walls. A prescribed spin that turns on is not angular-momentum conserving, so it radiates in a way GR would not allow.
- **Detector add-on.**
  - **P8e PASS.** The 7% residual between the momentum coupling and the exact resync coupling falls with width as predicted from sin k − ½ sin 2k. The ratios are 0.59 at σ = 4 and 0.39 at σ = 5, against predictions of 0.56 and 0.36.
  - **P9 PASS.** At NB = 64, links-only and PS agree to 4.6 × 10⁻¹⁰, which confirms 60's face-wrap diagnosis. k_sh(σ = 5) is unchanged to 10⁻⁵.
  - **P10.** A tilted spin in a general uniform gradient follows the isotropic vector prediction to 0.36% (Σ) and 0.30% (J). At the off-symmetry source site Σ holds (0.9%), but J misses by 41%.

| | What was scored | Result | Verdict |
|---|---|---|---|
| P1 | events vs RK4 < 1e-3, ratio 4 ± 1; shuffle < 1e-12 | 8.4e-5, 4.00; 2e-15 | PASS |
| P2 | t½ − d constant within ±0.75, slope 1 ± 0.03 | axis spread 0.91, slope 0.89; plane 1.38, 0.83 | FAIL |
| P3 | rms vs n = 193 static < 1% | 0.22% | PASS |
| P4 | 97 = 129 to 1e-4 for t ≤ 70; first > 1e-3 in [76, 92] | 3.0e-5 ✓; 74.1 ✗ | FAIL |
| P5 | axis Σ, J / GR = 1 ± 0.015; / 59 in [1.010, 1.035]; n = 129 same to 0.2% | 1.0012, 1.0072; 1.025, 1.031; 1e-13 | PASS |
| P6 | diagonal / frozen prediction within 1%; J / rotation alone within 1% | 1.0002, 0.9947; 1.0004 | PASS |
| P7 | late axis J rate differs > 0.5% | −80% | PASS |
| P8 | matched readout (a) < 1e-3, (b) 4 ± 1, (c) 5%, (d) 1e-3, (e) σ⁻² | 0.063, 2.0, 32%, 2e-3 ✗; 0.59, 0.39 ✓ | FAIL (a–d) |
| P9 | NB = 64 P = PS < 1e-7; k_sh unchanged < 0.1% | 4.6e-10; 1e-5 | PASS |
| P10 | vector rates within 3% (flat and source site, Σ and J) | flat 0.36%, 0.30%; source Σ 0.9%, J 41% | FAIL (source J) |

## Checks

| | | |
|---|---|---|
| C1 | zero source stays zero | 0 — PASS |
| C2 | field energy conserved after the ramp | 4.2e-6 — PASS |
| C3 | 61C probe code with 59's static field reproduces 59's axis rate | identical (0) — PASS |
| C4 | ∇·j / \|j\| of the source current | 6.8e-4 (reported) |

**Recompute.**
- The axis and diagonal probes were rerun with the field and the probe at dt = 0.05, spin read from summed densities and slopes from polyfit. They agree with the main runs to 2 × 10⁻⁶ (axis) and 1.6 × 10⁻⁴ (diagonal Σ).
- Arrival times from a cubic-spline crossing agree to 9 × 10⁻⁴.
- Shared dependencies: the source builder, the field class, the probe Hamiltonian and the frozen coefficients.

## Why the four failures happened (post hoc, not re-scored)

- **P2 and P4: the source is not a point.** At the 10⁻³ level the field first appears 5–10 time units before the light-front time from the centre. That is the packet's tails (σ = 3), which reach to about 10 sites.
  - The half-rise time also contains the switch-on terms of a time-varying dipole (ṁ/r² on the axis, and ṁ/r² plus m̈/r in the plane). These make it earlier the farther out you go, which is why the slope came out below 1.
  - A point-dipole retarded model gets the axis half-rise times right to 0.7–2 time units, but not the plane.
  - The first departure between the boxes, at 74.1, is the image path (80) minus the source's extent.
  - The 10⁻¹³ agreement of the two boxes during the probe window is the clean causality test, and it passes.
  - A clean arrival-speed test needs a compact source or much larger distances.
- **P8a–d: the preparation was not transformed.** The registered D1 started every shear-description run from the free preparation ψ₀ and transformed only the detector. The review asked for one free preparation, which in the shear description is Uψ₀.

  Rerun post hoc from Uψ₀:
  - the transformed detector reads the pure-gauge H_g run as the free one to 1.8 × 10⁻⁹ of the signal;
  - the momentum-coupled shear run reads k_P − k_Hg to 8 × 10⁻⁹;
  - the sum field reads k_rot + (k_P − k_Hg) to 2.5 × 10⁻¹⁰.

  So, with the preparation and the detector both transformed, the shear drift of raw Σ is exactly the coordinate relabelling at first order. What is left is the lattice current difference, which P8e shows scales as σ⁻².

  The registered (a) residual of 6.3%, which is first order in s (doubling s doubles it, hence (b) = 2.0), is the effect of starting from the untransformed state.
- **P10 at the source site: J's orbital part.** L alone does not follow the isotropic two-coefficient form even in the flat field (17% off), and is 113% off at the source site. J = L + Σ/2 inherits that where L's share is large.

  Candidates are the cubic anisotropy of the orbital response, and the field's curvature across the packet (the source field is not uniform). This is open, and J's diagonal-site agreement (0.5%) does not cover off-symmetry sites with tilted spins.

## Deviations

1. **Amendment A1, before any run.** The scored runs use flat weights, so that the static limit is exactly 59's law. The N̄-weighted law (53's) was run as a variant and reported only:
   - it arrives 0.2–0.4 time units later;
   - at t = 50 it is 1.7% from its own zero-box weighted static solve, about the same as the flat field's 1.8% difference from its zero-box solve (the box images that the travelling field has not yet felt).
2. **The 53 comparison in 61A uses an RK4 of 53's shift right-hand side with the matter held fixed,** not `evolve_b` itself, which also evolves matter.
3. **A bug in the job's spin key** ("dn" against "down") crashed the first C97dn job at start-up. It was fixed and rerun; no data were affected.
4. **The Taylor series for U was cut at 6 terms** (60 used 11). The truncation is about 10⁻¹⁴.
5. **The plan's flaws:**
   - P2 and P4 treated the source as a point;
   - D1 did not transform the preparation;
   - P8b expected a second-order residual, but with the untransformed preparation the residual is first order.

## What it means

- **The gravitomagnetic chain is now causal and local end to end:** source current, then local events, then the probe. The shift travels by two-link coincidences, settles behind its front, and reads as free-space linearised GR at the probe. In this model the static box solve was the only thing that made 59 fall short of GR by 2%.
- **The detector question is settled at first order:**
  - with preparation and readout transformed together, the coordinate part of Σ's shear drift is removed exactly, and rotation is kept;
  - the remainder is the lattice momentum-versus-energy-current difference, and it falls as σ⁻²;
  - J is the clean rotation readout in symmetric geometries, but not yet at tilted, off-symmetry sites.
- **A new physical point:** switching on a prescribed spin radiates through the shift, and the radiation dominates the probe once it reflects. GR forbids dipole gravitational radiation because angular momentum is conserved. Our prescribed source breaks that conservation, so this radiation is an artefact of the source, not of the rule.
- **Still put in by hand:** λ = 4, the source j = P + ¼Cᵀs, and the wave law itself. 61 shows these are consistent as local events, not why they hold.
- **Suggested next (not approved):**
  1. Make the source dynamical and conserving (matter evolving in its own shift, with J conserved). That removes the switch-on radiation and tests back-reaction.
  2. Look for what fixes λ = 4 in the event rule (a counting reason, as 37B gave the factor 2).
  3. The off-symmetry J response (P10), if J is to be the benchmark detector.
