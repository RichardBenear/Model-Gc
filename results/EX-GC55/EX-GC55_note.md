# EX-GC55 — one interaction for the link shift: momentum, not energy current; and a clock that reads the shift exactly once it is synchronised

2026-11-28 | EXPLORATORY | Claude | PI CLOSED 2026-11-28 (Richard: "close 55 and 56") | [Plan](../../spec/EX-GC55_plan.md) (pushed before any run) | About 0.75 core-h on the cloud container (Linux, Python 3.13.18, numpy 2.5.3). Scoring: `src/exgc55_summary.py` → `summary.json`; recomputation from raw records: `src/exgc55_recompute.py` → `recompute.json`; post-hoc diagnostics: `src/exgc55_posthoc.py`, `src/exgc55_posthoc2.py`. This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**Of the rules where one operator both sources and feels the shift, only the momentum rule (53's law) keeps the co-moving pair Lorentz-invariant.** 54's energy-current idea fails once it is made reciprocal.

- **Energy bookkeeping** (N = 1, one moving packet, T = 8). The registered test separates the three couplings cleanly:

  | Coupling | Energy drift / interaction energy | |
  |---|---|---|
  | P sources, P felt (53) | 3.7 × 10⁻⁵ (falls 32× at dt/2: time-step error) | conserves |
  | J_E sources, P felt (54's mixed rule) | 4.2 × 10⁻³ | does not conserve |
  | J_E sources, J_E felt (reciprocal, RJE) | 4.0 × 10⁻⁵ | conserves |

  So the Codex concern was right: 54's rule was not a single interaction.
- **Co-moving pair, full dynamics** (n = 49, T = 8, V = 0.2):
  - Momentum rule: R = 0.9986 (53 at n = 65, T = 16: 0.9975).
  - Reciprocal energy current: R = 1.0469, so it fails Lorentz invariance by about 5%.
  - The registered prediction for the difference was +0.030 ± 0.012. The measured +0.048 is right in direction but outside the band, so it is a FAIL.
- **Why the energy current fails.** In the packets' own clock-rate wells J_E is only 0.84 (pair) and 0.88 (single) of its uniform-field value. Its hops and velocity carry the bond weights N_xN_y, while P has none. A reciprocal J_E puts that loss into both the source and the coupling.
- **Reading.** The operator that sources and feels the shift is the lattice momentum P, the same structure as GR's ADM coupling −β·p. 54's counting argument can still supply the strength λ = 4. But in this lattice, the event count should be read as transporting momentum, not lab-time energy.

**The clock (registered fails; post hoc fixed).**
- **Registered result.** The covariant light rule I added (the (∂β)u_t term) changes nothing measurable. The smooth-gauge residuals stay at up to 0.097 of the reference (P1 FAIL), the one-sided pattern is 0.929 of ¼, and the ring wave delays are 0.7–2.3% off readout × 0.935 (P2 FAIL; uniform passes at 1.002). My plan-stage diagnosis was wrong.
- **The residual grows with laps** (one lap, two laps, three laps: −1.88, −3.75, −5.62 × 10⁻³). So it is propagation, not readout. It also grows for narrower pulses. The pulse was launched unsynchronised: the same u, u_t as at β = 0, regardless of the shift around the clock. That mismatch persists every lap.
- **Post hoc fix.** Prepare the pulse in the clock's locally synchronised frame: u → u + χ_loc u_t, with χ_loc the shift summed outward from the clock (Einstein synchronisation near the clock). The sign is fixed on the gauge cases. Then:
  - every smooth-gauge case is below 9 × 10⁻⁵ of the reference (all Λ, phases and amplitudes);
  - the one-sided pattern gives ¼ to 1.2 × 10⁻⁵;
  - uniform β gives wave delay / 2Σβ = 1.0036;
  - the rotating ring gives wave delay / readout = 0.9987, 1.0031 and 1.0047 on the three loops, flipping exactly with the rotation.
- **What changes from 53.** The 0.935 "lattice factor" in 53A and 53D was an artefact of the unsynchronised launch. A synchronised one-clock loop reads ∮β/N⁴ directly, to 0.5%. These numbers are post hoc and are not re-scored.

| | What was scored | Result | Verdict |
|---|---|---|---|
| P1 | smooth gauge < 10⁻³ of (a), all cases | up to 0.097 (post hoc synchronised: 9 × 10⁻⁵) | FAIL |
| P2 | uniform / 0.935 within 0.5%; one-sided ¼ within 2%; ring / 0.935 within 1% | 1.002; 0.929; 1.007–1.023 | PASS; FAIL; FAIL |
| P3 | P/P and RJE drift < 10⁻⁴; mixed > 10⁻³ | 3.7e-5, 4.0e-5; 4.2e-3 | PASS |
| P4 | R_RJE − R_P = +0.030 ± 0.012 | +0.048 | FAIL (direction right) |
| P5 | R_P (n = 49, T = 8) within 1% of 0.9975 | 0.9986 (ratio 1.0012) | PASS |

## Checks

| | | |
|---|---|---|
| C1 | covariant code at β = 0 vs 53 loop | both 0 to 2 × 10⁻¹²; difference 2.5 × 10⁻¹² (different centroid density), not bit-identical as the plan said |
| C2 | P/P drift falls ≥ 10× at dt/2 | 32× — PASS |
| C3 | RJE coupling Hermitian; ⟨H_β⟩ = −Σβs | 1.0 × 10⁻¹⁵; 0 — PASS |

The recompute (independent lstsq on the raw P₁(t)) agrees to 10⁻¹⁵. Fitting t ≥ 2 only gives R_P = 1.0003 and R_RJE = 1.0482, a difference of +0.048.

## Deviations

1. My plan-stage claim that the 53 light rule's missing covariant term caused the clock residual was wrong. Covariant and advected rules give the same residuals to within 2%. The working fix (synchronised preparation) was found after the scored runs and is post hoc.
2. The synchronisation sign was chosen on the gauge cases, then applied unchanged to the controls and the ring.
3. C1's wording ("exactly") was too strong. The two codes use different centroid densities.

## What it means

- **The moving-matter sector is now one consistent interaction:** link shift β, coupling −β·P for matter (and −β·π∂u for light), source P, strength 4, field action with GR's sign.
  - It conserves energy.
  - It keeps a co-moving pair Lorentz-invariant (0.9986).
  - With a synchronised clock, it reads a rotating ring's frame dragging exactly as the loop integral of the shift.
- **What is still put in.** The strength 4 has a counting rationale (54: two-site × Doppler-squared), but the momentum operator is chosen because it is the conjugate of the coupling, not derived from coincidences.
- **Suggested next (not approved):**
  1. Derive −β·P from directed coincidence weights on matter's hops. A hop weighted N_xN_y e^{±iθ} by direction couples to Im(ψ†ψ′) ∝ P at leading order. Test whether the natural θ per unit of event offset gives λ = 4.
  2. Codex's well-depth and speed generalisation, now with the momentum rule: halve κ, add V = 0.15 and a different packet width (about 1 core-h).
  3. Codex's rotation step 3 with the synchronised clock: an absolute ring delay plus a matter probe, under reversal.
