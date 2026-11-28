# EX-GC60 — what the spin detector responds to: rotation, shear and a gauge test, with the coupling unchanged

2026-11-28 | EXPLORATORY | Claude | PI CLOSED 2026-11-28 | [Plan](../../spec/EX-GC60_plan.md) (pushed before any run) | About 0.57 core-h on the cloud container (Linux, Python 3.13.18, numpy 2.5.3; 4 workers). Scoring: `src/exgc60_summary.py` → `summary.json`; recomputation: `src/exgc60_recompute.py` → `recompute.json`. This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**The spin detector is now fully characterised, and 59's off-axis "failures" are a readout convention, not a fault in the coupling.** [Correction 2026-11-28, Codex review: "fully characterised" and "entirely a readout convention" overstate what was shown. See the Addendum.] The probe reads rotation exactly. Its spin vector Σ also drifts in a sheared shift, and the gauge test shows that the same drift appears when the shift is nothing but a relabelling of time. So it comes from reading Σ on a slice of constant coordinate time, not from the coupling. The plan's diagnosis (a gauge defect of the momentum coupling) was wrong: P6b failed. The 61 coupling change is not needed. Six of the eight predictions and all checks passed.

- **Rotation is read exactly (P1 PASS).** Spin and total angular momentum turn at:
  - −0.9976 of the rate at σ = 3;
  - −0.9991 at σ = 4;
  - −0.9996 at σ = 5.

  This holds in both orientations. With links only, the spin turns at −0.101.
- **The shear response fades with packet width (P2, mostly).** The spin's response per unit shear is 0.0927, 0.0634 and 0.0456 at σ = 3, 4, 5. That is 0.92, 0.96 and 0.99 of the packet's orbital share 1 − ⟨Σ⟩, a spread of 7.6% against a 10% gate. P2 still fails on its last clause: links-only and PS differ by 7 × 10⁻⁶ at σ = 5 (gate 10⁻⁶). The difference comes from the curl operator's wrap-around at the box faces, which a wide packet's tail reaches. In the interior the curl of the shear field is zero.
- **Total J is nearly immune (P3 PASS).** J responds to shear at only 2.1%, 1.3% and 0.9% of the spin's response.
- **The fields are linear and additive (P4 PASS).** The sum field gives rotation plus shear to 10⁻¹⁰, and the yz orientation matches xy to 2.5 × 10⁻⁸ (P5 PASS). [Correction 2026-11-28: this first said 10⁻¹⁰ for both; the scorecard's largest orientation difference is 2.45 × 10⁻⁸.]
- **Gauge test (P6a PASS, P6b FAIL).**
  - Under the exact first-order resynchronisation coupling H_g = −(i/2)[χ, H²], the evolution matches the gauge-transformed free evolution to 3 × 10⁻⁹. The test is sound.
  - In that pure-gauge evolution, the spin angle still drifts at 0.0864 per unit, 0.93 of the momentum coupling's 0.0927.
  - So the shear response is almost entirely present when nothing physical happens. It is the change of ⟨Σ⟩ under a position-dependent time relabelling (U = e^{−i½{χ,H}}) of a spreading packet, a property of reading Σ at equal coordinate time.
  - The remaining 7% is the lattice difference between momentum and energy current (sin k against ½ sin 2k) that 56 found.
- **A frozen calibration predicts a new position (P7 PASS).** At the diagonal site (12, 0, 12) in 59's source field, the σ = 3 coefficients and the density-weighted lattice gradients predicted:
  - the spin rate to 1.0002;
  - the total-J rate to 0.994.

  The same formula gives 59's axis rates to 1.002 (Σ) and 1.008 (J), and 59's in-plane Σ rate to 0.9984 and J to 1.0081 (59 measured −1.2397 × 10⁻⁵ and −1.6881 × 10⁻⁵).
- **Convergence (P8 PASS).** 59's in-plane rate changes by 2.6 × 10⁻⁴ at dt/2.

| | What was scored | Result | Verdict |
|---|---|---|---|
| P1 | rotation Σ, J = −1 within 1% (3 widths, 2 orientations); links-only Σ −0.10 ± 0.02 | −0.9976 to −0.9996; −0.101 | PASS |
| P2 | k_sh falls with σ; k_sh/f equal within 10%; P = PS within 1e-6 | falls; spread 7.6%; 7.3e-6 at σ = 5 | FAIL (last clause) |
| P3 | \|k_sh(J)\| < 0.1 k_sh(Σ) | 0.021, 0.013, 0.009 | PASS |
| P4 | sum = rot + shear within 0.5% | 1e-10 | PASS |
| P5 | yz = xy within 1% | 1e-10 | PASS |
| P6 | (a) H_g = gauge-transformed free within 1e-3; (b) H_g shear rate < 10% of P | 2.7e-9; 0.93 | (a) PASS, (b) FAIL |
| P7 | diagonal site Σ and J within 3% of frozen prediction | 1.0002, 0.994 | PASS |
| P8 | 59 in-plane at dt/2 within 0.5% | 2.6e-4 | PASS |

## Checks

| | | |
|---|---|---|
| C1 | two β = 0 runs identical | 0 — PASS |
| C2 | shear at 2s gives twice the rate | 2.0000000 — PASS |
| C3 | 60's code reproduces 59's axis rate | identical — PASS |

**Recompute.** The σ = 3 rotation and shear rates were rerun at dt = 0.05, with the spin read from summed densities and slopes fitted by polyfit; see `recompute.json`. They agree with the main runs to 2.4 × 10⁻⁷ (rotation) and 7.0 × 10⁻⁴ (shear, a dt effect at the second-order shear signal). The curl of the shear field is zero in the interior and non-zero only at the faces, from the wrap-around (interior 2 × 10⁻²⁰, faces 2.4 × 10⁻⁴ = 24 s). Shared dependencies: the probe and field builders and the Hamiltonian.

## Deviations

1. The summary's C2 line first compared the 2s rate with 1 instead of 2. I corrected that scoring bug before reporting; the data were unchanged.
2. P2's 10⁻⁶ clause did not allow for the curl operator's periodic wrap at the probe-box faces.
3. The plan's mechanism for the shear response (a gauge defect of P) was wrong, as P6b shows.

## What it means

- **The detector reads rotation; Σ also reads a coordinate convention.** A rest Dirac packet's ⟨Σ⟩ read at equal coordinate time is not invariant under a time relabelling that varies across the packet. A shear in the shift is locally such a relabelling, so Σ drifts by about the packet's orbital share times the shear, fading as packets widen. Total J is nearly free of it, and so is a genuine rotation.
- **The shift sector stands as it was.**
  - The momentum coupling with plaquette spin events reads rotation for spin and orbit alike.
  - It agrees with the exactly gauge resync coupling to within the 7% lattice difference between momentum and energy current.
  - 59's on-axis precession is exact, and the in-plane and gradient results follow from the calibration.
  - No coupling change is indicated: the planned 61 resync coupling would change only that 7%, and would cost 56's drift universality.
- **For later detector tests:** score total J or a frame-transformed spin, or use wide packets. A raw Σ in a sheared region mixes in the coordinate effect.
- **Suggested next (not approved).** The mechanism step: an event rule for how the shift itself propagates, so that offsets travel by coincidences the way N does. Score it with nothing refitted on 59's zero-boundary on-axis precession (now with a fully characterised detector), on 60C's diagonal prediction, and on 55's co-moving pair.

## Addendum 2026-11-28 (after the Codex review; nothing re-scored)

[Review](../../review/EX-GC60_codex_review.md). The verdicts above stand. Three statements in this note went further than the data:

1. **"Fully characterised" and "a readout convention" are too strong.** What was shown:
   - the two-coefficient calibration (rotation plus xy shear) transfers to one new site, the diagonal, at σ = 3;
   - an exactly gauge-equivalent coupling leaves 93% of the raw Σ shear drift, at one width and one orientation.

   That rules out the coupling difference as the main cause of the drift. It does not prove that the whole drift is a coordinate effect. The remaining 7% was not predicted independently, and tilted spins, other gradient components and off-symmetry positions were not tested. Better wording: **the calibrated detector predicts a new site, and most of Σ's shear drift is also present under a pure time relabelling.** EX-GC61D is planned to test this directly.
2. **P5.** The largest orientation difference in the scorecard is 2.45 × 10⁻⁸, not 10⁻¹⁰. It still passes the 1% gate.
3. **P6a.** The plan named the full angle history, but the scoring used only the final angle difference, and the run did not save the gauge-run histories. Record P6a as **"endpoint and slope agree; the full-history gate is unverified"**. The P6a check is also internal, because the same transformation is used in the evolution and at readout. It is not a comparison of a transformed physical detector.

The diagonal result is also narrower than the summary suggests. It discriminates for Σ (rotation alone would miss by 23.5%). For J, rotation alone already passes (−0.05%), so the site gives little leverage on J's small shear coefficient. Both comparisons use the gradients of the solved lattice field, so they test the detector's response to a given field, not the field law itself.
