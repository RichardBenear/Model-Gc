# EX-GC59 — spin source to spin probe, end to end

2026-11-28 | EXPLORATORY | Claude | PI CLOSED 2026-11-28 | [Plan](../../spec/EX-GC59_plan.md) (pushed before any run) | About 0.23 core-h scored, plus about 0.2 core-h of post-hoc and recompute runs, on the cloud container (Linux, Python 3.13.18, numpy 2.5.3, numba 0.68.1). Scoring: `src/exgc59_summary.py` → `summary.json`; recomputation: `src/exgc59_recompute.py` → `recompute.json`; post-hoc diagnostics: `src/exgc59_posthoc.py` → `posthoc.json`. This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**A polarised packet's own count sources a shift that makes a second spin precess. On the source's axis the rate is linearised GR's to 1.3–2.3%, and it equals the curl of the solved field to 0.03–0.1%.** Four of the seven registered predictions passed and three failed, and one check failed. All three failures trace to one effect the plan missed, or to the box: **the probe's spin also responds to the shear of the shift, not only its curl.** That response is about 0.09 per unit shear, the packet's orbital share of its angular momentum, carried by the canonical link coupling. GR's point gyroscope has no such term.

- **On axis (P1, partly).** Measured rate over predicted rate is 0.987 at d = 12 and 0.977 at d = 16, inside the 3% band. Recomputed against the lattice field the probe actually sits in, the ratios are 0.9990 and 0.9998. The shortfall from free-space GR is the field in a 97³ box with zero boundaries (lattice curl / free GR = 0.988 and 0.977).
- **In the plane (P1 FAIL).** The ratio is 0.775. Post hoc, two flat-field calibrations of the probe were made:
  - in a rigid-rotation shift its spin turns at −0.9976 of the rate;
  - in a curl-free xy shear it turns at +0.093 per unit shear.

  Combined with the lattice gradients of β at the probe, these predict the measured in-plane rate to 0.04% (axis: 0.03%). On the axis the shear vanishes by symmetry; in the plane it is 1.5 times the curl and opposes it. The probe's total angular momentum J = L + ½Σ responds to shear only at 0.002 per unit, and its in-plane precession is 1.055 of GR.
- **Reversal and unpolarised source (P2 PASS).** The odd part flips to 3 × 10⁻¹¹, and an unpolarised source leaves the probe alone to 2 × 10⁻¹³.
- **Links only (P3 PASS).** Without plaquettes the probe turns at 0.0099 of the PS rate: about 0.1 of the source's J times about 0.1 of the probe's response.
- **Box and boundary (P4 FAIL).** On the axis at d = 16:

  | Box and boundary | Ω / main | Ω / GR prediction |
  |---|---|---|
  | n = 97, zero boundaries (main) | 1 | 0.977 |
  | n = 129, zero boundaries | 1.014 | 0.991 |
  | n = 97, dipole boundary | 1.025 | 1.001 |

  The zero-boundary deficit shrinks by 2.5× from n = 97 to 129, close to the (129/97)³ = 2.35 expected of boundary images. So it is a finite-box effect that vanishes for large boxes. The 2% box gate was too tight for n = 97.
- **Gradient null (P5 FAIL).** β = ∇χ moved the spin angle by 7.7% of the source's effect. A pure gradient has no curl but does have shear. Here the gradient field's amplitude was matched to the largest β in the probe's 48³ box, which includes the region near the source, so at the probe it is 1.2 × 10⁻³, against zero for the source on its axis. The shear coefficient above predicts 13% from the point values, against the measured 8%. Gradient shifts are therefore not pure gauge for this spin readout. That is the same lattice fact 56 noted (P ≠ J_E), now seen in spin.
- **Light on the same region (P6 PASS).** The calibrated light clock on the 16 × 16 loop around the probe gives 0.988 of the free-space GR loop prediction, with the same sign as the spin. Against the solved lattice loop integral it gives 1.023.
- **Asynchronous events (P7 PASS).** With every link and site event in a new random order each tick, the precession rate agrees with the Hamiltonian run to 1.2 × 10⁻⁶, and the two seeds agree to 3 × 10⁻¹¹.

| | What was scored | Result | Verdict |
|---|---|---|---|
| P1 | Ω_meas/Ω_pred within 3% (axis 12, axis 16, plane 16); opposite signs | 0.987, 0.977, 0.775; signs opposite | FAIL (plane) |
| P2 | reversal and unpolarised within 1% | 3e-11, 2e-13 | PASS |
| P3 | links-only / PS < 5% | 0.0099 | PASS |
| P4 | n129 and dipole within 2% of main and 3% of prediction | 1.014 ✓, 1.025 ✗; 0.991 ✓, 1.001 ✓ | FAIL |
| P5 | gradient null < 1% | 7.7% | FAIL |
| P6 | light / loop prediction within 2%, same sign | 0.988, yes | PASS |
| P7 | async vs Hamiltonian, seeds within 1e-3 | 1.2e-6, 3e-11 | PASS |

## Checks

| | | |
|---|---|---|
| C1 | β = 0: up and down source N give the same probe | 2.5e-16 — PASS |
| C2 | source angmom(j_PS) = L + ½⟨Σ_z⟩ | 1.6e-15 relative — PASS |
| C3 | light calibration f within 1% of 1.045 | 1.0451 — PASS |
| C4 | probe norm < 1e-8 (RK4), < 1e-12 (event runs) | 1.07e-7 for both — FAIL |

C4: both kinds of run use RK4 for H_m at dt = 0.1, and RK4 is not exactly unitary. The bound in the plan was wrong for the event runs, which also contain an RK4 step. The norm loss is the same 1.07 × 10⁻⁷ with and without the shift, so it does not enter the odd part.

**Recompute.**
- The axis runs, repeated at dt = 0.05 with the spin read by a different routine and slopes by polyfit, agree to 2 × 10⁻⁵ and 3 × 10⁻⁶.
- The probe-density average of the lattice curl of the solved shift is matched by the measured rate to 0.9990 and 0.9998.
- The light loop's lattice integral is 0.966 of the free-space loop prediction.
- Shared dependencies: the source builder, the static solvers and the probe builder.

## Deviations

1. **Crash and rerun.** The first `main` run crashed on its β = 0 control: the curl was taken of a missing field. I fixed the code and reran `main` in full. The axis d = 12 value before the crash was identical (0.98674).
2. **Light predictor bug.** The stored light predictor was −∞, because a loop link coincided with a source link carrying the source's exponential tail. I fixed it by dropping that singular self-term (`exgc59_fixpred.py`) before scoring. The wave delays were not rerun.
3. **Gradient amplitude.** It was matched over the probe's whole 48³ box, as in 57, not at the probe as the plan worded it. The null field was therefore far stronger at the probe than the source's field.
4. **C4's bound** did not allow for the RK4 step inside the event runs.
5. **The plan's assumption** that the precession is a rigid rotation by Ω·(L + S) ignored the shear of β. That is the main lesson of this study.

## What it means

- **The chain works for what GR calls a gyroscope.** Polarised source, the shift field law, then probe spin: on the axis, where the shift is a pure local rotation, the probe precesses exactly at the curl of the field the source made. That field is linearised GR's up to the box. Light round the same region reads the same sense and size, and event order does not matter.
- **A new lattice-specific effect.** Off the axis, the probe's spin vector also turns under the shift's shear, at about 9% per unit shear. This comes through the canonical link coupling acting on the packet's orbital part, so it would be the same under links alone. The probe's total angular momentum is nearly immune to it (0.2%). GR has no shear coupling for a point gyroscope, because the symmetric part of ∂β is a coordinate effect. So a lattice Dirac spin's Σ is not a clean gyroscope off the axis: it mixes in a gauge-dependent part, as does the gradient null.
- **Suggested next (not approved).**
  1. Remove the shear response at the source. If GR's symmetric momentum density is used throughout (canonical P plus the ¼∇×s term, already in), the leftover shear coupling must come from P ≠ J_E. A shift coupling built from the time-resynchronisation generator would make β = ∇χ exactly pure gauge on the lattice. That would fix P5 and the plane together. Test: the same three uniform calibrations, then 59's plane and gradient cases.
  2. Then 60 as discussed: an event rule for the shift's own transport, scored on 59's on-axis observable with nothing refitted.
