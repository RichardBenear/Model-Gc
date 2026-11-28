# EX-GC64 — the per-local-tick count transfers; reciprocity costs about 7% of strength; the hybrid spin source holds on axis

2026-11-28 | EXPLORATORY | Claude | PI closure pending | [Plan](../../spec/EX-GC64_plan.md) (pushed before any run; PI go-ahead 18:28 UTC, cap 2.5 core-h) | 1.81 core-h on the cloud container (Linux, Python 3.13.18, numpy 2.5.3; 4 processes). Scoring: `src/exgc64_summary.py` → `summary.json`; recomputation: `src/exgc64_recompute.py` → `recompute.json`. Responds to `review/EX-GC63_codex_review.md`. This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**The frozen 63 rule transfers.**
- Against the momentum rule P at the same setting, its offset is ΔR = +0.0054 to +0.0064 at three well depths, two widths and both boost directions.
- The offset is depth-independent to 0.0004.
- It falls by the expected factor of about 4 at half the speed, and does not move at half the time step (10⁻⁶).
- The counted strength is 3.852 at every depth. 54's N-weighted count falls from 3.60 to 3.37 to 2.95 as the well deepens.

That settles the cause of 54's deficit: it was the clock-well weighting.

**Making it one interaction costs strength.** When the same flat current both sources the shift and is felt (RJE1), R = 1.0093, which is +0.0106 above P (predicted +0.0086 ± 0.004). A reciprocal counted interaction is consistent, but carries λ_eff ≈ s²/4 ≈ 3.71 instead of 4. To be Lorentz-invariant it would need its multiplier at about 4 × (4/s)² ≈ 4.3. Alternatively, the coupling stays canonical momentum P (55's choice), and then the counted source does not conserve energy (55: 4.2 × 10⁻³).

**Hybrid spin (counted orbital current + 58's ¼Cᵀs).** It precesses 59's probe at 0.9918 (axis, d = 12) and 0.9923 (axis, d = 16) of the PS rate, against the frozen 0.993 ± 0.003. In the plane at d = 16 it gives 0.9857, which is a FAIL.

| | What was scored | Result | Verdict |
|---|---|---|---|
| P1 | counted strength 3.85 ± 0.02 at κ₀/2, κ₀, 2κ₀ | 3.852 at all three (σ = 4: 3.788, reported) | PASS |
| P2 | ΔR = 0.02943 (4 − s) ± 0.003: base, shallow, deep, narrow, reversed | +0.0055, +0.0054, +0.0057, +0.0064, +0.0054 against +0.0044 (narrow +0.0063) | PASS |
| P3 | V = 0.1: ΔR = ΔR(0.2)/4 ± 0.0015 | +0.0012 against +0.0014 | PASS |
| P4 | dt/2 changes R by < 0.001 | 1.4 × 10⁻⁶ | PASS |
| P5 | RJE1 ΔR = +0.0086 ± 0.004, \|R − 1\| < 2% | +0.0106, 0.93% | PASS |
| P6 | hybrid / PS precession 0.993 ± 0.003 (axis 12, axis 16, plane 16) | 0.9918, 0.9923, 0.9857 | FAIL (plane) |

## Checks

| | | | |
|---|---|---|---|
| C1 | 63's C2 revisited (reported) | see below | — |
| C2 | F₁(0): D2-E1 = P to 1e-3 | 5.1 × 10⁻⁵ | PASS |
| C3 | N = 1: RJE1 = 55's RJE, D2-E1 = 55's MIX source | 0, 0 | PASS |
| C4 | hybrid moment / PS = 0.9930 to 1e-4 | 0.99304 | PASS |

**C1. 63's failed velocity-flip check was the comparison, not rounding.**
- A mirror-reflected packet (x → −x, spinor transformed by α_yα_zβ) gives the mirrored source link by link to 1.8 × 10⁻¹⁶. So the counted source is exactly mirror-symmetric.
- The mirror of a spin-up packet is spin-down (a reflection is parity times a π rotation). 63 compared spin up at −k with spin up at +k, which is not a symmetry pair.
- With compensated summation the totals still differ: by 7 × 10⁻⁷ in the n = 49 box, against 1.5 × 10⁻¹¹ at n = 65. So it is a small spin-dependent lattice part of the total current that shrinks with the box, not summation error.

**Recompute.**
- An independent least-squares fit reproduces every R to 10⁻¹⁵.
- Fitting only t ≥ 2 shifts every R up by about 0.0015, with every ΔR within 0.0003 of the full-window value. The offsets are robust to the window, as the Codex review asked.
- Precession ratios from a separate fit of the saved angles agree to 10⁻¹².
- Shared dependencies: records, packet builders, operators, evolve routines, probe Hamiltonian.

## Deviations and caveats

1. **ΔR runs about 0.001 above the line prediction everywhere.** It is inside the band, but it is systematic. A likely cause (not checked) is the shape difference between the counted current and P, which the line ignores.
2. **The momentum rule's own R depends on well depth**: 1.0031 (κ₀/2), 0.9986 (κ₀), 0.9897 (2κ₀). This was not predicted and is not scored. At fixed λ = 4, the P rule's Lorentz invariance in this box shifts by about 1.3% per doubling of depth. Candidates are post-Newtonian terms of order κ² (the field law is linear in the source but the dynamics are not) and the n = 49, T = 8 box. It does not affect ΔR, which is depth-independent. But it means "λ* = 3.91" belongs to one depth, and the next pair test should separate κ² terms from the box.
3. **P6's plane failure is not diagnosed.** On axis the shear of β vanishes by symmetry; in the plane, Σ responds to shear at 0.093 per unit (59/60). The hybrid changes the orbital, near-field part of the source, and so the shear at the probe, by more than its dipole. Reading J (total angular momentum), which is nearly shear-blind (60), would test that. 59's pipeline records only Σ.

## What it means

- **The orbital frame-dragging source is a count of directed signals in local ticks.** That is now established as a transferable rule in this lattice, not a one-setting coincidence.
- **The cost of one interaction is now quantified.** Reciprocal counting gives about 3.71, not 4. Either the counting multiplier is not exactly 4, the coupling is not the same count, or "count × velocity" should be a different operator (p rather than E·v). This sharpens the λ question: λ = 4 fits the momentum operator; the energy-current count falls about 4% short at each end.
- **Spin coexists with the counted orbital source** at the measured normalisation on axis, with an undiagnosed 0.7% in the plane.
- **Suggested next (not approved):**
  1. A directed-coincidence rule that tallies events, so that 4 is measured, not multiplied in. Test whether its natural operator is p (λ = 4 reciprocal) or E·v (3.71).
  2. Separate the P rule's depth dependence into κ² terms and box (a larger box at 2κ₀).
  3. Re-read P6's plane with J.
