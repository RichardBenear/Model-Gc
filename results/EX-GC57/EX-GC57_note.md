# EX-GC57 — frame dragging read by light, by orbiting matter and by spin; and 56's light miss closed

2026-11-28 | EXPLORATORY | Claude | PI closure pending | [Plan](../../spec/EX-GC57_plan.md) (pushed before any run) | About 0.014 core-h scored (A 1 s, B 49 s), plus about 0.03 core-h of post-hoc and recompute runs, on the cloud container (Linux, Python 3.13.18, numpy 2.5.3). Scoring: `src/exgc57_summary.py` → `summary.json`; recomputation: `src/exgc57_recompute.py` → `recompute.json`; post-hoc diagnostics: `src/exgc57_posthoc.py` → `posthoc.json`. This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**A rotating source drags light loops and orbiting matter the same way, and the drag is GR's. But an elementary spin feels only about 5–11% of it.** Three of the six registered predictions failed. All three failures were my choice of instrument or comparator, not the physics. The post-hoc checks below show this, and they are not re-scored.

- **56's light miss is closed (P1 PASS).** With an exact one-way lattice mode, light's energy centroid drifts by the mode-energy-weighted −β⟨cos k⟩ to 1.4 × 10⁻⁵, in all six cases (k₀ = 0.25 and 0.5; widths 8, 16 and 32). The envelope centroid matches the |û|² weighting to 3 × 10⁻⁵. 56's launch, propagated exactly, reproduces its 0.9942 and 0.9972. So the 0.6% was the launch (it is not a pure one-way mode) and the readout (an energy centroid weights by ω²|û|²). Light and matter share the drift exactly.
- **Orbiting matter reads the local frame dragging.** A small probe ring (radius 5, L ≈ ±0.97) gives a precession frequency Ω_orb = ΔE/ΔL:
  - Against the lattice curl −½∇×β at its centre: 0.960, 1.070, 0.973 and 1.033 (on axis and in the plane, d = 24 and 32). The differences are the probe's size.
  - Post hoc, against the loop integral of β round a circle of its own radius: 0.997 to 1.003.
  - The torque i⟨[H_β, L]⟩ gives the same Ω to within 0.5% and 3.2% (P4 PASS).
  - It flips exactly with the ring's rotation.
- **Light reads it too.** The synchronised light clock on a 16 × 16 loop reads its own loop integral of β. Post hoc, after dividing by the small clock's instrument factor (1.045, from C3), the odd part agrees with the readout to 0.993–0.998.
- **The shift is GR's frame dragging for the actual source.** Post hoc, the lattice curl at each probe equals the continuum linearised-GR prediction from the ring's own momentum density to 0.04–0.5%. The axis/plane ratio is −1.561 on the lattice against −1.565 in the continuum. My registered −2 was the point-dipole value, which is wrong for a ring of radius 10 seen from 32 (53 had already flagged this).
- **Spin is left out (P5 PASS).** A resting Dirac packet's spin precesses at:
  - 0.109 of the local Ω at σ = 3;
  - 0.048 of it at σ = 5.

  That is 1.12 and 1.06 times 2c_B, where c_B is the small spin part of the lattice momentum density, fitted on the same packet. GR predicts 1. So the momentum coupling is the canonical one. An elementary spin feels only a lattice remnant of frame dragging, and that remnant shrinks as the packet widens. Spin and orbital gyroscopes disagree.
- **Pure-gradient control (P6 PASS).** In β = ∇χ, orbit, spin and light all read below 6 × 10⁻⁴ of their ring values.

| | What was scored | Result | Verdict |
|---|---|---|---|
| P1 | one-way launch: E-centroid / ε-predictor and envelope / \|û\|² within 1e-3; 56's launch reproduces 0.9942, 0.9972 within 5e-4 | max dev 1.4e-5 and 3.0e-5; 0.99424, 0.99721 | PASS |
| P2 | Ω_orb/Ω_curl, Ω_light/Ω_curl, Ω_orb/Ω_light within 10% (4 positions); odd/even in ℓ > 10⁶ | orb 0.960–1.070 ✓; light 0.898–1.253 ✗; orb/light 0.854–1.069 ✗; light odd/even 15–20 in the plane ✗ | FAIL |
| P3 | Ω_orb(axis)/Ω_orb(plane) at d = 32 = −2 within 15% | −1.471 (lattice curl −1.561) | FAIL |
| P4 | Ω_torque/Ω_orb within 10%, on axis | 1.005, 0.968 | PASS |
| P5 | Ω_spin/Ω_curl < 0.3, falls with σ, within 25% of 2c_B | 0.109 (σ = 3), 0.048 (σ = 5); 1.12, 1.06 × 2c_B | PASS |
| P6 | gradient control below 1% of ring values | ≤ 5.9 × 10⁻⁴ | PASS |

## Checks

| | | |
|---|---|---|
| C1 | β = 0 gives zero drift; Fourier propagation of 56's launch vs 56's RK4 | 0; 6.7 × 10⁻⁹ — PASS |
| C2 | rigid rotation: −Σβ·P = −ω·L | 1.3 × 10⁻¹³ — PASS |
| C3 | small wave clock (L = 64, width 4) in uniform β: delay / 2Σβ within 1% | 1.045 — FAIL |
| C4 | ±ℓ give the same q and N; ℓ = 0 shift is zero | 1.4 × 10⁻¹⁵, 8.9 × 10⁻¹⁶, 7.5 × 10⁻¹⁹ — PASS |

**Recompute.**
- The one-way launch, time-stepped by RK4 with 56's equations, against a continuum-Fourier predictor: 0.999991 and 0.999993.
- The probe and spin energies from ⟨ψ|H_βψ⟩ with 53's operator, and the curl from 53's loop code on 1 × 1 plaquettes: every P2 and P5 ratio agrees to 10⁻¹³.
- Shared dependencies: 53's ring sources, static N and shift solver, and the probe builders.

## Why the registered fails are instrument choices (post hoc, not re-scored)

1. **The small light clock is off by 4.5%.** A width-4 pulse on a 64-site loop disperses, and its energy wraps round the loop. The same instrument reads uniform β at 1.045 (C3). That factor, not the physics, is most of the light-versus-curl gap.
2. **The in-plane light loop picks up a part that does not flip with the ring.** Rerun with flat N, that part falls from about 5% to 7 × 10⁻⁵ of the odd part. So it comes from the clock-rate well across the loop combined with the synchronised preparation, and it cancels between ±ℓ.
3. **Light and matter sample different areas.** The 16 × 16 light loop and the radius-5 probe average a curl that varies as 1/r³ over different regions. Each matches the loop integral of its own path to within 0.8%.
4. **P3 used a point dipole.** The ring's own near field gives −1.565, and the lattice gives −1.561.

## Deviations

1. P3's comparator (point dipole −2) ignored the ring's finite size, which 53 had already flagged. The plan should have used the source-integrated continuum prediction.
2. The light-clock size (L = 64, width 4) was chosen without a prior instrument test. C3 shows its 4.5% factor.
3. The plan said the point curl averages 8 cells. The code uses the 4 plaquettes that touch the site (stated in the code header before the run).
4. The plan-stage spin fraction (c_B ≈ 0.10 at σ = 3) was 2c_B in the run's normalisation (c_B = 0.049). The registered prediction used 2c_B fitted in the run, so it is unaffected.

## What it means

- **Rotation is now covered.** The shift sector passes three things:
  - translation: 56's universal drift, with light now exact;
  - curl: one local Ω for light loops and orbiting matter;
  - the field: the lattice shift is GR's linearised g₀ᵢ for the actual source, to 0.5%.

  Frame dragging in this model is the substrate's local rotation, and everything that circulates through it reads that rotation as the loop integral of the drift.
- **The spin gap is a concrete requirement.** GR couples g₀ᵢ to the symmetric (Belinfante) momentum density, which includes ½∇×(spin). Our lattice P is the canonical density, so spin is nearly blind to frame dragging. That is a lattice violation of the equivalence principle between spin and orbit. An event rule for the shift (58) should be checked for which momentum density it counts. If it counts canonical hops, spin stays out. If it counts something like the symmetric density, spin joins.
- **Suggested next.** 58 as approved, the directed-coincidence rule. Score it as planned, and also on whether its momentum density has the spin term (a resting spin packet, as here).
