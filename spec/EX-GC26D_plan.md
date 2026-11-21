# EX-GC26D — moving knot with its own field: is there drag? (EXPLORATORY)

Claude, 2026-11-21. PI go-ahead (Richard, 2026-11-21): "I sign off on EX-GC26C. Go ahead with 26D on my PC." Runs
on the PI's PC (≤ 8 cores; PC allowance about 29–31 of 40 core-hours left). Base: closed 26C code. Massless Skyrme
reference model with a dynamical Nordström scalar, c = 1. Not Model Gc; the construction is ours.

## Question

Model Gc failed because a consumed carrier defines a rest frame, which produces drag (EX-GC17, EX-GC24A). In the
frame-free replacement, a knot moving with its own self-generated field should coast with **no extra drag**. Its field
should be carried along with it, Lorentz-contracted, with no lag behind the knot.

## Initial data

1. **Static coupled state.** Rerun 26C for B = 1 (G_N = 0.01, T = 12) on the PC into `local_raw/EX-GC26D`. This
   regenerates the knot in its own settled well. Its record goes to `results/EX-GC26D/prep_B1.json`; the committed 26C
   record is not overwritten.
2. **Static field.** Take Φ from that state for r ≤ 9.5, blended smoothly over 9.5–10.5 into the fitted −A/r outside.
   This removes the outgoing front. Φ = 0 is kept on the outer wall at 17; the mismatch there emits a wave, but it
   cannot reach the knot before t ≈ 15.
3. **Exact Lorentz transform** at t = 0 of the static knot *and* field. Speed u along +x, starting at x₀ = −uT/2, so the
   knot ends as far right of the centre as it started left. φ and v are transformed as in 26B; for
   the field, Φ(x) = Φ_rest(x′) and Φ̇ = −γu ∂′Φ_rest, using cubic interpolation.

## Runs (h = 0.2, matter box L = 7, field box L_Φ = 17, dt = 0.01, T = 8)

| Run | Description |
|---|---|
| 1 | coupled, u = 0.3 |
| 2 | uncoupled control, u = 0.3: G_N = 0, Φ ≡ 0, boosted 26E1 rest state, same x₀ |
| 3 | coupled, u = 0.1 |
| 4 | uncoupled control, u = 0.1 |

Runs 2 and 4 need no prep, so they run alongside it.

## Readouts and predictions

- **Speed.** From the energy centroid: early (t ∈ [0.5, 2]) and late (t ∈ [6, 8]) fits, giving the relative speed
  change. Also the baryon centroid.
- **Speed ÷ (P/H)**, as in 26B2.
- **Primary: field-induced drag.** The relative speed change of the coupled run minus that of its control. Prediction:
  0 within the run-to-run noise.
  - A drag of the Model Gc type would show as a steady deceleration that grows with G_N and is absent in the control.
  - As a scale: Model Gc particles lost their speed within a coasting time of about 0.04 (EX-GC24A).
- **Field co-motion** at t = 8, on the axis through the knot's centroid X:
  - fore/aft ratio Φ(X + d)/Φ(X − d) for d = 2–5. Prediction 1; a lagging field gives a ratio below 1.
  - axial/transverse ratio Φ(X + d·x̂)/Φ(X + d·ŷ). Prediction 1/γ: 0.954 at u = 0.3, 0.995 at u = 0.1.
- **Conservation.** Total energy (knot + field) to ≲ 1e-9 relative, since the action is time-independent. Total
  momentum (knot + field) drift is reported; it is not exact on the lattice. Geometric degree at start and end.

**Cost:** about 6 runs-worth at about 0.5 core-hours each on the PC: roughly 3 core-hours, about 1.5 hours wall.
