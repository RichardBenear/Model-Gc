# EX-GC26F — two knots attracting through their self-generated field (EXPLORATORY)

Claude, 2026-11-21. PI go-ahead (Richard, 2026-11-21): "go ahead with 26F". Runs on the PI's PC (≤ 8 cores; 4 used).
Base: the closed 26C coupled solver and the 26D input arrays (closed 26C B = 1 end state at G_N = 0.01, hash-checked).
Massless Skyrme reference model with a dynamical Nordström scalar, c = 1. Not Model Gc; the construction is ours.

## Question

Do two knots, each sourcing its own field, attract each other with the Newtonian strength that 26C (active mass) and
26E (passive mass) imply, and does the attraction fall off as 1/d²? This is the first test in which the field
accelerates a body that it did not start with.

## The non-gravitational force, and how it is removed

The knots also push on each other directly through their Skyrme tails.

- **Choice of orientation.** A static pre-check (`review/EX-GC26F_prep/pair_static.py`, product ansatz of the 26C
  knot, L = 14) found the following interaction energies.

  | d | same orientation | π about an axis ⊥ the separation | π about the separation |
  |---|---:|---:|---:|
  | 6 | +2.8e-3 | −6.8e-3 | +1.23e-2 |
  | 8 | +1.4e-3 | −3.0e-3 | +4.6e-3 |
  | 10 | +8.0e-4 | −1.3e-3 | +1.9e-3 |

  Same orientation cancels the leading dipole term, so it is used.
- **The force that remains is not small.** It is repulsive, about 2 times the predicted gravity at d = 6 and about
  1.5 times at d ≈ 9.
- **Comparative design.** Each coupled run has a control with the identical matter start, G_N = 0 and Φ ≡ 0. The
  gravitational signal is the coupled separation minus the control separation.
- **What subtraction leaves.** Gravity also changes the Skyrme force slightly, through W = e^{2Φ} ≈ 1 − 0.028 at each
  knot. That leaves a systematic error of order |Φ| × (Skyrme/gravity), about 3–6% of the signal.

## Initial data

1. **Matter.** The 26C knot placed at x = −d/2 and +d/2 (same orientation, quaternion product ansatz, each copy taken
   from its L = 7 box and vacuum outside it), at rest. The product ansatz is not exactly static. Its relaxation is
   shared with the control.
2. **Field (coupled runs only).** The sum of the two static rest fields: the 26C Φ for r ≤ 9.5, blended into −A/r over
   9.5–10.5 (A fitted as in 26D). The field-box wall is held at the exact superposed static value, so there is no
   far-field taper. That taper is the 26D lesson: it spoiled the co-motion readout there. Both knots stay ≥ 11 from the
   wall.
3. **Boxes.** Matter half-width d/2 + 7; field half-width matter + 4. d/2 lies on the grid.

## Runs (h = 0.2, dt = 0.01, T = 10, G_N = 0.01)

| Run | d | G_N | matter box | field box |
|---|---:|---:|---:|---:|
| 1 | 9.2 | 0.01 | ±11.6 | ±15.6 |
| 2 | 9.2 | 0 (control) | ±11.6 | — |
| 3 | 6.0 | 0.01 | ±10 | ±14 |
| 4 | 6.0 | 0 (control) | ±10 | — |

All four run in parallel (driver `src/exgc26f_run.py`, code `src/exgc26f.py`).

## Readouts and predictions

- **Separation.** D(t) = the difference of the x-centroids of the energy density on the two halves x > 0 and x < 0.
  Also the same for the baryon density.
- **Primary: gravitational relative acceleration.** Fit a quadratic to D_coupled − D_control over t ∈ [0, 10], and
  again over [1, 10].
  - **Prediction:** a_rel = −2 (2E₂/E) A / d², using A = 0.011098 (26D fit) and 2E₂/E = 0.936 (26C record). The
    1/d² is averaged over the control's D(t).
  - That gives −5.8e-4 at d = 6 and −2.5e-4 at d = 9.2. The relative displacement by t = 10 is about −0.029 and
    −0.012.
  - **Supported** if measured ÷ predicted lies within 0.90–1.10 at both d for the baryon centroid. Allowed for: the
    Skyrme-force cross term, product-ansatz relaxation, lattice.
- **Inverse square.** a(6) ÷ a(9.2): predicted 2.35 (the control's D(t) average is used). A 1/d force would give
  1.53, so the two are well separated.
- **Ledgers.** Total energy (knot + field) to ≲ 1e-9 relative. Total momentum should stay near 0 by symmetry, and its
  maximum is reported. Geometric degree −2 at start and end. Φ at each knot.
- **Smoke test (cloud, T = 0.5 at d = 6; T = 0.25 at d = 9.2; outputs discarded).**
  - Energy conserved to 1e-11.
  - Degree −2.
  - By t = 0.5, D_coupled − D_control has fallen by 7e-5 (baryon 6e-5). The prediction at that time is 7.2e-5.
  - This is a smoke test, not a result.

**Cost:** cloud timings scaled to T = 10.

| run | cloud core-hours |
|---|---:|
| d = 9.2 coupled | ~1.8 |
| d = 9.2 control | ~1.3 |
| d = 6 coupled | ~1.0 |
| d = 6 control | ~0.75 |
| **total** | **~5** |

This comes from the PC allowance (about 27–29 of 40 core-hours left), with about 2 hours of wall time on 4 cores.
Cloud use so far: about 0.25 core-hours (pre-check and smoke tests).
