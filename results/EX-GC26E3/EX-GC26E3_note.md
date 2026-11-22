# EX-GC26E3 — does a spinning knot's internal energy fall and resist force like rest energy? (EXPLORATORY)

Claude, 2026-11-22. PI go-ahead (Richard, 2026-11-22): "Go ahead with the GC26E3 test." **PI CLOSED 2026-11-22** ("You can close GC26E3.").

- **Plan:** `spec/EX-GC26E3_plan.md`. Plan and LOG row 4c286fe, pushed about a minute after stage 1 started (see
  Protocol). Stage 2 predictions and criterion c7e5563, pushed before the fall runs.
- **Code:** `src/exgc26e3.py`, which reuses `src/exgc26e1.py` unchanged.
- **Separate recomputation and figure:** `src/exgc26e3_figure.py` → `recompute.json` and
  `paper/figures/exgc26e3_spin.png`.

Massless Skyrme reference model with the 26E1 Nordström coupling, not Model Gc. The construction is ours.

## Conclusion

**The registered readout fails. The momentum ledger and the energy centroid say the spinning knot's extra energy has
both weight and inertia.**

- **Registered readout (baryon centroid): NOT SUPPORTED.** The spinning knot's baryon centroid falls 8% (ω = 0.5) and
  13% (ω = 0.7) faster than the static knot's. That lies closer to the alternative "internal motion adds weight but no
  inertia" than to the energy-coupling prediction, at both spin rates.
- **The momentum readouts say the opposite:**
  - **Weight.** The force on the spinning knot equals the energy-coupling prediction to 0.1%. Its gravitational charge
    rises with its spin energy, as von Laue's theorem implies.
  - **Inertia.** Momentum divided by the velocity of the energy centroid equals the knot's total energy, spin
    included, to 0.5%.
  - The energy centroid falls at 1.045 and 1.060 of the static rate, closer to the energy-coupling prediction (1.016,
    1.014) than to either alternative.
- **Why the two disagree: the baryon centroid separates from the energy centroid.** Once the knot spins in the field,
  the baryon (topological-charge) centroid moves away from the energy centroid, increasingly with time: by 0.017
  (ω = 0.5) and 0.027 (ω = 0.7) at t = 5. Without spin the separation settles at 0.007.
  - The baryon centroid is not the knot's centre of inertia once it spins. That is a post-hoc reading, not tested.
  - Under the registered rule the result stays "not supported".

## Stage 1: the spinning knot without a field (predictions fixed before the fall runs)

All three spin rates are stable to T = 5. The degree stays −1, energy is conserved to ≤ 3e-11, and the knot drifts by
≤ 2e-7. The spin-up is sudden, so the knot also breathes; K₂ varies by 6–17%.

| ω | kinetic fraction | R_lat (energy coupling) | R_nw (no weight) | R_ni (no extra inertia) | R_ctl (control) |
|---:|---:|---:|---:|---:|---:|
| 0.35 | 3.2% | 1.011 | 0.968 | 1.045 | 1.201 |
| 0.5 | 6.4% | 1.016 | 0.936 | 1.085 | 1.385 |
| 0.7 | 11.8% | 1.014 | 0.882 | 1.149 | 1.678 |

The time-averaged lattice source 2⟨E₂ − K₂⟩ tracks the total energy. The predicted energy-coupling ratio stays within
1.6% of 1, while the kinetic energy rises to 12%. This is the continuum virial (von Laue) result, realised dynamically on
the lattice.

## Stage 2: fall (ratio R to the static 26E1 knot; static a/g = −0.920 baryon, −0.887 energy centroid, t ≥ 1)

| ω | readout | R, fit t ≥ 1 | R, fit t ≥ 2 | R_lat | R_ni | nearer to |
|---:|---|---:|---:|---:|---:|---|
| 0.5 | **baryon centroid (registered)** | **1.080** | 1.091 | 1.016 | 1.085 | R_ni |
| 0.5 | energy centroid | 1.045 | 1.071 | 1.016 | 1.085 | R_lat |
| 0.7 | **baryon centroid (registered)** | **1.134** | 1.145 | 1.014 | 1.149 | R_ni |
| 0.7 | energy centroid | 1.060 | 1.099 | 1.014 | 1.149 | R_lat |
| 0.7 control | baryon / energy centroid | 1.879 / 1.698 | 2.040 / 1.725 | R_ctl = 1.678 | | R_ctl (check passes) |

The control coupling gives the large, predicted change, so the method detects a non-energy coupling.

## Momentum ledger (separate recomputation, `recompute.json`)

| run | force ÷ (g E_total) | predicted | dP/dt ÷ force | momentum ÷ (energy-centroid velocity × E_total) | same, with baryon centroid |
|---|---:|---:|---:|---:|---:|
| static (26E1) | 0.935 | 0.936 | 0.963 | 1.015 | 0.978 |
| ω = 0.5 | 0.950 | 0.951 | 0.980 | 1.005 | 0.936 |
| ω = 0.7 | 0.948 | 0.949 | 0.990 | 0.998 | 0.899 |

- **The force is exactly what the energy coupling gives.** The spin energy counts in the weight.
- **The momentum gained per unit impulse is 0.96–0.99.** This is the known lattice transfer factor of 26B3/26E, a little
  better with spin.
- **The inertia read with the energy centroid equals the total energy, spin included.** Read with the baryon centroid it
  comes out 6–10% low.

## What this does and does not show

- **Shown.** In this energy-coupled scalar theory, a spinning knot's gravitational charge rises with its internal
  kinetic energy, time-averaged, as von Laue's theorem implies. This holds to 0.1% on the force. Its momentum per
  energy-centroid velocity is its total energy.
- **Not shown under the registered readout.** The registered readout was the baryon centroid, as in 26E1/E2, and it
  disagrees. The cause is untested. Two candidates:
  1. The rotating knot's charge distribution is redistributed in the field gradient.
  2. Sudden spin-up radiates massless waves, which carry energy, and so shift the energy centroid, but carry no
     baryon number.
  A slow spin-up, or a relaxed (stationary) isospinning state, would separate these.
- **One resolution** (h = 0.2), with a sudden spin-up (breathing included).

## Protocol

- **Late push.** The plan and LOG row were committed (4c286fe) before stage 1 started. Because of a merge with
  ChatGPT's concurrent 27A/B commits, they reached GitHub about a minute after the stage 1 runs had begun. Stage 1 makes
  no fall measurement; the stage 2 criterion was pushed before the fall runs.
- **Cost.** Each run cost about twice the plan's estimate, because of the extra measurements.

## Compute

Six runs of 697–746 cpu-s, 4,336 cpu-s in total = **1.20 core-hours** (cloud; Python 3.13.18, NumPy 2.5.3, SciPy
1.18.1, Linux). End states are kept in `local_raw/EX-GC26E3`, with hashes in the records.
