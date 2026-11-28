# EX-GC59 — spin source to spin probe, end to end: a polarised packet's frame dragging makes a second spin precess

2026-11-28 | EXPLORATORY | Claude | PI go-ahead 2026-11-28 04:32 UTC ("you can close 57 and 58. Go ahead with 59."). Cap: 2 core-h (expected about 0.2). This plan is pushed before any run. Our construction; not LaViolette's Model G and not Chain Theory.

## Why

The Codex review of 57/58 found three gaps:
- **No measured precession.** 58's spin results were energy splittings in a ring field, and its spin-source test compared the field with a sum of the same field. The solver was also handed the source's angular momentum at the boundary.
- **No real light–matter comparison.** 57's light versus orbit comparison used different apertures and an uncalibrated small clock.
- **Order independence untested.** It was shown only for one split step, not for asynchronous updates.

59 joins source, field and detector in one experiment, with the prediction fixed before the run:
- a resting, polarised spin packet sources the shift under 58's rule (links plus plaquettes, label PS);
- a separate probe spin, pointing sideways, is evolved in that field;
- its precession angle is read from its spin components over time.

## Setup

**Source.**
- A resting positive-energy packet (m = 0.6, σ = 3, r = 0), spin up or down along z, at the centre of the box.
- Its charge sources N (51's static solver, κ = 4π·12·0.01/Q).
- Its conjugate count sources the shift: j = P + ¼Cᵀs under PS, j = P under links only. The static shift is solved by 53's `poisson` with **zero boundary values**, so the solver is not given J.
- The fields are frozen (the probe is a test particle).

**Probe.**
- The same packet with spin along +x (χ = (1, 1, 0, 0)/√2), built in a 48³ box centred on the probe site. The frozen N and β are cut to that box.
- It is evolved by RK4 under H_m(N) + H_β (+ H_pl for PS), dt = 0.1, T = 10.
- The angle is φ(t) = atan2(⟨Σ_y⟩, ⟨Σ_x⟩). The rotation of the source enters only through the odd part, φ_odd = [φ(+S) − φ(−S)]/2. This removes free spin drift, spreading and the N well, which are all even in the source spin.
- Ω_meas is the least-squares slope of φ_odd(t).

**Positions:** on the source's axis at d = 12 and 16, and in its plane at d = 16.

**Prediction (fixed before the run).**
- Linearised GR from the source's own measured count, in free space:
  - β_GR(x) = −(λκ/4π) Σ j(r')/|x − r'|;
  - Ω_GR(x) = −(λκ/8π) Σ (x − r') × j(r')/|x − r'|³.
- Averaged over the probe's initial density: Ω_pred = Σ_x |ψ_probe(x)|² Ω_GR,z(x).
- Under PS the precession is a rigid rotation by Ω·(L + S), so the probe's spin should turn at Ω_pred. The point-dipole value 2GJ/d³ (axis) and −GJ/d³ (plane), with J = angmom(j), is reported beside it.

**Controls.**
- Source reversed (±S).
- An unpolarised source: the average of the up and down densities, giving a source current of ½(j_up + j_down).
- A β = 0 run.
- Links only (P) for both source and probe.
- A pure gradient β = ∇χ (Gaussian χ with peak |β| at the probe matched to the source's).
- Box and boundary: n = 97 with zero boundaries (main), n = 129 with zero boundaries, n = 97 with 53's dipole boundary from angmom(j).

**Light loop on the same region.**
- A synchronised wave clock (55's `wave`, sg = −1, width 4, k₀ = 0.5) on the square loop of half-side 8 centred at the axis probe site d = 16, normal to z.
- It is calibrated before any source run: the delay in uniform β on the same 64-site loop gives the factor f. 57 measured f = 1.045.
- Ω_light = −(Δt_odd/f)/(4A), with A = 16².
- Predicted value: Ω_pred,loop = −∮β_GR·dl/(2A) on that loop path.

**Asynchronous events.** The PS axis d = 16 probe is also stepped by 58's events, with a new random order of all links each tick.
- Half the events go in one random order, then RK4 for H_m, then the same order reversed.
- Spin events go at random positions in the sequence.
- Two seeds, compared with the Hamiltonian run. Numba is used for the sequential loop.

## Predictions (scored)

- **P1 (spin precesses at GR's rate).** Under PS, with n = 97 and zero boundaries, Ω_meas/Ω_pred = 1 within 3% at axis d = 12, axis d = 16 and plane d = 16. Axis and plane have opposite signs.
- **P2 (reversal and unpolarised source).**
  - φ(+S) − φ₀ = −(φ(−S) − φ₀) within 1% of |φ(+S) − φ₀| at t = T, where φ₀ is the β = 0 run.
  - The unpolarised source moves φ by less than 1% of that at t = T.
- **P3 (links only leave spin out).** |Ω_meas(P)| < 5% of |Ω_meas(PS)| at axis d = 16. Expected about 1%: the source's J and the probe's response are each about 0.1 of PS.
- **P4 (box and boundary).** n = 129 with zero boundaries, and n = 97 with dipole boundaries, are each within 2% of the main Ω_meas at axis d = 16, and within 3% of Ω_pred.
- **P5 (gradient null).** |φ_grad − φ₀| < 1% of |φ(+S) − φ₀| at t = T.
- **P6 (light on the same region).** Ω_light/Ω_pred,loop = 1 within 2%. Its sign matches the spin probe's.
- **P7 (asynchronous order).** For both seeds, Ω_meas(async) agrees with the Hamiltonian Ω_meas within 10⁻³ relative, and with each other within 10⁻³.

## Checks

- **C1.** With β = 0, φ(+S) and φ(−S) are identical (the N field does not depend on the spin direction to 10⁻¹²).
- **C2.** The source's angmom(j_PS) = L + ½⟨Σ_z⟩ (58: 0.4988) to 10⁻¹⁰.
- **C3.** The light calibration factor f is within 1% of 57's 1.045.
- **C4.** The probe norm is conserved to 10⁻⁸ in the RK4 runs, and to 10⁻¹² in the event runs.

## What the outcomes mean

- **All pass.** Spin gravitomagnetism works end to end in this model. A polarised particle's own count sources a shift, which propagates under the shift field law to a second particle and turns that particle's spin at the linearised-GR rate. Light round the same region reads the same rotation. The result does not depend on how the events are ordered.
- **P1 fails with P4 passing.** The coupling and source are consistent with each other, but the precession of a real spin is not the energy-splitting rate (non-rigid response).
- **P4 fails.** 58's 0.1% source agreement depended on handing the solver J at the boundary.
- **P7 fails.** The event rule is not order-free at a finite tick for a spin observable.
- **Not tested:** the shift's own transport law (still 53's field equation), the strength 4 and the ¼. That is 60.

## Outputs

`src/exgc59.py`, `src/exgc59_summary.py`, `src/exgc59_recompute.py`, `results/EX-GC59/`, note, LOG row and HANDOFF section.
