# EX-GC26D — moving knot with its own field: is there drag? (EXPLORATORY)

Claude, 2026-11-21. PI go-ahead (Richard, 2026-11-21): "I sign off on EX-GC26C. Go ahead with 26D on my PC." **PI CLOSED 2026-11-21.** Plan
`spec/EX-GC26D_plan.md` (df01310, amended before any PC run in ef866e8 and 820ac4f). Code `src/exgc26d.py`, driver
`src/exgc26d_run.py`; records from the PI's PC in 6dec330. Post-run field check `review/EX-GC26D_check/field_only.py`;
figure `src/exgc26d_figure.py` → `paper/figures/exgc26d_drag.png`.

Massless Skyrme reference model with a dynamical Nordström scalar, not Model Gc. The construction is ours.

## Conclusion

**A knot moving with its own gravitational field coasts with no field-induced drag, and the field moves with it.**

- **No drag.** The coupled knot's speed changes between early and late times exactly as the no-field control does.
  The difference is −7e-5 at u = 0.1 and +5e-6 at u = 0.3, over 8 time units. A Model Gc-type drag would have removed
  the speed within about 0.04 time units (EX-GC24A). The rest-frame problem that sank Model Gc is absent here.
- **The field is carried along, Lorentz-contracted** (from a clean field-only check, see below). Fore/aft ratio
  1.000–1.013 (prediction 1). Axial/side ratio 0.949–0.957 at u = 0.3 (prediction 1/γ = 0.954) and 0.992–0.997 at
  u = 0.1 (prediction 0.995).
- **A flaw in my plan spoiled the pre-registered co-motion readout at d ≥ 3.** The PC runs show the field weaker ahead
  of the knot than behind (fore/aft down to 0.76 at d = 5, u = 0.3). That is the far-field taper I added, not lag.
  Details below. The drag result is not affected.
- **Bookkeeping is clean.** Total energy (knot + field) conserved to ≤ 5e-10. The degree stays −1 in all runs.

## Speed (primary readout; energy centroid, fits t ∈ [0.5, 2] and [6, 8])

| u | run | early → late relative change | field drag (coupled − control) | mean speed | speed ÷ (P/H) |
|---|---|---:|---:|---:|---:|
| 0.1 | coupled G_N = 0.01 | −2.99e-3 | **−6.9e-5** | 0.100113 | 0.9865 |
| 0.1 | control (no field) | −2.92e-3 | | 0.100085 | 0.9857 |
| 0.3 | coupled G_N = 0.01 | −2.36e-3 | **+5.0e-6** | 0.296903 | 0.9816 |
| 0.3 | control (no field) | −2.36e-3 | | 0.296761 | 0.9808 |

- **The baryon centroid agrees.** Its drag is +4.4e-5 at u = 0.1 and +1.5e-5 at u = 0.3.
- **The −0.2% to −0.3% speed change is the lattice, not the field.** It is identical in the control and comes from
  grid pinning and the stencil (26A/26B).
- **The coupled knot is 0.03% (u = 0.1) and 0.05% (u = 0.3) faster than its control.** Its speed ÷ (P/H) is
  slightly higher as well. The field carries 0.34% of the total momentum. This is a small offset in the start state,
  not a trend.
- **A separate recomputation** (endpoint finite differences instead of fits) gives drags of −8e-5 and 0.

## Field co-motion at t = 8 (Φ on lines through the knot's centroid)

| u | d | PC run fore/aft | field-only, PC setup | field-only, exact far field | axial/side (exact) | 1/γ |
|---|---:|---:|---:|---:|---:|---:|
| 0.1 | 2 | 1.007 | 1.006 | 1.006 | 0.997 | 0.995 |
| 0.1 | 3 | 1.005 | 1.004 | 1.005 | 0.992 | |
| 0.1 | 4 | 0.959 | 0.959 | 1.001 | 0.995 | |
| 0.1 | 5 | 0.836 | 0.836 | 1.000 | 0.995 | |
| 0.3 | 2 | 0.986 | 0.984 | 1.006 | 0.957 | 0.954 |
| 0.3 | 3 | 0.916 | 0.914 | 1.009 | 0.952 | |
| 0.3 | 4 | 0.826 | 0.824 | 1.013 | 0.949 | |
| 0.3 | 5 | 0.759 | 0.756 | 1.001 | 0.955 | |

**What went wrong.**

- The plan tapered the −A/r tail to 0 between r = 12 and 16 (rest frame) so the field met the Φ = 0 wall smoothly.
- That taper is not a static solution, so it sends a disturbance inward at the speed of light.
- I checked only that the disturbance could not reach the knot itself before t ≈ 11. By t = 8 it has reached the
  readout points at d = 3–5.
- Boosted, the disturbance is not symmetric front to back, so it shows up as a fake fore/aft asymmetry.

**The check (post-run, not registered).**

- **Setup.** The field equation alone, with the knot's source moving rigidly at u and Lorentz-contracted. There is no
  back-reaction, which is about 1e-4 here.
- **"PC setup" case.** It uses the PC run's initial data and reproduces the PC readouts to within 0.003. That shows the
  knot's own dynamics play no role in the asymmetry.
- **"Exact far field" case.** It drops the taper and holds the wall at the exactly boosted −A/r. The asymmetry then
  disappears, and the axial/side ratio matches 1/γ at every d.

## Conservation

| run | energy drift | total momentum drift ÷ E | degree start → end |
|---|---:|---:|---|
| u = 0.1 coupled / control | 9.4e-11 / 8.6e-11 | 2.2e-4 / 2.1e-4 | −1 → −1 |
| u = 0.3 coupled / control | 5.3e-10 / 5.1e-10 | 8.5e-4 / 7.8e-4 | −1 → −1 |

- **Momentum is not exactly conserved on the lattice**, as in 26A/B. The drift is almost the same with and without the
  field.
- **The field's own energy falls by 0.6% (u = 0.1) and 1.2% (u = 0.3) over the run** (from 6.2e-3 and 6.3e-3, with
  E ≈ 1.26). This is most likely the taper's disturbance. It is included in the conserved total.

## What this does and does not show

- **It shows the property Model Gc lacked.** A particle together with its self-generated gravitational field has no
  preferred rest frame. It coasts, and its field is simply the boosted static field. This is what the action's Lorentz
  invariance implies; the run shows the lattice dynamics realise it to 1e-4.
- **Only weak coupling** (G_N = 0.01, field energy 0.5% of E) and one resolution (h = 0.2).
- **No acceleration.** A uniformly moving knot does not radiate. Two bodies (26F) remain untested.
- **The co-motion result at d ≥ 3 comes from the field-only check.** It is not from the registered PC readout, which the
  taper spoiled.
- **Scalar gravity limits remain**: no light bending.

## Compute

- **PC:** 1.82 core-hours. Four runs of 1109–2171 cpu-s in parallel, 37 minutes of wall time. Python 3.13.2, NumPy
  2.2.4, SciPy 1.16.1, Windows 11. MSVC-built kernel.
- **Cloud:**
  - smoke tests about 0.05 core-hours;
  - field-only check 0.83 core-hours (four runs of about 745 s). Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1, Linux.
- **Inputs:** the closed 26C end state and the 26E1 rest state, hash-checked (`inputs.json`). End states stay on the PC
  in `local_raw/EX-GC26D`, with hashes in the records.

## Next (PI decision)

- **26F:** two knots attracting through their self-generated field. This is the first test with acceleration.
- **26E3:** isospinning or breathing B = 1, as an equivalence check with internal motion.
