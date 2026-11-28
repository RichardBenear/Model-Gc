# EX-GC60 — what the spin detector responds to: rotation, shear and a gauge test, with the coupling unchanged

2026-11-28 | EXPLORATORY | Claude | PI go-ahead 2026-11-28 13:13 UTC ("yes, start with 60"). Cap: 2 core-h (expected about 0.4). This plan is pushed before any run. Our construction; not LaViolette's Model G and not Chain Theory.

## Why

59's probe spin precessed at the solved field's curl on the source axis. In the plane it ran at 0.775, and a pure-gradient shift moved it by 7.7%. A post-hoc calibration (rigid rotation −0.9976, xy shear +0.093 per unit) predicted both rates, but that calibration was made after the data were seen.

The Codex review asks for a decisive detector test, done before touching the coupling:
- three fields with controlled gradients, at several packet widths and orientations, tracking spin, orbit and their sum;
- a proper gauge (resynchronisation) test with a rotation as positive control;
- a frozen calibration that predicts a new position before it is run;
- the in-plane convergence check.

## The idea (worked before any run)

**1. Time resynchronisation on the lattice.** Relabelling time by t → t + χ(x) changes the state by U = e^{−iG}, with G = ½{χ, H}. To first order it adds the coupling

H_g = −i[G, H] = −(i/2)[χ, H²].

This is the exactly gauge-equivalent coupling for β = −∇χ.
- In the continuum Dirac case H² = p² + m², so H_g = −β·p: the momentum coupling **is** the resync, and a gradient shift is pure gauge.
- On the naive lattice H² = Σ sin²k_a + m², so for uniform β it is H_g = −β sin k cos k = −β·½ sin 2k. That is a two-site hop, which is 55/56's energy current ½{H, v} (−(i/2)[x, H²] = ½{H, i[H, x]}).
- Our coupling −β·P uses sin k. 56 chose it because its drift −β cos k is universal and equals light's. The resync coupling drifts by −β cos 2k, the mass-dependent energy current that 56 excluded.

**2. So the lattice cannot have both.** The coupling with universal drift (P) is not exactly gauge for gradient shifts, and the exactly gauge coupling (the energy current) breaks drift universality. They differ by sin k − ½ sin 2k ≈ k³/2, which matters only at k ~ 1/σ. A gradient's symmetric second derivative is the shear, so under P a probe can respond to shear. 59's 0.093 would then be this gauge defect, and it should fade for wide packets.

**3. Why the spin vector and not total J.** Both couplings are spin-blind, so they act on the packet's orbital part. A rest Dirac packet locks a small orbital share f(σ) = 1 − ⟨Σ⟩ (0.101 at σ = 3, 0.046 at σ = 5; 57, 58) to its spin. The shear defect reshapes the orbital part, and the spin follows through that locking, while J moves little (59: 0.002). The expected scaling is k_sh(σ) ≈ f(σ), about 0.093 → 0.04.

## Parts

**60A. Calibration in uniform-gradient fields** (flat N, 48³ box, probe at the centre, dt = 0.1, T = 10; strength s = 10⁻⁵).
- Probe: rest packet, m = 0.6, r = 0, widths σ = 3, 4 and 5.
- Fields, which have the same gradients everywhere and so are matched over the probe density by construction:
  - rotation β = s(−y, x, 0), with curl 2s and Ω = −s;
  - xy shear β = s(y, x, 0), which is curl-free and equals ∇(s xy);
  - their sum.
- Orientations:
  - xy fields with the spin along x (angle about z);
  - the same fields rotated into the yz plane with the spin along y (angle about x).
- Readouts: the angle of ⟨Σ⟩, of L = angmom(P), and of J = L + ½⟨Σ⟩ in the field's plane. Each is a rate per unit s, relative to a β = 0 run of the same probe.
- PS rule (58) throughout. Links only (P) for shear at all widths and for rotation at σ = 3.

**60B. Gauge test** (σ = 3, xy plane, χ = −s xy so that β_g = −∇χ = s(y, x, 0)).
- (i) Evolve under H_m + H_g, with H_g ψ = −(i/2)(χ H²ψ − H²(χψ)) and H² from applying H_m twice.
- (ii) The gauge prediction: evolve U†ψ₀ freely (β = 0), then apply U = e^{−iG} (Taylor series) at each readout time. The preparation, detector axes and readout times are then all transformed with the field.
- (iii) The P-coupled shear run from 60A.
- Positive control: rigid rotation, which no χ can produce. 60A's rotation runs stand.

**60C. Frozen-calibration prediction at a new position.**
- Setup: 59's source field (n = 97, zero boundaries, PS) at the diagonal site (c + 12, c, c + 12), probe spin along x, ±S odd part.
- The prediction is made before the run from the σ = 3 xy coefficients of 60A, and the probe-density-weighted lattice gradients of β at the site (∂_yβ_x, ∂_xβ_y):
  - Σ rate = k_rot·(∂_xβ_y − ∂_yβ_x)/2 + k_sh·(∂_xβ_y + ∂_yβ_x)/2;
  - the same with the J coefficients.
- The same formula is also applied to 59's axis and plane cases with density-weighted gradients, as a consistency row.

**60D. Convergence.** 59's plane d = 16 case (±S) at dt = 0.05, against dt = 0.1.

## Predictions (scored)

- **P1 (rotation is read exactly).** Under PS, the Σ and J rates in rotation are −1 within 1% for σ = 3, 4, 5 in both orientations. Under links only at σ = 3, the Σ rate is −0.10 ± 0.02.
- **P2 (shear response is the packet's orbital share).**
  - Under PS, k_sh(Σ) falls with σ, and k_sh(Σ)/f(σ) is the same at σ = 3, 4, 5 within 10%, with f = 1 − ⟨Σ_z⟩ of the spin-up rest packet.
  - k_sh under links only equals k_sh under PS within 10⁻⁶ relative.
- **P3 (total J is nearly immune).** |k_sh(J)| < 0.1|k_sh(Σ)| at each width.
- **P4 (linear and additive).** Sum field = rotation + shear within 0.5% for Σ and J.
- **P5 (cubic symmetry).** yz coefficients equal xy within 1%.
- **P6 (the shear response is the gauge defect of P).**
  - (a) The H_g run and the gauge-transformed free run give the same Σ angle history, within 10⁻³ of the P-shear signal at t = T. This checks the implementation of the first-order resync.
  - (b) The Σ rate under H_g is below 10% of the P-shear rate.
- **P7 (prospective).** At the diagonal site, the measured Σ and J rates equal the frozen-calibration predictions within 3%.
- **P8 (convergence).** 59's in-plane Σ rate at dt = 0.05 is within 0.5% of dt = 0.1.

## Checks

- **C1.** Two β = 0 runs of the same probe are identical.
- **C2.** Shear at 2s gives twice the rate within 10⁻³.
- **C3.** 60C's code reproduces 59's axis d = 16 Σ rate to 10⁻¹⁰.

## What the outcomes mean

- **P1–P6 pass.** The probe is a faithful gyroscope for rotation, and its spin vector picks up a finite-packet response to shear. That response is exactly the lattice momentum coupling's departure from a pure resynchronisation, and it fades as packets widen.
  - On this lattice, drift universality (56) and exact gauge invariance of a gradient shift cannot both hold. They coincide in the continuum.
  - 59's plane and gradient failures are then explained, and changing the coupling (the planned 61) would trade them for 56's universality failure.
  - The way forward is either to accept it (the defect is O(1/σ²)), or to find a coupling that is both, for example a P-like hop with a matching spin term.
- **P2 fails (no width dependence).** The shear response is not the orbital-share mechanism, and the diagnosis is open.
- **P6b fails.** The resync coupling also gives shear precession, so the effect is not the gauge defect of P. It may be a readout convention.
- **P7 passes.** The calibrated response predicts precession in a new geometry, so the frozen detector model is complete at first order.

## Outputs

`src/exgc60.py`, `src/exgc60_summary.py`, `src/exgc60_recompute.py`, `results/EX-GC60/`, note, LOG row and HANDOFF section.
