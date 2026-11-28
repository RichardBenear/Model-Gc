# EX-GC57 — frame dragging read by light, by orbiting matter and by spin, and the light-drift readout closed

2026-11-28 | EXPLORATORY | Claude | PI go-ahead 2026-11-28 ("Your 57 plan and 58 plan are very good. go ahead with them."). Cap: 2 core-h (expected under 0.3). This plan is pushed before any run. Our construction; not LaViolette's Model G and not Chain Theory.

## Why

- 56 showed that the link shift is a drift of the substrate: −β·P is a real directed transfer, and in a uniform shift it carries every kind of matter and light by the same −β⟨cos k⟩. One registered miss stayed open: light's drift was 0.9942 of the prediction at k₀ = 0.25. Codex reproduced it in Fourier space and traced it to the launch (not an exact one-way mode) and to the weighting of the predictor.
- A uniform shift is only the translation part. The physical content of g₀ᵢ near a rotating body is its curl, which drags the local frame: a gyroscope precesses at Ω_GR = −½∇×β (Lense–Thirring; for our dipole β = −2G(J×r)/r³ this is G[3(J·r̂)r̂ − J]/r³).
- If the shift is one drift of the substrate, then light going round a small loop, an orbiting matter rotor, and a spinning particle should all read the same local Ω. That is the rotational form of 56's equivalence test, and it can fail.

## The idea (worked before any run)

1. **A rotor reads the loop integral of β.** For matter, the coupling is E_int = −Σ_links β·P. In a rigid-rotation shift β = ω×(r − X₀), this is exactly −ω·L, with L = Σ (r_link − X₀)×P_link (53's `angmom` about X₀). So a rotor's energy splitting gives Ω·L with Ω = −½∇×β, and for a ring of current it is the loop integral of β along its own path (Stokes). The light clock on a loop reads 2∮β/N⁴ (55, synchronised). So light and an orbiting rotor measure the same thing by construction of the drift, up to the size and shape of what each samples.
2. **Spin should not precess fully.** The coupling −β·P uses the lattice momentum density Im(ψ†ψ′), which is the canonical density. GR couples g₀ᵢ to the symmetric (Belinfante) density, which adds ½∇×s (s = ψ†Σψ, the spin density). With P = c_B·½∇×s for a packet at rest, summation by parts gives E_int = c_B Σ s·Ω, so the spin splitting is E_up − E_down = 2c_B·Ω, against GR's Ω (Dirac spin ½). A plan-stage fit gave c_B ≈ 0.10 at σ = 3 and ≈ 0.045 at σ = 5, a lattice remnant that shrinks with width. So the prediction is that **spin feels about 2c_B of the GR frame dragging, falling with packet width**: an equivalence-principle gap between spin and orbit unless the coupling gains the spin term.
3. **The light-drift readout.** For 56's covariant rule with uniform β, the right-moving branch has ω(k) = 2 sin(k/2) − β sin k, so v_g = cos(k/2) − β cos k exactly, and the exact right-moving eigenvector is ŵ = −iΩ₀(k)û with Ω₀ = 2 sin(k/2). An energy centroid moves at the mode-energy-weighted mean of v_g, and the mode energy is about Ω₀²|û|², not |û|². At k₀ = 0.25, width 16 this weighting moves ⟨cos k⟩ by about 0.4%, the size of 56's miss.

## Parts

**57A. Light readout (exact, Fourier space).** 1D periodic chain L = 512, N = 1, uniform β = 10⁻³, T = 100, widths 8, 16 and 32, k₀ = 0.25 and 0.5.
- Launch (a): the exact one-way eigenpacket (û = Gaussian on k > 0 only, ŵ = −iΩ₀û, real parts).
- Launch (b): 56's launch (u = env·cos k₀x, u_t = Ω₀(k₀)·env·sin k₀x).
- Propagate exactly by the 2×2 mode matrix in Fourier space (no time-stepping). Readouts with and without β, over βT:
  - the energy centroid (56's density ½w² + ¼(gradient terms) + βw∂u);
  - the envelope centroid (|analytic signal|²).
- Predictors: −β⟨cos k⟩ weighted by |û|² and by mode energy ε_k. One RK4 cross-check through 56's `light_drift` for (b) at width 16.

**57B. Frame dragging by the rotating ring (static, frozen fields).** 53's n = 129 spin-balanced ring (m = 0.6, radius 10, ℓ = ±2, u = 0.01), its static N and the shift `shift_static` (momentum rule, λ = 4). Everything is first order in β and uses the ℓ = +2 and −2 fields (odd part).
- **Point comparator:** the lattice curl of β at each probe centre (plaquette sums, averaged over the 8 cells around the site), Ω_curl = −½∇×β, and the GR dipole formula Ω_LT for reference.
- **Light.** Synchronised wave clock (55's `exgc55_posthoc2.wave`, sg = −1) on a square loop of half-side 8 (L = 64), width 4, k₀ = 0.5, centred at each probe position in the plane normal to z. Ω_light = −Δt/(4A), A = 16². The readout 2∮β/N⁴ is reported beside it.
- **Orbiting matter.** A spin-balanced probe ring (m = 0.6, radius 5, width 1.5, ℓ_p = ±1, axis z), built in flat N and rolled to the probe positions. Ω_orb = (E_int(+ℓ_p) − E_int(−ℓ_p)) / (L_z(+) − L_z(−)), with E_int = −Σβ·P_probe and L about the probe centre.
- **Torque.** The same probe with axis x at the on-axis positions: dL/dt from i⟨[H_β, L]⟩ (a symmetric finite difference of angmom along −iH_βψ). Then Ω_torque = (dL_y/dt)/L_x. GR: dL/dt = Ω×L.
- **Spin.** A rest Dirac packet (packet_k, m = 0.6), spin up and down along z, σ = 3 and 5, at the on-axis positions. Ω_spin = E_up − E_down. In the same run, c_B is fitted (P_link against ½(∇×s)_link, least squares on the rest packet).
- **Positions:** on axis (0, 0, d) and in the plane (d, 0, 0), d = 24 and 32. The dipole gives Ω_z(axis)/Ω_z(plane) = −2.
- **Pure-gradient control:** β = ∇χ with χ a Gaussian bump (width 8, same peak |β| as the ring field at d = 24) centred near each probe: Ω_orb, Ω_spin and the light delay should all vanish.

## Predictions (scored)

- **P1 (light drift closed).** With the exact one-way launch (a), the energy-centroid drift / (−β⟨cos k⟩_ε) = 1 within 10⁻³ for all six cases. With the |û|² predictor, the envelope centroid drift matches within 10⁻³. With 56's launch (b) at width 16, the energy-centroid ratio to the |û|² predictor reproduces 56's 0.9942 (k₀ = 0.25) and 0.9972 (k₀ = 0.5) within 5 × 10⁻⁴.
- **P2 (one frame dragging for light and orbiting matter).**
  - Ω_orb / Ω_curl within 10% at all four positions;
  - Ω_light / Ω_curl within 10% at all four positions;
  - Ω_orb / Ω_light within 10% at all four positions;
  - every Ω flips sign with ℓ (odd part / even part > 10⁶).
- **P3 (dipole shape).** At d = 32, Ω_orb(axis)/Ω_orb(plane) = −2 within 15%.
- **P4 (torque = precession).** Ω_torque / Ω_orb within 10% at the two on-axis positions.
- **P5 (spin feels little of it).** Ω_spin / Ω_curl below 0.3 at both widths and smaller at σ = 5 than at σ = 3, and within 25% of 2c_B fitted on the same packet. (GR predicts 1. A pass means P-coupled matter has a spin–orbit gap in frame dragging.)
- **P6 (gradient control).** In β = ∇χ, |Ω_orb| and |Ω_spin| below 1% of their ring values at d = 24, and the light delay below 1% of the ring delay at d = 24.

## Checks

- **C1.** 57A with β = 0 gives zero drift (to 10⁻¹²). The Fourier propagation of launch (b) matches 56's RK4 run at width 16, k₀ = 0.25 within 10⁻⁴ relative.
- **C2.** Rigid-rotation identity: in β = ω×(r − X₀), −Σβ·P_probe = −ω·L to 10⁻¹².
- **C3.** Wave-clock instrument check on L = 64, width 4: uniform β gives wave delay / 2Σβ within 1% (55's L = 256 value: 1.0036).
- **C4.** The ring's q and N for ±ℓ agree (53: identical), and the ℓ = 0 shift is zero.

## What the outcomes mean

- **P1 passes.** 56's light miss was the readout, not the physics: light and matter share the drift exactly.
- **P2–P4 pass.** The shift's curl drags light loops and orbiting matter by the same local Ω, with GR's dipole shape: frame dragging in this model is the substrate's local rotation, read the same way by everything that moves through it.
- **P5 passes.** Spin is left out. The momentum coupling is canonical, so an elementary spin feels only the lattice remnant 2c_B of the drag. Matching GR's spin precession (and the equality of spin and orbital gyroscopes) would need the symmetric momentum density in both the coupling and the source. That is a concrete requirement for 58's event rule.
- **P2 fails for matter but not light (or the reverse).** The drift is not one substrate rotation, and the 56 reading breaks down for curls.
- **P5 fails with spin near 1.** The lattice P carries the spin term after all, and the plan-stage fit was wrong.

## Outputs

`src/exgc57.py`, `src/exgc57_summary.py`, `src/exgc57_recompute.py`, `results/EX-GC57/`, note, LOG row and HANDOFF section.
