# EX-GC52 — a moving source in the 51 model: does a single coincidence-carried scalar pick out a preferred frame?

2026-11-27 | EXPLORATORY | Claude | Plan pushed before any run. PI go-ahead 2026-11-27 19:54 UTC ("for now, let's do the moving source run first"). Cap 2 core-h (expected about 1.3). Our construction; not LaViolette's Model G and not Chain Theory.

## Why

51 passed the static tests. A theory built on a single field can match all of them, and the theories that did so historically failed on moving matter. This study asks whether gravity in the 51 model depends only on the relative motion of bodies (no preferred frame), as in GR. The alternative is that the lattice's own rest frame shows up in the gravitational force. The PI expects something beyond the single scalar to be needed, for example torsion in the field lattice or clock, or a non-symmetric rewriting rule. 52 measures how large the gap is, and in which part of the force it sits.

## The model

The model is exactly as declared in 51 and is not changed here.
- **Matter.** Wilson–Dirac matter with bond weights N_vN_w and site mass weight N_v. The charge is q = ∂H/∂N.
- **Field.** (1/N)∂_t[(1/N)∂_t f] = Σ N_xN_y(f_y − f_x) − κq, with f = −1/(2N²).
- **Units.** The lattice rest frame is the frame of the coordinate t.

## Continuum prediction (worked before any run)

Set ψ = ln N. For a weak field and r = 0, three things hold.

1. **The field is a wave.** Linearised, the rule is ψ_tt − ∇²ψ = −κq.
2. **Matter follows geodesics.** A point particle has H = √(N⁴p² + N²m²). Its Lagrangian is −mN√(1 − v²/N⁴), which is exactly geodesic motion in g₀₀ = −N², g_ij = N⁻²δ_ij.
3. **A moving body's charge grows with speed.** For a body with velocity V, the charge is Q(V) = ∂H/∂N = E(1 + V²) = γm(1 + V²).

**The field of a uniformly moving source** (point charge, steady state) is

ψ_V(X, ρ) = −κQ(V)γ / (4π√(γ²X² + ρ²)), where X = x − Vt.

This is the Lorentz-contracted field, with charge Q(V). The model's g₀₀ is −e^{2ψ}. At linear order it equals the boosted Schwarzschild g₀₀ exactly, because GR's h₀₀ for a moving mass also carries γ²(1 + V²).

The difference from GR is that the model has no g₀ᵢ, the gravitomagnetic part. That part cannot be removed by a change of coordinates. It cancels the extra attraction between bodies that move together, just as the magnetic force between parallel moving charges cancels part of their electric force.

**Prediction for a co-moving pair.** Two bodies move side by side with the same velocity V along x, at transverse separation ρ.
- In the model, the lab-frame transverse force on the probe is ṗ_y = −(∂H/∂N)∂_yψ ∝ γ³(1 + V²)².
- Lorentz invariance requires ṗ_y(V) = ṗ_y(0)/γ.
- The model therefore breaks it by the factor R(V) = γ ṗ_y(V)/ṗ_y(0) = γ⁴(1 + V²)² ≈ 1 + 4V². This is 1.174 at V = 0.2.
- In PPN language, g₀ᵢ = 0 means α₁ = −8 and α₂ = −1. Observation bounds these at about 10⁻⁴ (α₁) and below 10⁻⁷ (α₂).

**Probe at rest beside a moving source.**
- The transverse force at closest approach is the same as GR's at linear order.
- The longitudinal force differs from GR's through the ∂_t g₀ₓ term. Its relative size is about 4V².

So the expected outcome is a preferred-frame effect of order one in the V² terms. 52 checks whether the lattice model really does this, and measures the lattice corrections.

## Parts

All parts use r = 0 (naive fermion, with packets at small k so the doublers are unoccupied) as the primary case, because the Wilson r·m term distorts the dispersion at order one (48). r = 1 appears in 52A only.

**52A. Charge of a moving packet (k-space, cheap).**
- **Packet.** A positive-energy packet with mean momentum k₀x̂ chosen so that its group velocity is V = 0, 0.1, 0.2 and 0.3, with m = 0.6 and σ = 5, in uniform N = 1, L = 96.
- **Measure.** Q = Σ q_x and E.
- **Compare.** The lattice-exact prediction ⟨∂E/∂N⟩ averaged over the packet's k distribution, and the continuum value γm(1 + V²).

**52B. Field of a moving source (steady comoving solution).**
- **Source.** The rigid 52A source density q(x − Vt) for each V.
- **Equation.** Solve the full rule in its steady comoving form, with ∂_t → −V∂_x: (V/N)∂_x[(V/N)∂_x f] = Σ N_xN_y(f_y − f_x) − κq.
- **Method and size.** Picard iteration with a DST preconditioner (anisotropic Laplacian ∇² − V²∂_x²), Dirichlet faces set by the linear moving solution, n = 129.
- **Strength.** κ chosen so that U(ρ = 12) ≈ 0.01 at V = 0.
- **Readouts.**
  - The fitted contraction: the ratio of the half-width along x to that along ρ, as a function of V.
  - The amplitude at X = 0, compared with the predicted κQ(V)γ/(4πρ) over 12 ≤ ρ ≤ 36.

**52C. Forces from the 52B fields (semiclassical, lattice-exact dispersion).**
- **Co-moving probe.** m = 0.6 at transverse separation ρ = 16, with packet k distribution as in 52A, so its velocity is V.
  - ṗ_y = −⟨∂E/∂N⟩ ∂_yN, read from the 52B field at the probe and averaged over the probe's density.
  - R(V) = γṗ_y(V)/ṗ_y(0).
- **Probe at rest at ρ = 16.** Transverse and longitudinal forces as the source passes (X from −32 to +32), compared with GR's linearised boosted-Schwarzschild force.

**52D. Full dynamic check (one speed).**
- **Setup.** Two packets with m = 0.6 and σ = 5 (r = 0), transverse separation 14, both with group velocity V = 0.2 along x, n = 81, T = 16, dt = 0.1. The field is fully dynamic, starting from the 52B-type steady comoving solution of the frozen rigid pair at the same V.
- **Control.** The same with V = 0.
- **Measure.** The transverse canonical-momentum rates dP_y/dt of each packet, fitted as in 51D, and the dynamic R_dyn = γ (dP_y/dt)_V / (dP_y/dt)₀.

## Predictions (scored)

- **P1 (charge).** Q(V) equals the lattice-exact ⟨∂E/∂N⟩ to 1e-6 relative at every V (r = 0 and 1). For r = 0, Q/(γm(1 + V²)) lies within 2% of 1 at V ≤ 0.2.
- **P2 (moving field).** At V = 0.1, 0.2 and 0.3, the 52B amplitude at X = 0 matches κQ(V)γ/(4πρ) within 2% over 12 ≤ ρ ≤ 36. The fitted contraction ratio is within 2% of 1/γ. This would show no preferred frame in g₀₀.
- **P3 (preferred frame, semiclassical).** R(V) from 52C equals the model's lattice-exact prediction within 1% at V = 0.1, 0.2 and 0.3. It also exceeds 1 by more than 50% of 4V² at V = 0.2 and 0.3, which is a preferred-frame effect at order V². Lorentz invariance (GR) would give R = 1.
- **P4 (dynamic).** R_dyn at V = 0.2 equals the 52C semiclassical R(0.2) for the same packets within 3%.
- **P5 (probe at rest).** The transverse impulse on a probe at rest, integrated over the source's pass, matches GR's linearised value within 3% at V = 0.2. The longitudinal force is reported against GR's but not scored, because I have not fixed the gauge of that comparison.

## Checks (separate paths)

- **C1.** The 52B solver at V = 0 reproduces the 51 static solver (`exgc51.static` with the same frozen source) to 1e-10.
- **C2.** At small κ, the 52B linear part matches an independent FFT solution of the anisotropic Poisson equation (periodic box, image-corrected) to 1e-3 over the fit window.
- **C3.** 52D at dt/2 for the V = 0.2 case: R_dyn changes by less than 1%.

## What the outcomes mean

- **P2 and P3 pass, as predicted.**
  - The single coincidence-carried scalar gives the right field shape for a moving mass: the g₀₀ part is Lorentz-contracted correctly.
  - It also gives a preferred-frame force between co-moving bodies of order V², with α₁ ≈ −8 and α₂ ≈ −1, which observation excludes by many orders of magnitude.
  - So the 51 model is a static-gravity mechanism only. Moving matter needs a gravitomagnetic part (g₀ᵢ) from somewhere.
  - The size and form of the missing term (−4V² in R, the ∂_t g₀ₓ force) are then the target for the extra ingredient the PI proposes: torsion of the lattice or clock, or a non-symmetric rewriting. The next study would ask whether a coincidence rule with a handedness, or a shift of the firing pattern carried by moving matter, supplies g₀ᵢ = −4Vψ.
- **R(V) ≈ 1.** Something in the lattice rule supplies the cancellation that the continuum analysis misses. That would be a real surprise, and I would trace it.
- **P2 fails.** The comoving field is not the Lorentz-contracted one. Then the coincidence rule itself has a preferred frame already in g₀₀, which is worse.

## Outputs

`src/exgc52.py`, `src/exgc52_recompute.py`, `results/EX-GC52/`, note, LOG row and HANDOFF section.

## Amendment A1 (2026-11-27, written with the code, before any scored run)

A smoke test of the code at small size (n = 48 and 33, not a registered run) showed the following. At σ = 5 the packets carry enough internal momentum that Q(0) = 1.22 m (r = 0). So the continuum clause of P1, Q/(γm(1 + V²)) within 2% of 1, measures the packet's spread, not the moving-source physics. The changes below were all made before the scored runs.

1. **P1.** Only the lattice-exact clause is scored: Q equals ⟨∂E/∂N⟩ to 1e-6. Q(V)/Q(0) against γ(1 + V²) is reported, not scored.
2. **P2 contraction readout.** The contraction is read as an amplitude ratio. For a point source, ψ·ρ at X = 0 is −κQγ/4π and ψ·|X| on the axis is −κQ/4π, so their ratio is γ. The score is that ratio over γ, within 2%, averaged over 12 ≤ ρ, |X| ≤ 36 (both signs of X).
3. **52D size.** n = 65 instead of 81, to keep within the cap. The packets sit at y = c ± 7, about 25 sites from the faces.
4. **52D self-field.** As in 51D, a single-packet run with the same V is subtracted. Only packet 1's force is scored.
5. **C2 method.** The linear part (κ = 10⁻⁶) is checked against an independent conjugate-gradient solve of the same discrete anisotropic problem (scipy sparse), to 1e-8, instead of the periodic FFT with images.
6. **C3.** The check is F₁(V = 0.2) at dt/2, within 1%, rather than R_dyn, to keep within the cap.
