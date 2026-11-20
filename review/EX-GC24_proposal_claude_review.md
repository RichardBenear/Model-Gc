# Review of the EX-GC24 proposal (ChatGPT): does the G medium generate inertial mass?

Reviewer: Claude, 2026-11-20. Reviews `spec/EX-GC24_proposal.md` at c2b3ce6. EXPLORATORY proposal; nothing has been
run.

## Verdict
**The question is the right one, and the design is careful.** Stage A, the prescribed-sink impedance, is the cheap
and correct first step. The S² versus S discriminator is the key insight.

However, the proposal's own linear-response formula already answers most of Stage A in closed form, and the answer is
discouraging for medium-generated inertia. **I suggest shrinking Stage A to a quick spectral confirmation**, then
spending the remaining effort on the comparison with head-on drag, or on the no-preferred-frame route.

## What the proposal's formula already implies
Expanding the bracket for small ω gives F̂/ξ = (2/c²) Σ |ŝ|² k_z²/k⁴ · (ω² + iω/τ_J). That means:
- **m_add = (2/c²) ∫ |ŝ|² k_z²/k⁴, and b_wake = m_add/τ_J.**
- **In the weak-damping limit, mass/drag from the medium equals τ_J exactly**, independent of the source. The
  self-wake "coasting time" is just the medium's relaxation time.
- **Spot check.** On a 64³, h = 0.5 periodic box (Gaussian σ = 0.8, S = 5.32, c = 2):
  - τ_J = 20 gives m = 0.2308, b = 0.01158 and m/b = 19.93, at ω = 0.01 and 0.003.
  - τ_J = 2 gives m/b = 1.31, below 2 because the damping is no longer weak relative to the box's lowest modes.
  - In both cases m > 0.
- **Gaussian, infinite box:** ∫ |ŝ|² k_z²/k⁴ d³k/(2π)³ = S²/(12π^(3/2)σ), so m_add = S²/(6π^(3/2)c²σ). That is 0.265
  for the values above, consistent with the finite-box 0.231.
- **Scaling with S.** As the proposal says, m_add ∝ S² at fixed shape, while the gravitational source and push scale
  as S. So acceleration under a given flux would go as 1/S, which breaks universal free fall badly, unless the shape
  scales to compensate (for example σ ∝ S). This is the classical "field-energy mass", like electromagnetic mass
  e²/(rc²). It is real but non-universal (Nordtvedt-type).

## Comparison with the head-on absorption drag (my estimate; to be checked in Stage A item 1)
- **Bare drag in this force normalization.** The EX-GC17 terminal speed is v_T = 2J/G. With the force proxy
  F = 2S·J/G − b_bare·v, that gives b_bare = S. This assumes the EX-GC17 lattice result carries over to the wave
  medium, which Stage A should verify.
- **Coasting time from medium mass alone** is then m_add/b_bare = S/(6π^(3/2)c²σ) ≈ 0.05 time units for the K = 0.2
  particle. For weaker particles it is smaller still, ∝ S.
- **Consequence.** In a medium with a rest frame, field-generated inertia is overwhelmed by absorption drag by orders of
  magnitude, unless c is tiny. A medium-added mass cannot make a consuming particle coast. Inertia would have to come
  from the particle's internal content, or the medium must have no rest frame.

## Suggested changes
1. **Stage A: confirm the analytic formula cheaply.** Use the spectral solver over ω, τ_J, c, σ, S and box size, plus
   one direct time-evolution check. Report m_add/b_wake → τ_J, the S² scaling, and the low-mode/box dependence. That
   is minutes of compute.
2. **Make the head-on drag comparison the decisive Stage A item.** Define b_bare in the same normalization (my estimate
   b_bare = S) and report m_add/b_bare. If it is ≪ 1 as estimated, stop the added-mass route there, with a clear
   negative result.
3. **Stage B (momentum ledger) only if item 2 is not decisive.** Otherwise move to the alternative:
4. **A no-preferred-frame (Lorentz-invariant) G field as EX-GC25 candidate.** In a relativistic field theory a soliton's
   inertia is its total energy/c², automatically, with no drag for uniform motion. If gravity couples to the same
   energy, universality follows naturally. This is a bigger model change, and the PI's decision.
5. **Keep the h = 0.125 slip check** (EX-GC23's 2.5% slip). It is independent and cheap.
6. **Repo size.** `paper/figures/ModelGc_topological_particle.png` is 1.8 MB. That is fine if Richard wants it kept;
   otherwise a compressed copy would keep the repo small.

Orbits: agree, deferred. Budget: Stage A as revised needs well under 0.5 core-hours.
