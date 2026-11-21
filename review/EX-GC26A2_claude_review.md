# EX-GC26A2 — Claude review (EXPLORATORY)

Reviewer: Claude, 2026-11-20. Reviewed commits 30ac594 (plan) and 2faf866 (results), by ChatGPT. Plain massless Skyrme
model; the numerical construction is ours.

## Verdict

Sound, and it answers the questions 26A left open. The knot's topological degree is exactly conserved (the
non-integer values in 26A were the estimator, not the particle), the lopsided start reaches the same state, and the
momentum leak is mostly a box effect rather than grid pinning. I agree with its dispositions (PASS / CONDITIONAL /
PARTIAL) and with proceeding to 26B as a gated exploratory study under ChatGPT's five conditions, with the additions
below.

## What I checked

1. **Validation rerun.** `python src/exgc26a2_qualify.py validate` reproduces `validation.json` exactly (energy 3.9e-16,
   gradient 2.8e-17, C++ vs NumPy acceleration 2.5e-14; checkerboard energy positive for the higher-order stencil).
2. **Boost with my original NumPy solver only** (no C++ kernel), h = 0.2, box ±5, axis, the same initial velocity
   −0.1 ∂ₓφ, T = 1.5, dt = 0.01, from my own saved 26A state: max |ΔP|/E₀ = 4.1858e-4 (theirs 4.1851e-4), energy drift
   4e-14, initial transverse momentum −0.000842 per component (theirs −0.000842). Same code author for the solver, so
   this checks their kernel and bookkeeping, not the physics.
3. **Geometric degree of the 26A state** with their tetrahedra routine: −1 for all three regular values. The hemisphere
   certificate held; the many degenerate tetrahedra are vacuum-region cells, as expected.

## Comments

- **Geometric degree is the important new result.** It shows the particle number is exactly conserved in every state
  tested, including after the boosts. 26A's C1 degree criterion was set on the wrong estimator; it stays "not met" as
  registered, but the physical question is now answered.
- **Momentum leak is a box effect.** Halving it by enlarging the box (L = 5 → 7), with no change from halving dt and
  little from refining h, points to the clamped wall acting on the massless Skyrmion's long tail (f ∝ 1/r²), like an
  image force. *(My reading, not tested.)*
  - Option for 26B: add the standard pion-mass term (m = 1 in Skyrme units). The tail then decays exponentially and the
    wall effect becomes negligible at L ≈ 6. It is a small, standard change to the model, but it needs a new radial
    reference and rest state.
  - Otherwise use L ≥ 7 to 9 and include a boundary momentum-flux term, as ChatGPT asks.
- **Pinning.** Agree with "conditional". The Fourier-taper test is post hoc and changes the tail, so it bounds the
  pinning of a slightly different field. The dynamic momentum ledger in 26B is the measurement that matters.
- **Higher-order stencil.** Promising (Derrick ratio 1.001 at h = 0.15) but not converged and without dynamics. Agree
  to keep it separate.
- **Protocol.** The plan reached main 34 seconds before the results (both pushed by the PI). The note cites a local
  prospective commit (a0b5c58) that is not on main, so the pre-run order cannot be checked from the repository. The
  execution amendments are clearly labelled, and the taper test is labelled post hoc. Compute 0.37 core-hours, within
  the cap.

## On 26B

I agree with ChatGPT's five conditions: a larger box, measured (not imposed) velocity, axial vs diagonal mass, the
boundary momentum flux, and convergence toward E₀/c². Additions:

1. **Framing.** For any Lorentz-invariant action, M_i = E₀/c² in the continuum is a theorem (von Laue). So 26B measures
   how far the lattice and box approximation departs from it, not whether inertia emerges. Its real purpose is error
   bars for 26E (gravitational vs inertial mass). I'd target about 1% and not spend PC hours beyond that.
2. **Early indication.** From the existing boost: initial P/v_imposed = 0.994 E₀, and using the measured centroid speed
   (about 0.097) gives P/v ≈ 1.03 E₀. That is already within a few percent, before any refinement. *(Rough, not a
   measurement.)*
3. **Use exact Lorentz-boosted initial data** at several speeds (e.g. 0.1, 0.3, 0.5) to test E = γM and P = γMv. The
   small-velocity translation used in A2 is only valid as v → 0.
4. **Track the energy-density centroid** for velocity, and close the momentum ledger with the stress-tensor flux through
   a surface around the particle.

**Who runs it.** Either works. ChatGPT has the faster C++ kernel and the A2 diagnostics in hand, which argues for it
running 26B with me reviewing.
