# Review of EX-GC22 (ChatGPT) and evaluation of the proposed EX-GC23

Reviewer: Claude, 2026-11-20, at the PI's request. EX-GC22 is EXPLORATORY. This is a same-programme review, not an
independent external replication. Branch reviewed: `ex-gc22-chatgpt` at 028f637, based on df7e73a.

## Verdict
**EX-GC22 is sound, and its main finding is correct and important.**

The coasting loss I reported in EX-GC19/20 was mostly a defect in my flow model, not a physical drag. My velocity field
V was not carried along with the moving particle. Because the particle's own inward push builds a small radial V inside
the core, the consumption-weighted "core velocity" drops as the pattern moves through that fixed radial field, even
with zero net force. My EX-GC20 "accretion drag" hypothesis is therefore not supported. Neither is the EX-GC20 reading
of the push-off control as evidence for a real backward push: switching the push off also removes the radial field, so
that control could not tell the two explanations apart. I have added errata to the EX-GC19 and EX-GC20/21 notes.

I recommend merging EX-GC22 into main (done with this review) and running EX-GC23, with the additions below.

## What I checked
1. **Regression tests.** `python tests/test_regression.py` passes on the branch.
2. **Prescribed-sink arm reproduces exactly.** Rerunning `src/exgc22_drag.py fixed` on this container gives a maximum
   relative difference of 0 against the committed `fixed.json`.
3. **Analytic drag formula, derived independently.** For a Gaussian sink, angular averaging gives
   ∫ e^(−k²σ²) k_z²/k⁴ d³k/(2π)³ = 1/(12π^(3/2)σ), so γ = M/(120π^(3/2)σD), the note's formula. The finite periodic boxes
   sit below this infinite-box value (0.0087 at side 32 and 0.0093 at side 64 for K = 0.2, D = 1, with G = 1 in the
   force, against 0.00995), as expected from the missing low modes. The full J/G force is larger (0.0114) because G is
   depleted at the core.
4. **Full-particle arm and start-of-kick balance reproduce.** Rerunning `full --workers 2` (3.5 min wall), then
   `exgc22_initial_balance.py` and `exgc22_summary.py`, reproduces every numeric field of `summary.json` to 3 × 10⁻¹²
   relative. It also reproduces `initial_balance.json` to all printed digits. This uses the same code, so it checks
   reproducibility, not independence. The committed records were kept unchanged.
5. **Independent consistency with my own runs.** If the pattern moves at v through an inward radial field V ≈ −αr, the
   reweighting term gives a loss rate of about α. In EX-GC19/20 the recorded flow spread was about 0.018 over a core of
   radius about 1, so α ≈ 0.018. That matches the loss rates I measured, 0.017–0.023. My own logs thus support
   ChatGPT's mechanism, without using their code.
6. **Algebra of the balance.** d/dt(∫sV/∫s) = ∫s ∂ₜV/∫s + ∫(∂ₜs)(V − V_core)/∫s is exact. The transport check
   −v₀∂_zV is a diagnostic only, as the note says.

## Points to note (none blocks merging)
- **Physical motion is still unresolved.** The W centroid moves at only 0.38–0.58 of the flow speed at h = 0.5, even
  for the weakest source. Grid pinning at h = 0.5 was seen before (EX-GC14/15). EX-GC23 must measure centroid motion at
  h = 0.25.
- **Commit hash mismatch.** The note cites the pre-run plan as local commit 890c21e, but in the pushed history the plan
  is f64aa4e, the first commit on the branch, ahead of all results. The pre-run order holds; the hash changed in
  transfer.
- **Out-of-date sentence.** The note says no remote write succeeded. The branch is now on GitHub, so that sentence is
  out of date. No action needed.
- **Assumptions still in place.** The weak-source results assume inertia proportional to absorption (assigned, not
  derived), and the EX-GC17 head-on absorption drag is still switched off. The note states both.
- **Optional C++ path.** The C++ stencil is optional and needs a C++ compiler. On Richard's Windows PC it would need
  MSVC or MinGW, or the NumPy path.

## Evaluation of EX-GC23
I agree with the proposal. It is the right next step, and it should come before any Lorentz-invariant redesign.
Suggested refinements:
1. **Primary observable: the W (or s-weighted) centroid trajectory**, with velocity-field averages as secondary.
   Coasting means the particle keeps moving, not that a field average stays constant.
2. **Material derivative first, then a conservative form.** Add (V·∇)V, rerun the strong wave-like kick (K = 0.2,
   D_eff = 80) with and without it at h = 0.25, and run a weak case (K = 0.01).
   - A momentum density P = ρV needs ρ to be carried with the pattern. Here ρ ∝ s is set by the chemistry, so P is not
     conserved whenever s changes. Either carry ρ as its own advected field, or report the momentum ledger
     (recoil, cohesion, creation and loss of ρ) explicitly.
   - The cohesion term as written does not conserve ∫ρV. Replace it with a momentum-conserving exchange, for example a
     symmetric pairwise or viscous-flux form, or record its net contribution.
3. **Check survival.** The inward self-push still exists, so the core must survive with transport switched on (M, W
   peak).
4. **Flux response.** Repeat the uniform-flux trials to confirm the fall acceleration (2J/(Gτ)) and the near-type
   independence survive transport (K = 0.2 against 0.4, as in EX-GC19).
5. **Walls.** One larger box (±9) for the wave case, to separate wall reflections (the EX-GC20 ringing) from the
   physics.
6. **Then the bare absorption drag (EX-GC17).** This is the real remaining obstacle in a medium with a rest frame. Treat
   it as a separate study.
7. **Budget.** About 2 core-hours with the C++ kernel at h = 0.25, if limited to roughly six full runs of 20 time units.
   Ownership is the PI's decision. Either AI can do it; ChatGPT has the working kernel.
