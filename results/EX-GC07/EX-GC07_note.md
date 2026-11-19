# EX-GC07 (EXPLORATORY): 3D time simulation of the clocked particle

Claude, 2026-11-18, at the PI's request ("go ahead with the 3D time simulation"). Exploratory track. Gc3 (Gc1 plus a
slow core inhibitor, EX-GC06) is our construction, not LaViolette's.

## Answer
**The breathing does not settle into a steady pulse. It grows until the core collapses and the particle dies.** The
core never started circling. So in Gc3, as set up, the internal clock destroys the particle instead of running
it, and the 2:1 test cannot be done yet.

| Run | Linear prediction (EX-GC06) | 3D result |
|---|---|---|
| κ = 1, τ_V = 4 | breathing, growth 0.0512, ω 0.308 | grows at 0.0525, ω 0.314; collapses at t ≈ 182 |
| κ = 0.9, τ_V = 7 | breathing, growth 0.0228, ω 0.292 | grows at 0.0236, ω 0.286; reaches ±30% in M, then collapses at t ≈ 425 |
| κ = 0.9, τ_V = 6 | just past onset | still growing (0.0078) at t = 800, amplitude 0.8% of M; not settled |
| κ = 0.9, τ_V = 50 | breathing (0.10) + slow l = 1 | collapses at t ≈ 82 |
| κ = 1, τ_V = 20 | real l = 0 growth (0.30) + l = 1 oscillation (0.006) | collapses at t ≈ 48 |
| κ = 0.9, τ_V = 5.5, 50× larger kick | stable | kick decays (amplitude in M 1.3e-3 → 7.5e-4 over 300) |

- **Validation.** Early growth rates and frequencies agree with the EX-GC06 eigenvalues to 2–4%, and the unperturbed
  particle stays steady to 1e-15. This checks the solver, not the physics.
- **No saturation seen.** In every run that got past about 30% of M the particle died, and the growth stayed
  exponential right up to collapse. That looks like a subcritical (or very high) saturation. A stable small pulse
  might still exist just past onset, but the τ_V = 6 run had not settled by t = 800 (an e-folding time of about 130).
- **No circling.** The l = 1 content stayed at the seeded level (≤ 1e-3) and the core centre moved less than 1e-3.
  Where l = 1 is linearly unstable (κ = 1, τ_V ≥ 19), the round breathing or collapse mode is much faster and kills
  the particle first.

## What it means for the 720° route
Route 1 needs a particle that both rotates and keeps a clock. Gc3 gives a clock only in the linear sense: the
pulsation it starts does not stop growing. The simplest fix would be a brake that saturates, so the oscillation caps
itself, for example a V term like −κVW/(1 + V/V_s). That is another new model choice. Each step in this direction
adds a hand-tuned ingredient to make the particle do something it does not do by itself. That is worth weighing
against routes 2 and 3, where the 720° property comes from structure (an orientation field, or ribbons in the graph)
rather than from tuned dynamics.

## Method
- `src/exgc07_sim3d.py`: fields = steady background + δ, with δ expanded in complex spherical harmonics up to
  l = 4 (all m) on the radial finite-volume grid of the earlier radial studies (R = 30, dr = 0.05). Nonlinear terms
  are evaluated on a Gauss–Legendre (8) × Fourier (16) angular grid. Diffusion is implicit per (species, l), using
  the same radial operators as the eigenvalue code; reactions are explicit (IMEX Euler, dt = 0.004). Fixed points are
  exact steady states. Seeds: random l ≤ 2 perturbations of W with relative size 1e-3 (5e-2 for the "kick" run),
  shaped like the core.
- Diagnostics every 0.5: consumption M, W-weighted centre, and l-resolved W amplitudes. Time series are in
  `ts_*.json.gz`. The final states (about 0.9 MB each) are kept out of git, with SHA-256 in
  `final_states_sha256.txt`; copies are in the project files (`model-gc/EX-GC07_states/`).
- `src/exgc07_summary.py` → `summary.json`, `paper/figures/EX-GC07_breathing.png`.
- Platform: Linux cloud container, Python 3.13.18, numpy 2.5.3, scipy 1.18.1. Compute: about 1.8 core-hours, within
  the 2 core-hour cap.

## Deviations and limits
1. **Not a full 3D grid.** The angular truncation at l = 4 is accurate for small and moderate deformations. Large
   off-centre motion of the core would not be represented, but none occurred. The collapse is radial (l = 0), so it
   is well represented.
2. **Two runs were added:** τ_V = 6 (nearer onset) and a 50× kick at τ_V = 5.5 (a test for sudden collapse below
   onset). Dead runs (κ = 1, τ_V = 4 and κ = 0.9, τ_V = 50) were stopped early to save compute.
3. **The outer boundary is at R = 30**, not 120. That shifts M by about 0.2%.
4. **The τ_V = 6 run is unresolved.** Settling at small amplitude (supercritical) is not ruled out. A longer run, or
   a weakly nonlinear calculation, would settle it.
