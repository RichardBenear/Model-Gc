# EX-GC02 (EXPLORATORY): inertia in Model Gc through a slow core species

Run by Claude, 2026-11-18, under the PI's go-ahead ("yes, go ahead"). Well under 1 core-hour. Radial 3D, linear dipole
(l = 1) theory plus a direct linear time-response check. Exploratory: this supports no claim on its own.

## Result
- **Inertia appears, and it is real.** We gave the core species W a slower response (time constant τ_W; the resting
  particle is unchanged). The particle's drag γ then falls linearly and an effective mass m appears. At τ_W = 25 the
  drag is 0.12 of its baseline, the mass is m = 34 (in units where the baseline drag is 1), and the velocity takes
  about 290 time units to build up.
- **The direct test confirms it.** We switched on a G gradient and integrated the linear dipole equations. The
  displacement follows m·X'' + γ·X' = f to within 3–4%, against 8–23% for the inertia-free law γ·X' = f
  (the figure's third panel). At τ_W = 1 both agree with the simulation to 0.2%.
- **Drift bifurcation.** The drag reaches zero at τ_W ≈ 28.3 for the reference particle. Beyond that a resting
  particle is unstable and starts moving on its own. No other instability was found up to τ_W = 27 (l = 0, 1, 2,
  shifts along the real and imaginary axes).
- **The cost: free fall becomes non-universal.** The threshold depends on how much each particle's core contributes
  to its drag (the W weight), so it differs by particle type: τ_W ≈ 28.3 (reference), 40.1 (wider core, D_W 0.2),
  17.5 (doubled consumption, k 0.4). Near threshold the mobilities therefore spread apart. At τ_W = 25 they are
  −0.214, −0.066 and (k 0.4) already self-propelled, against −0.026, −0.025 and −0.046 at τ_W = 1. The EX-GC01
  equivalence-principle problem gets much worse, not better.

## Why (the structure)
- **Drag.** For time constants τ_s, the drag is γ(τ) = Σ_s τ_s·c_s, where c_s = ⟨ψ†_s, U'_s⟩ is fixed by the resting
  particle. For the reference particle c = (X 1.130, Y −0.096, W −0.037, G 0.003), normalised so γ = 1 at τ = 1.
  Only species with negative weight can lower the drag.
- **Why W and not Y.** Slowing Y would reach zero drag at τ_Y ≈ 11.4, but the ether's own X–Y oscillation sets in
  first, at τ_Y = A²/(B − 1) = 8.64. W is zero in the ether, so slowing it is safe there.
- **Mass.** The effective mass is the next term in the expansion (Jordan chain): m = −⟨ψ†, T w⟩ with
  J₁w = −(T − α)U'. It is carried almost entirely by W (35.7 of 34.2 at τ_W = 25). It agrees to 0.6% between grid
  spacings 0.05 and 0.025 and between domain sizes 60 and 120.
- **Mobility.** μ = −⟨ψ†, P⟩/γ, so as γ → 0 the mobility grows in proportion to 1/γ. That is how the
  type-dependence of the threshold becomes type-dependence of free fall.

## What it means for the physics mapping
- **Inertia.** A particle can acquire inertia from its own internal response time. There is a regime where it is
  nearly drag-free and keeps its velocity for a long time: the closest this model has come to Newton's first law.
- **Self-propulsion is not Galilean.** Beyond the threshold, motion is self-propelled at a speed fixed by the
  chemistry, so there is no freedom to move at any speed.
- **What Newtonian behaviour would need.** Gravitational acceleration a_i ∝ M_j requires the inertial mass m_i to be
  proportional to the gravitational mass M_i, with the drag negligible, for *every* particle type. Here m and γ are
  set by core structure that differs between types, so this regime gives inertia at the price of universality.
- **Revised design target.** A coupling in which the same quantity sets both a particle's G consumption and its
  inertia, so that the ratio is fixed by structure. The natural candidate is to let G itself (or its consumption)
  supply the slow mode, rather than an independent core species.

## Method and records
- `src/exgc02_inertia.py`: weights, Jordan-chain mass, generalized eigenvalues (J₁ψ = λTψ), stability, and the time
  response.
- `src/exgc02_figure.py` produces `paper/figures/EX-GC02_inertia.png`.
- `results.json` holds the weights, scans for the three types, stability (basic and extended), and response curves.
- The l = 1 grid "pinning" (translation eigenvalue −1.8e-4 at dr = 0.05, −4.4e-5 at 0.025) is a discretisation
  effect. It scales as dr², and it is included as a weak spring in the response comparison.

## Limits
- Linear theory only (small speeds). No nonlinear moving-particle or pair run.
- One slowed species (W). Combined τ_Y and τ_W and other routes were not scanned.
- The mass is extracted near threshold by the leading-order reduction. Away from threshold it is a small correction.
