# EX-GC34B — what does each part of the shared structure respond to? (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "go ahead with 34B". Cloud, cap 2 core-hours (expected
< 0.5). Second stage of EX-GC34 (`docs/EX-GC34_outline.md`). Our construction, not LaViolette's Model G.

## Question

In 34A-r every field moved through one shared structure: a lapse N (clock rate per site) and a spatial metric γ_ij
(carried by link weights). If gravity is that structure responding to matter, the source of each part is fixed by how
the matter action depends on it:

- **lapse source** s_N = −∂L/∂N;
- **spatial source** T^{ij} = 2 ∂L/∂γ_ij (the stress tensor).

GC33A's universality test asks whether a source equals κ × (energy) for every kind of matter: rest, binding, motion,
radiation.

## Analytic part (written up in the note, checked numerically below)

1. **Lapse.** Every field built on the shared structure has kinetic terms ∝ 1/N and all static terms (gradient, mass,
   and the knot's quartic Skyrme term) ∝ N. So at N = 1, s_N = kinetic + static = the energy density, exactly, for
   every field, at every binding and motion. The GC33A identity holds for the lapse by construction.
2. **Spatial metric.** The source is the stress tensor. For a static bound body its integral vanishes in the
   continuum (von Laue / Derrick). For radiation moving along x it is ∫T^{xx} = E. For a moving body it is
   ∫T^{xx} = E v².
3. **Consequence (prediction, not a test).** The far field of a static bound body therefore has a 1/r lapse part and
   no 1/r spatial part *if each part of the structure obeys its own local wave/Poisson equation* (the "simplest"
   action in the outline). That gives PPN γ = 0: light bends at the lapse-only value, half of GR (34A-r B2). GR's
   γ = 1 needs the spatial metric to be driven by the energy too, through the trace structure of the linearised
   Einstein (Fierz–Pauli) action. **So the shared-structure idea fixes the coupling (universal fall, equivalence
   principle) but not γ.** γ is set by the gravity field's own action, which Model Gc does not yet determine. That is
   the decision for 34C.

## Numerical checks

**Knots.** B = 1 and B = 2 rest states regenerated with the EX-GC26E1 static optimiser (h = 0.2, L = 7, maxiter
300). Hashes and energies are compared with the committed `results/EX-GC26E1/static_B*.json`. B = 1 is also done at
h = 0.15 (the 26E2 settings) for the lattice-offset scaling. Outputs go to `results/EX-GC34B/`, with
26E1/26E2 records untouched.

- Metric dependence uses the 26E1 orientation stencil: e = √γ [γ^{ij} D_ij + ½ γ^{ik}γ^{jl}(D_ij D_kl − D_ik D_jl)].
- Sources are central differences in γ_ij at γ = δ (step 1e-5).
- The lapse source is computed as −∂L/∂N with a static-plus-kinetic split. For a rest knot the kinetic part is 0, so
  s_N = E is a trivial check here. The boosted case is covered by the analytic identity.

**Packets.** Massless and massive (μ = 1) scalar packets from the 34A solver (flat, k = 3 along x, σ = 2, forward
initial data), at t = 0 and after T = 5 of evolution. Their sources come from the same numerical γ- and
N-derivatives of the 9-link lattice Lagrangian.

## Registered expectations and criteria

| item | expectation | criterion |
|---|---|---|
| lapse source, every case | = total energy H | relative difference ≤ 1e-6 (implementation check) |
| B1, B2 rest: trace of ∫T^{ij} | = E4 − E2 (lattice Derrick offset; algebraic identity) | ≤ 1e-6 relative to E |
| B1, B2 rest: each component of ∫T^{ij} | ≈ 0 up to the lattice offset | all \|∫T^{ij}\| ≤ 1.2 × \|E4 − E2\| |
| offset scaling, B1 h = 0.2 → 0.15 | → 0 in the continuum (∝ h² gives 0.5625) | ratio of \|E4 − E2\|/E in [0.40, 0.75] |
| binding: B2 vs 2 × B1 | lapse source tracks binding exactly; spatial monopole stays ≈ 0 | reported |
| massless packet | ∫T^{xx} = E, others ≈ 0 | ∫T^{xx}/E within 5% of the lattice expectation E·v_g/v_p (reported both ways); \|∫T^{yy}\|, \|∫T^{zz}\| ≤ 0.05 E |
| massive packet | ∫T^{xx} = E v_g v_p … ≈ E v² | ∫T^{xx}/E within 5% of v_g v_p from the lattice dispersion |

The packet's lattice expectation follows because momentum flux = P v_g and P = E/v_p for a narrow packet. It reduces
to E v² in the continuum (massive) and to E (massless).

## Separate recomputation

The knot's E2 and E4 are recomputed from the regenerated states with an independent closed-form stress (per
orientation, T^{ij} = 2D_ij + 2(trD D_ij − (D²)_ij) − δ_ij(e2 + e4), from T^{ij} = −2∂E/∂γ_ij for a static field; its trace is e4 − e2), and
compared with the numerical γ-derivatives. The packet expectations come from the dispersion relation (closed form)
rather than from the field.

## Compute

Static optimisations: about 2 min (B1 h = 0.2), about 5 min (B2), about 6 min (B1 h = 0.15). Derivatives: seconds.
Total < 0.5 core-h.
