# Projected Y1 reciprocal scalar benchmark

Use this recipe for the compensated periodic weak-field scalar benchmark. It carries the first scalar correction back into the second-order metric and controls the synchronous constraints at the accepted times.

## Set up the reference calculation

1. Use L=12 and a standing massive packet with m=1, nominal carrier k=.6 in the generic direction, width 2 and initial velocity equal to the initial field. Keep scalar input Fourier indices through 2 per axis and both metric orders through 8 per axis, excluding the homogeneous metric. Use the original GC34O source, the 34R tangent scalar evolution, and the symmetric densitized second-order source of 35A/C. Run the matched periodic continuum spectral control with continuum derivatives and the same Fourier support.
2. Keep the compensated zero-mode prescription of 35A/C. A spatially uniform auxiliary stress cancels the leading mean source and the auxiliary energy and current are evolved at the next order. The benchmark therefore describes nonzero modes on a particular periodic background.
3. Use physical spacing h=.125 for the refined check, with h=.25 and .5 for spatial controls. Evaluate band-limited products on the 18³ quadrature, using physical translations exp(ikₐh) and the Y1 derivative Kₐ=2sin(kₐh/2)/h. Keep these two grids conceptually separate. Use RK4 with Δt=.05 through t=20 and a smaller timestep control at matched physical times. The source-strength scale is ε=10⁻³.

## Project and record the corrections

For each nonzero mode and each accepted step, form the instantaneous second-order sourced residual

```
r = RC h₂ + RB v₂ + RA f₂,spatial − f₂,constraint.
```

Let W be diagonal with weights (1,1,1,2,2,2) in the six-component tensor basis, with the same weights for velocities. In coordinates z=(W¹ᐟ²h₂,W¹ᐟ²v₂/|K|), form C=[RC W⁻¹ᐟ², |K| RB W⁻¹ᐟ²]. Compute Δz=−C⁺r using the per-mode pseudoinverse and transform back to metric and velocity. Apply the correction to the second-order state only. This tensor norm preserves the transverse-traceless parts to roundoff. An equal weighting of the six packed entries does not have that property.

Record the normalized residual immediately before and after projection, the metric and velocity correction norms at each step, and their accumulated sums. Normalize the residual by the original constraint-source norm at that time and multiply by ε. Keep the free run with the same RK stage sources as a control. A smaller timestep shortens the interval in which a defect can build up between resets. The pre-reset value therefore depends on timestep even for an unchanged finite-spacing source.

## Read out the curvature

At t=4, 8, 12, 16 and 20, compare all six components and the transverse-traceless (TT) part against the matched continuum tensor. Reconstruct −∂²ₜh₂/2 from the accepted positions to include the projection resets. Use centered second differences at interior samples and a second-order backward difference at the endpoint. Record the nominal RHS acceleration separately. For the synchronous coordinate curvature R₀ᵢ₀ⱼ at second order, add (∂ₜh₁)²/4 as a matrix product. Also report the velocity reset divided by twice the timestep, in units of the continuum tidal norm.

Report the results as a finite-time, compensated, perturbative benchmark. A roundoff residual after projection follows from enforcing the constraint surface. Source conservation, curvature accuracy, and the precision of the unprojected equations are separate checks. The measured window and source do not establish an isolated mass, arbitrary-packet stability, or exact nonlinear lattice symmetry.

At h=.125 the pre-reset peak was 4.83e−4 and the post-reset peak was at roundoff. The reconstructed tidal error was 2.29% at t=16 and 2.25% at t=20. It missed the original 2% accuracy gate. The result supports the projected benchmark at a stated few-percent scale for this profile and window.

## Reproduce

Run from the repository root. The drivers generate coefficient snapshots and per-step JSON locally.

```text
python src/exgc35e.py coarse
python src/exgc35e.py fine --n 48
python src/exgc35e.py continuum --continuum
python src/exgc35e.py half --dt .025 --end 4
python src/exgc35e_curvature.py
python src/exgc35e_analysis.py
python src/exgc35e.py finer --n 96 --output EX-GC35G
python src/exgc35g.py
```

[GC35C free-evolution comparison](../results/EX-GC35C/EX-GC35C_note.md) · [GC35E projection setup](../results/EX-GC35E/EX-GC35E_note.md) · [GC35G refinement](../results/EX-GC35G/EX-GC35G_note.md)
