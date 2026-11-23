# GC34J follow-up review and recommended next experiment

2026-11-23 | Assessment after GC34K and GC34L

**Recommendation.** Keep J1 as the linear vacuum-wave benchmark. Next, couple a simple dynamical matter field to J1 and test whether its energy, momentum, and stress source the same geometry that governs its motion.

## What to keep from GC34J

J1 made re-slicing an exact symmetry of the linear lattice equations, retained GC34D’s static field, and gave two near-light-speed waves with no growing mode found in the sampled vacuum spectrum. GC34J’s registered time test failed because it seeded initial conditions that violated the constraints. Its corrected run was post hoc. GC34K then passed a registered time-evolution test with constraint-compatible seeds. J1 is a credible weak-field vacuum benchmark.

## What still limits it

GC34L corrected the moving source’s re-slicing mismatch and used a gauge-invariant curvature readout. The worst curvature error still doubled when grid spacing doubled, reaching 1.55% at kh = 0.05. The bounded search over exact kinetic kernels made little improvement, and the optimized candidates developed short-wavelength instabilities. This supports moving on from stencil tuning for now. It leaves open whether a different lattice action could achieve second-order convergence.

## Recommended next experiment

1. **Start with a scalar wave packet before returning to the knot.** Vary one discrete matter action with respect to J1’s lapse, shift, and spatial variables to obtain energy, momentum, and stress sources. The first gate is whether those sources satisfy the lattice conservation identities as the packet moves. GC34L used a prescribed moving source.
2. **Test whether changing energy changes gravity.** Compare two preparations with the same initial field profile but different kinetic energy. Track the full stress and momentum alongside the energy difference. Check whether curvature agrees with the prediction from independently calculated sources. Use more than one spacing to measure the known first-order error alongside the physical signal.
3. **Let that result guide any geometry redesign.** If the coupled test works, use J1 as a first-order reference model and the knot as a next source. If the intended matter coupling cannot satisfy the lattice conservation identities, redesign the lapse and shift couplings together. GC34K’s centered-gradient obstruction rules out fixing the existing formulation by simply widening that stencil.

The experiment would test whether evolving matter can source the same geometry that governs its motion. The origin of that shared geometry in the original Model Gc rules remains open.

## Records

- [GC34J results](../results/EX-GC34J/EX-GC34J_note.md) and [prior independent review](EX-GC34I-34J_chatgpt_review.md)
- [GC34K registered tests and centered-gradient obstruction](../results/EX-GC34K/EX-GC34K_note.md)
- [GC34L source audit and kernel search](../results/EX-GC34L/EX-GC34L_note.md)
