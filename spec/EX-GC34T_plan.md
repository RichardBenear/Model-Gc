# EX-GC34T — continuum next-order identity and the periodic mean (EXPLORATORY)

2026-11-24 | ChatGPT | PI authorization in the request to verify both items. Cap 2 core-hours.

## Scope

Verify the conservation identity a second-order gravity source must satisfy, and test a continuum prescription for the homogeneous mode. This is a prerequisite audit. Y1 currently supplies a linear lattice operator, so this does not certify a second-order Y1 action or a reciprocal packet evolution.

## T1. Next-order identity

In mixed-index components and with the coupling parameter factored out, let the first-order Einstein equation be G1 = κ T0. Expand the contracted Bianchi identity and covariant matter conservation. Verify exactly with symbolic differentiation that the divergence of G2 equals the connection correction built from G1. Use a periodic synchronous manufactured metric perturbation with h_xx = (1+t) cos x and h_yy = h_zz = (t+t²) sin x. Build the Einstein coefficients independently of the conservation identity, by expanding the connection and Ricci tensor. Set the manufactured T0 to G1/κ solely for testing the identity. Pass if the full divergence identity cancels identically in every component, and the connection correction is nonzero for at least one component.

As a negative control, subtract a constant mixed tensor C^0_0 = -1 from T0 only and keep the same metric. Determine whether the resulting conservation mismatch survives removal of its spatial mean. Pass if the control shows nonzero spatially varying mismatch. This is a consistency test of source-only subtraction, not a physical model of a compensating sector.

## T2. Homogeneous mode

Derive the spatially averaged linear Hamiltonian constraint about a static flat periodic background and test it against positive density. Use a spatially flat expanding homogeneous metric on a three-torus and a homogeneous massive scalar as a separate continuum benchmark. Set κ=1, m=1, a(0)=1, φ(0)=φdot(0)=1. For coupling ε=1e-4 and 2.5e-5 initialize H using the expanding Friedmann branch. Evolve a, H, φ and φdot to t=2 using the acceleration and scalar equations, and check the Friedmann constraint independently over 101 equally spaced times with relative residual <=1e-8. Check H(0) ratio 2 within 1e-10. Confirm choosing H(0)=0 would leave a unit relative initial constraint residual for either positive coupling.

## Verdict

Pass if T1 and T2 pass. Distinguish a verified continuum identity and admissible homogeneous background from an actual lattice construction. Report how either an expanding background or an explicitly conserved compensation prescription changes the earlier periodic flat-box setup.
