# ChatGPT review of EX-GC34P

2026-11-24 | Independent review of the exploratory results

**Verdict.** GC34P passes. Y1 can take energy, momentum, and stress from an evolving scalar wave and produce the expected weak gravitational field. The next test should keep the same initial scalar field profile, change its kinetic energy, and measure the change in gravity. GC34P reversed the packet's direction at equal energy and answered a different question.

## What GC34P established

- **Packet normalization.** The old packet lost about a factor of 64 in energy each time the grid spacing halved. Before explicitly matching the energies, the new packets differed by only 0.21% between the scored spacings. The curvature comparison across grids is meaningful.
- **Source compatibility.** With GC34O's corrected stress, halving the timestep cut the constraint residual by 16×, as expected for RK4. With GC34N's original stress, the constraints missed by 4–8% and timestep refinement barely changed the result. This control distinguishes source inconsistency from time-stepping error.
- **Tidal curvature.** At \(h=0.25\), the committed relative errors were 0.48–1.18%, with about a factor of four reduction on halving the spacing from 0.5. This held for both directions of travel and their difference. At the reported times, the difference includes the packets having moved to different places. It tests the time-dependent field from the full source rather than isolating one gravitomagnetic component.

## Recommended next experiment

**Use two standing scalar packets with the same initial \(\phi\) and different initial \(\dot\phi\).** For example, set the initial field velocity proportional to the field profile in both cases, with different constants. The initial mass and gradient energy terms then match while kinetic energy changes and total momentum is zero. Include the local momentum and stress, which can still differ, in the continuum comparison.

Register the difference in gauge-invariant tidal curvature as the primary observable. Calculate the energy difference independently from the flat scalar Hamiltonian, preserve that difference across spacings, and check the gravity constraints before interpreting the signal. Reuse GC34P's matched periodic initial conditions and evolve only the modes needed for the comparison.

**Then test feedback as a controlled weak-field correction.** GC34P lets the scalar source gravity while the scalar evolves on the flat grid. After the energy-change test passes, use GC34O's matter action to calculate the leading change in scalar motion caused by its own field. Keep track of perturbative order before using the corrected scalar as a gravity source. GC34O's exact source balances were derived on the flat background. Repeatedly inserting a gravity-modified source into the unchanged linear equations could violate the constraints GC34P validated.

## Scope and records

I reviewed the plan, implementation, and committed data. I checked the recorded curvature errors and did not rerun GC34P.

- [GC34P plan](../spec/EX-GC34P_plan.md) and [results](../results/EX-GC34P/EX-GC34P_note.md)
- [Normalization](../results/EX-GC34P/norm.json) and [constraint comparison](../results/EX-GC34P/compat.json)
- [Axis curvature](../results/EX-GC34P/curv_axis.json) and [off-axis curvature](../results/EX-GC34P/curv_generic.json)
- [GC34O source derivation](../spec/EX-GC34O_plan.md) and [GC34P solver](../src/exgc34p.py)
