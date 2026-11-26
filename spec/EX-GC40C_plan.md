# EX-GC40C — is a shared-entanglement "hub" a mass, or only a local defect?

2026-11-26 | EXPLORATORY | Claude | Predictions written after 40A's results, as the [EX-GC40 plan](EX-GC40_plan.md) requires. They follow the [Codex review](../review/EX-GC40AB_codex_review.md) (its test 3). PI request: "yes, go ahead" (2026-11-26 16:38 UTC), for writing the plan. **Running needs a separate go-ahead.** Proposed cap 4 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

## Question

The Chain Theory compendium builds mass from reads through "hubs", systems shared by more than one link (its Part X). In 40A, matter registers in the entanglement only when pairs share entanglement with their neighbours (λ > 0).

Put the sharing in one place, a hub. Two questions:
- Does the hub carry energy and hold a localised (bound) excitation, as a mass would?
- Does it change the network around it with a long-range (1/r) profile, as a gravitational source would? Or is the change short-ranged, like a local defect?

## Model

- 40A's network: six fermion modes per site, one per link end, with a pair term on each link. Half filling.
- **Background:** uniform on-site coupling λ₀.
  - λ₀ = 0: pure pairs, the locality control.
  - λ₀ = 0.1 and 0.2: gapped, with gap 2 − 6λ₀ = 1.4 and 0.8.
  - λ₀ = 0.5: gapless, descriptive.
- **Hub:** λ₀ + Δλ on one site (main), and on a 7-site cluster (the site and its six neighbours) as a repeat. Δλ ∈ {0.1, 0.2, 0.4}.
- **Local-excitation control:** in the unchanged background, one particle–hole excitation of the background Hamiltonian, the lowest-energy one localised as a wave packet near the centre.

**Numerics.** The hub changes h on at most 6 × 7 = 42 modes, so it is a rank ≤ 42 update.
- Everything is computed from the translation-invariant Green's function by Woodbury updates on the imaginary axis, as in 40A's G4. This reaches box sizes L = 12, 16, 24 and 32.
- Bound states: zeros of det(I − G₀_RR(E) δh) for real E inside the gap.
- At λ₀ = 0.5 there is no gap, so a dense diagonalisation is used (L ≤ 12).

## Readouts

1. **Ground-state energy shift** ΔE₀(Δλ), from the Woodbury trace formula. Its sign and its scaling with Δλ.
2. **Bound states:** the number and energies of states the hub pulls into the gap.
3. **Far field.** At each distance r from the hub (sites along the three axes and the diagonals), the change in:
   - (a) the bond energy of each link;
   - (b) the single-site entropy (all six modes);
   - (c) the entropy carried by each link, measured as the entropy of the two end modes.

   Fits to A·e^{−r/ξ}/r^p and to B/r, over r = 3 to L/2 − 2. The decay is also compared across L = 16, 24 and 32 to separate it from finite-box effects.
4. **Second probe.** The interaction energy of two hubs a distance d apart, E_int(d) = E(both) − E(one) − E(other) + E(none), for d = 2–12, compared against exponential and 1/d.
5. **Linearity in the source.** Do far-field amplitudes scale with ΔE₀ (energy as the source) or with Δλ?

## Predictions

- **P1 (locality control).** In the pure-pair background (λ₀ = 0), every change is confined to the hub's own modes and the modes its links reach. Beyond distance 1, the far field and the interaction energy are exactly zero, below 1e-12.
- **P2 (energy).** At λ₀ = 0, ΔE₀ < 0 and quadratic in Δλ (within 5% over Δλ = 0.1–0.2).
  - Reason: in the pure-pair state, modes of different links on the same site are uncorrelated, so H₁ has zero first-order expectation. The second-order correction to a ground-state energy is never positive.
  - At λ₀ > 0 there is a first-order term. Its sign is recorded, not predicted.
- **P3 (bound state).** Expectation, recorded: at least one hub state appears in the gap for Δλ ≥ 0.2 at λ₀ = 0.1.
  - Reason: a one-site matrix perturbation of size comparable to the gap.
  - A hub that binds a state is the closest analogue of a "mass" this network has.
- **P4 (no long range in a gapped background; the main prediction).** For λ₀ = 0.1 and 0.2, every far-field profile and E_int(d) decays exponentially. The fitted ξ is within a factor 2 of the correlation length of the background correlator C(r), measured from 40A's FFT correlators. A pure 1/r fit is rejected (its residual exceeds the exponential fit's by more than 10×).
  - Expectation: **passes.** In a gapped free-fermion system, every local response is short-ranged. A hub is then a defect, not a source.
- **P5 (gapless background, descriptive).** At λ₀ = 0.5 the response may be a power law. Recorded: the exponent and any oscillation (Friedel-like). A 3D Fermi surface would give cos(2k_F r)/r³ for densities, not 1/r.
- **P6 (control).** A local excitation of the unchanged background leaves every far-field observable unchanged beyond the extent of its own wave packet. Free fermions have no interaction to carry its influence.

## What a result would mean

**If P4 passes, as expected.**
- A hub in this network is a local defect.
- Free (Gaussian) pair networks cannot make gravity, at any level of entanglement sharing. A long-range field needs a gapless mediator coupled to energy.
- In our programme that mediator is the lapse wave field of 37C and 39C (constrained, not free). In Chain Theory's language, a count that every read passes through.
- The next step would then be to couple a 40A-type network to such a mediator, not to search inside free pair networks.

**If P4 fails**, with a 1/r tail stable across L. A free pair network with hubs would then carry a long-range response, and a source-law study (does the amplitude scale with ΔE₀?) would be justified.

## Checks

- **C1.** Woodbury results equal dense diagonalisation at L = 8, for one site and seven sites, to 1e-10.
- **C2.** Δλ → 0 gives zero response. Δλ = 0 on a 7-site cluster at λ₀ = 0 reproduces 40A's exact one-ebit-per-link state.
- **C3.** Imaginary-axis quadrature is converged: 160 vs 320 points agree to 1e-10 for ΔE₀ and for the far field at r = L/2 − 2.

## Outputs

`src/exgc40c.py`, `results/EX-GC40C/`.
