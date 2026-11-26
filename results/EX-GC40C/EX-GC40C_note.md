# EX-GC40C — is a shared-entanglement hub a mass, or only a local defect?

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC40C_plan.md), with predictions written after 40A and following the [Codex review](../../review/EX-GC40AB_codex_review.md). PI closure pending. Cloud container, Python 3.13, numpy 2.5. About 0.2 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

**A hub is a local defect, not a gravitational source.**
- **It has some properties of a mass.** It lowers the ground-state energy, it traps bound states in the gap, and two hubs attract.
- **But every influence dies off exponentially.** The decay is faster than the background's own correlations (about half their decay length), with no 1/r tail, at box sizes 16, 24 and 32.
- So a network of entangled pairs, however much entanglement its sites share, cannot make long-range gravity by itself. A long-range field needs a gapless carrier coupled to energy.

## Results

Hub: on-site coupling λ₀ + Δλ on one site (or a 7-site cluster). Everything is computed by Woodbury updates of the background Green's function on the imaginary axis.

**Energy, bound states and the two-hub interaction.** These are the same to all printed digits at L = 16, 24 and 32.

| λ₀ (gap) | ΔE₀ at Δλ = 0.1 / 0.2 / 0.4 (one site) | ΔE₀ (7 sites, Δλ = 0.2) | In-gap states at Δλ = 0.2 (μ = midgap) | E_int(d) for two hubs, d = 2 / 3 / 4 / 5 |
|---|---|---:|---|---|
| 0 (2.0) | −0.0370 / −0.1430 / −0.5132 | −1.006 | −0.618, 0.908 (×5) | below 1e-13 at every d |
| 0.1 (1.4) | −0.107 / −0.271 / −0.722 | −1.905 | −0.443, 0.822 (×5) | −6.0e-6 / −3.5e-8 / −2.8e-10 / −2.5e-12 |
| 0.2 (0.8) | −0.166 / −0.377 / −0.894 | −2.646 | 0.723 (×2), 0.725 (×3) | −9.3e-6 / −5.2e-7 / −4.1e-8 / −3.0e-9 |

**Far field** (change of link energy along the x axis, one-site hub, Δλ = 0.2, L = 32):

| r | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|---|
| λ₀ = 0 | 0 (< 1e-33) everywhere beyond the hub's own links | | | | | | | | |
| λ₀ = 0.1 | 9.7e-5 | 4.8e-7 | 3.4e-9 | 3.1e-11 | 3.4e-13 | 4.3e-15 | 8.8e-19 | 2.2e-22 | 6.5e-26 |
| λ₀ = 0.2 | 3.1e-4 | 5.8e-6 | 2.7e-7 | 2.2e-8 | 2.1e-9 | 2.3e-10 | 3.2e-12 | 5.3e-14 | 9.9e-16 |

**Fits over r = 3 to L/2 − 2.** The model is A·e^{−r/ξ}/r^p. The background correlation length ξ_C is fitted the same way to |C(r)|.

| λ₀ | ξ_hub (1 site / 7 sites), L = 32 | ξ_C | 1/r residual ÷ exponential residual |
|---|---|---:|---:|
| 0.1 | 0.264 / 0.270 | 0.523 | 1.3e6 – 1.7e6 |
| 0.2 | 0.535 / 0.581 | 1.345 | 3.0e3 – 5.5e3 |

## Scoring

- **P1 (locality control at λ₀ = 0): PASS.**
  - Beyond the hub's own modes and their link partners, every change is below 1e-33 (Woodbury) or 8e-16 (dense, L = 8).
  - The two-hub interaction is −9e-14 at every d, a constant floor, so it is zero.
- **P2 (ΔE₀ < 0 and quadratic at λ₀ = 0): PASS.** ΔE₀(0.2)/ΔE₀(0.1) = 3.862, within 3.5% of 4.
  - At λ₀ > 0 the first-order term is negative too (recorded, not predicted).
- **P3 (a bound state in the gap at Δλ ≥ 0.2, λ₀ = 0.1): PASS.**
  - The hub pulls one state up out of the valence band (to −0.44, below μ = 0.2) and five down from the conduction band (to 0.82).
  - At λ₀ = 0.2 only the five conduction-side states appear inside the gap.
- **P4 (exponential decay in a gapped background, no 1/r): main clause PASS, ξ clause FAIL at λ₀ = 0.2.**
  - Every far-field profile and E_int(d) falls exponentially. The decay is the same at L = 16, 24 and 32, apart from periodic images at L = 16.
  - A 1/r form is rejected by 3e3 to 2e6 in residual, against the 10× gate.
  - The registered "ξ within a factor 2 of ξ_C" holds at λ₀ = 0.1 (ratios 1.98 and 1.94) but not at λ₀ = 0.2 (2.51 and 2.31).
  - The reason (inference, after the data): a change in the state is quadratic in the background propagator (G₀ … G₀), so it decays with ξ_C/2. The prediction should have said ξ ≈ ξ_C/2. The measured ratios, 1.9–2.5, are close to 2.
- **P5 (gapless λ₀ = 0.5, descriptive): unresolved.**
  - At half filling the dense spectrum has degenerate states at the Fermi level (splitting 2e-15) at both L = 8 and L = 10. So the ground state is not unique, and the profiles at the two sizes disagree in sign and size.
  - No exponent can be read off.
- **P6 (a local excitation leaves the far field alone): PASS, structurally.**
  - The link-energy change follows the packet's own weight and falls faster than it: 1.4e-19 at r = 1 to 2.8e-37 at r = 10, against packet weights 0.98 to 2e-19.
  - Free fermions carry no influence beyond the packet.
- **Readout 5 (do far-field amplitudes scale with ΔE₀ or Δλ): not run.** The far field was computed at Δλ = 0.2 only. This is a deviation from the plan.

## Checks

- **C1: PASS.** At L = 8, Woodbury matches dense diagonalisation:
  - ΔE₀ to 8e-12 (one site) and 6e-11 (7 sites);
  - the site blocks of δP to 6e-16.
- **C2:** structural. Δλ = 0 gives D = 0, so the response is exactly zero.
- **C3: FAIL, marginally.** 160 and 320 quadrature points agree to 2.4e-10 for ΔE₀ (gate 1e-10) and to 2e-15 for the far field.
  - This floor sets the E_int noise at large d (about 1e-12 to 1e-13). It does not affect anything above 1e-11.

## What it means (inference)

- **A hub has the local signatures of a particle but none of a gravitational source.**
  - It binds states, costs (or here releases) energy, and two hubs attract at short range.
  - But in a gapped pair network every effect is confined to about half the background correlation length. This holds at every level of sharing tested.
  - Free (Gaussian) networks also cannot let an excitation's energy reach anything beyond its own wave packet.
- **This closes the "pure pair network makes gravity" line.** The range of any influence is set by the gap, and the gap is what makes the entanglement area-law and local. A 1/r field needs a gapless mode that couples to energy.
  - In our programme, that is the lapse field: relativistic (37C), not free (39C), and acting on light by event counting (37B, 40B, 41A).
  - In Chain Theory's terms, it would be a count that every read passes through, which the compendium does not yet derive (its field equation is still open).
- **The next design question.** Couple a 40A-type pair network to a lapse field whose source is the network's energy density. Then ask whether a hub's ΔE₀ sources a 1/r lapse with a universal coefficient. That would be a source law, not a defect.

## Files

[Code](../../src/exgc40c.py) · `checks.json`, `main.json`, `gapless.json`, `control.json`
