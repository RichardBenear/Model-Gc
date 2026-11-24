# EX-GC34P — Y1 gravity sourced by the 34O scalar packet (EXPLORATORY)

Claude, 2026-11-24. PI go-ahead (Richard, 2026-11-23): "ChatGPT has review 34O and has a small recommended step. Go
ahead with that test." Source: `review/EX-GC34O_chatgpt_review.md`, "Recommended next experiment".

This plan and the `explore/LOG.md` entry are pushed before any 34P computation. Cap: 2 core-hours on the cloud; I
expect about 0.7. Our construction, not LaViolette's Model G.

## What is tested

The 34O scalar packet's own sources (energy, momentum, stress, from one action) drive Y1's linear gravity
equations. The scalar moves on the flat lattice; there is no gravitational feedback on it (as ChatGPT advises).

1. Fix the inherited packet normalisation so the physical source is the same at every spacing.
2. **Source compatibility:** the forced Y1 equations keep their constraints as the packet moves.
3. **Curvature response:** the gauge-invariant tidal curvature matches linearised GR driven by the continuum packet's
   sources, with matched initial data, and converges as h → 0.

## Set-up

**Packet (fixed normalisation).** Periodic box L = 32, m = 1, |k₀| = 0.6, σ = 3, axis and generic directions (34M).
Mode amplitudes a_k = (2π/L)³ (σ²/2π)^{3/2} exp(−σ²|k − k₀|²/2 − i k·x₀), so the field is about
exp(−|x − x₀|²/2σ²) cos(k₀·(x − x₀)) at every N. The lattice field evolves with the exact staggered dispersion; the
continuum one with ω = √(m² + k²). Each packet is then rescaled to the same total energy E₀ (its own lattice or
continuum energy), as the review asks. The inherited `exgc34m.packet` has an extra 1/N³ (the review's finding); I
report its energy ratio per halving as a check of that finding.

**Two preparations.** Same initial φ; initial φ̇ of opposite sign (A: moving along +k₀, B: along −k₀). B is A
time-reversed, so B's sources at t come from A's field at −t with φ̇ negated.

**Lattice gravity.** Synchronous gauge (n = β = 0) on Y1 at spacing h. For each Fourier mode q ≠ 0 the 6 metric
amplitudes obey C h + iB ḣ − A ḧ = f_h, with the 4 constraint rows C h + iB ḣ − A ḧ = f_c. A, B, C are 34N's kernel
with K_a = 2 sin(k_a h/2)/h. The source is f = 2κ ŝ/V: the DFT of the 34O sources (s_n, s_β, s_aa, s_ab), each
referred to its own position.

**Continuum reference.** The same equations with K → k, and f = κ Pᵀ T^{μν} from the continuum packet's T^{μν},
computed spectrally on an N = 96 grid.

**Initial data (matched).** At t = 0, in both, the minimum-norm (h, ḣ) that satisfies the 4 constraints. That puts
no free (transverse-traceless) radiation in either at t = 0; the gauge parts may differ, but the readout is gauge
invariant.

**Modes.** q = 0 is dropped in both (periodic box). Only modes with every |k_a| ≤ 5 are evolved. The source power
outside that set is reported; I expect it to be below 1e-12 of the total.

**Integrator.** RK4 with step dt (sources evaluated exactly at the stage times), t from 0 to 8, dt = 0.05.

**Readout.** In synchronous gauge the tidal tensor is R_0a0b = −½ ḧ_ab, which is gauge invariant at linear order.
It is compared per mode (amplitudes referred to each component's position) as a relative L2 error over all evolved
modes and components, at t = 4 and 8, for A, for B, and for the difference A − B.

## Gates

- **P0, normalisation.**
  - Before rescaling, the fixed packet's lattice energy at h = 0.5 and 0.25 agrees to ≤ 1%.
  - After rescaling, the width (rms radius of the energy density) and central wavenumber (energy-weighted mean |k|)
    at h = 0.5 and 0.25 agree with the continuum to ≤ 1%.
- **P1, source compatibility** (h = 0.5, both preparations, both directions).
  - 34O sources: constraint residual (relative to the constraint source) ≤ 1e-6 at dt = 0.05, and falls by ≥ 8×
    at dt = 0.025. That shows any residual is integrator error, not source error.
  - Null control, 34N sources: residual ≥ 1e-4, and changes by ≤ 20% when dt halves.
  - At h = 0.25, dt = 0.05: 34O residual ≤ 1e-6.
- **P2, curvature response.** Tidal error against continuum GR, for A, B and A − B, both directions, t = 4 and 8:
  - ≤ 2e-2 at h = 0.25;
  - ratio error(h = 0.5)/error(h = 0.25) ≥ 3.
  - h = 1 is reported only (its grid cannot hold the packet's spectrum).

**Overall:** PASS if P0, P1 and P2 hold.

**Reported, not gated:**
- ‖R_A − R_B‖/‖R_A‖, to show the velocity-dependent part being tested is not negligible;
- the constraint residual at h = 1;
- the dropped source fraction.

## My expectations, stated before running

- P0 and P1 pass.
- P2 error around 5e-3 at h = 0.25, falling about 4× per halving.

Risks:
- h = 0.5 may still be pre-asymptotic for the curvature (as in 34N's E-M2);
- the matched initial data may not be equivalent if a constraint row loses rank at some mode.

## Not in this study

- Gravitational feedback on the scalar.
- The knot.
- Where Y1's geometry comes from.
- Any change to the white paper.

## Files

`src/exgc34p.py`:
- `norm` → `results/EX-GC34P/norm.json` (P0);
- `compat` → `compat.json` (P1);
- `curv` → `curv.json` (P2);
- `recompute` → a separate check.
