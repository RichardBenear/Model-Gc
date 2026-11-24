# Handoff to ChatGPT: EX-GC34Q (does a change in kinetic energy change gravity as GR says?)

Claude, 2026-11-24. PI (Richard, 2026-11-23) asked for this handoff so that ChatGPT can run the next test, following
`review/EX-GC34P_chatgpt_review.md`. Ownership of the 34Q study passes to ChatGPT. Claude has run nothing for 34Q.
This is our construction, not LaViolette's Model G.

## Where things stand

- **34N:** the Y1 lattice. Fields are staggered by index parity; the vacuum is linearised Einstein–Hilbert with
  K_a = 2 sin(k_a h/2)/h. Both gauge symmetries are exact.
- **34O:** the scalar's metric couplings. The parallel normal stress is in product form,
  ¼ h_aa(v) G_a(v) G_a(v − e_a). Energy and momentum identities are exact.
- **34P:** the 34O packet sources Y1 gravity.
  - Constraints hold to RK4 error.
  - The tidal curvature matches linearised GR at 0.5–1.2% (h = 0.25), second order.
  - All four are PI CLOSED.
- **Correction to my 34P note.** As the 34P review says, by t = 4 and 8 the two preparations (±velocity) had moved to
  different places. So the A − B difference mixes position with velocity. It is not "the velocity-dependent part" as
  I wrote.

## The test the review recommends (34Q)

Two standing preparations with the **same initial φ** and **different initial φ̇**:
- φ̇(0) = c φ(0), with c = c₁ and c = c₂.
- The mass and gradient energies match; the kinetic energy differs; total momentum is zero.
- The local momentum p = −φ̇∇φ = −c∇(φ²/2) and the stress still differ between the two, and are included.

**Primary observable:** the difference in the gauge-invariant tidal curvature R_0a0b between the two runs, lattice
vs linearised GR driven by the continuum packet with the same initial data. Supporting checks:
- ΔE from the flat scalar Hamiltonian, computed independently, and held fixed across spacings;
- the constraints, checked before interpreting the signal.

**Then (separately, after 34Q passes):** feedback as a controlled weak-field correction, tracking perturbative order
(see the review's last section).

## What to reuse (`src/exgc34p.py`)

- **`Packet`:** the fixed-normalisation packet (amplitude independent of N; keep it).
  - `amp()` builds the mode amplitudes.
  - `field(t, s)` evolves a single travelling branch e^{−isωt}.
  - For 34Q you need general initial data: φ̂(t) = φ̂₀ cos ωt + φ̂̇₀ sin ωt / ω per mode, with φ̂̇₀ = c φ̂₀. Use the
    lattice ω for the lattice and √(m² + k²) for the continuum.
  - Only φ and φ̇ change; the sources then come out unchanged from `exgc34o.msources2` (lattice) and `cont_T`
    (continuum).
- **`lat_f`, `cont_f`:** the source vectors on the mode set, f = 2κ ŝ/V. Each component's phase is referred to its
  own position (`exgc34n.OFFS`).
- **`System`:** the batched Y1 kernel in synchronous gauge.
  - `init` gives the matched minimum-norm constraint initial data.
  - `acc` gives the evolution.
  - `cres` gives the constraint residual.
- **`evolve`:** RK4, with sources evaluated exactly at the stage times.
- **`tidal_err`:** the weighted L2 tidal error (off-diagonal entries counted twice).
- **`recompute`:** a separate per-mode path (`G_up` loop, direct DFT sums, DOP853), to adapt as the recomputation.

## Pitfalls I hit (please avoid)

1. **Constraint rows have a ḣ term.** B_hh = 0 but B_ch ≠ 0. Dropping iB_ch ḣ gave a 51% initial constraint
   residual. The current `cres` keeps it.
2. **Recompute on modes that matter.** My first recompute probed edge modes of negligible size, which is
   meaningless. Pick the largest modes by tidal norm.
3. **Normalisation.** Do not use `exgc34m.packet`: it loses 64× in energy per halving. For 34Q, fix φ₀ by `Packet`'s
   normalisation and check ΔE across h. On the lattice the kinetic energy is exactly ½ c² Σ V φ₀², so ΔE matches the
   continuum to quadrature accuracy without rescaling. Say which you hold fixed.
4. **Periodic box.** Drop q = 0 in both lattice and continuum.
5. **Mode cut.** |n_a| ≤ 25 (|k_a| ≤ 4.9). The dropped source was below rounding for the 34P packet; re-check it for
   the standing packet.

## Cost (cloud container, Linux, Python 3.13, numpy 2.5)

Cost is dominated by the 132,650 evolved modes, not the FFTs. Per 8-time-unit run at dt = 0.05:

| run | cost |
|---|---|
| lattice, h = 0.5 | about 110–150 s |
| lattice, h = 0.25 | about 250–380 s |
| continuum reference (N = 96) | about 240 s |
| h = 1 | about 90 s |

- 34P used 2.1 core-h, over the cap, partly from two runs I stopped.
- 34Q with 2 preparations, 3 spacings and one direction is about 0.4 core-h; both directions about 0.8.
- The dt-halving constraint check doubles the h = 0.5 cost. It may be enough at one direction, since 34P already
  showed it.
- The standing packet's spectrum is centred on ±k₀ with no drift, so a smaller mode cut may be justified. Check the
  dropped fraction first.

## Protocol reminders

- Exploratory track: an `explore/LOG.md` entry and a plan pushed before any run; the 2 core-h cap; the EXPLORATORY
  label.
- Separate recomputation, with shared dependencies stated.
- Never re-score; record deviations.
- Only Richard signs off.
