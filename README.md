# Model Gc — a conserved-G reaction–diffusion model of particles with long-range gravity

**PI:** Richard Benear. **AI collaborators:** Claude, ChatGPT / Codex (roles per `RESEARCH_PROTOCOL.md`).
Branched on 2026-11-18 from the Model G hypergraph study ([RichardBenear/ModelG](https://github.com/RichardBenear/ModelG)).

## Why this programme exists
Paul LaViolette's Model G (Subquantum Kinetics) forms stable particle-like dissipative structures. But in the ModelG
study we found that, as written, a stationary Model G particle has **zero net G**: the ether makes and destroys G
everywhere (−qG + gX), and it rebalances the particle's consumption, so the identity (q − gp)∫φ_G = 0 holds and the G
field is screened (it falls off exponentially). That gives no long-range gravity (EX002–EX007, summaries in
`docs/lineage/`).

**Model Gc** ("c" for *conserved G*) is our own construction, not LaViolette's. It changes one thing structurally:
**G has no source or loss in the ether**. G is supplied only from far away (a fixed G∞ background) and is consumed only
inside particle cores. Conservation then forces an unscreened field G = G∞ − M/(4πr), where M is the particle's
consumption rate. Particles drift toward lower G, so two particles attract with an inverse-square law.

## Model (Gc1, gravity only)
```
X_t = ∇²X + A − (1+B)X + X²Y + c_x·k·G·W
Y_t = 12∇²Y + B·X − X²Y + c_y·k·G·W
W_t = 0.1∇²W − W + Y·W² − W³            (core-only species: W = 0 is the ether state)
G_t = ∇²G − k·G·W,   G → G∞ far away       (no ether source or loss)
```
A = 140/9, B = 29, G∞ = 1. (c_x, c_y) = (1, 0) is "G feeds X"; (0, 1) is "G feeds Y". The X/Y pair is a Brusselator
(the same autocatalytic X–Y loop as Model G), and W marks the particle core. **Gc2** adds a conserved signed charge C
(EX010). The full definition is in `spec/MODEL_Gc.md`.

## Evidence so far
| Study | Track | Result |
|---|---|---|
| EX008 | exploratory | Stable 3D particle; far field G∞ − M/(4πr) to 1e-6; drift toward lower G (l = 1 linear response) |
| EX009 | exploratory | Direct axisymmetric pair solve: 1/d² attraction at 98.5% of prediction (d = 20–30); force ∝ k²; a charge analogue with the right signs but non-reciprocal |
| EX010 | exploratory | Tuned reciprocal charge chemistry (Q₊ = −Q₋, μ₊ = −μ₋); gravity near-universal (G mobilities within 1–2% across particle types) |
| **MG195** | **confirmatory, PI CLOSED** | P1, P2 SUPPORTED (v/prediction 0.98–0.99 for d = 22–40, exponent 1.98); P3 (weak k) UNRESOLVED because of the short-range Turing force; C2 screening control PASS |

These were run in the ModelG repository and copied here with their code and records (see `PROVENANCE.md`).

## Open questions (candidate next studies)
- **Mapping to physics.** The steady G equation is Poisson's equation (G∞ − G ↔ Newtonian potential, M ↔ mass).
  Open issues: the global G budget (a cosmological supply, or a slowly falling G∞), and the propagation speed (the
  field spreads diffusively, not as a wave at finite speed).
- **Discrete version.** A token / hypergraph rule set in 2D or 3D, where G conservation is exact by construction. There
  is no 1D particle for this chemistry, and 1/r needs 3D.
- **Robustness.** Existence window in (A, B, k), stochastic survival (minimum Ω), and P3 retested with the
  short-range part separated.
- **Charge.** Reciprocity without tuning (built-in symmetry), and the gravity-to-charge strength ratio.

## Layout
- `spec/`: model definition and frozen study specifications. `GATES.md`: registered criteria and outcomes.
- `STAGE_LOG.md`, `explore/LOG.md`: stage history and exploratory log. `HANDOFF.md`: current state.
- `src/`: solvers and study scripts. `tests/`: regression tests. `results/`: committed records.
  `paper/figures/`: all PNG figures.
- `docs/lineage/`: notes from the ModelG studies that motivated Gc.
