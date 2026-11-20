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
inside particle cores. Conservation then forces an unscreened field G = G∞ − S/(4π D_G r), where S is the particle's
consumption rate (written M in the earlier notes). Particles drift toward lower G, so two particles attract with an inverse-square law.

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

## Where the programme stands (2026-11-20)
The full account is the white paper **[`paper/Model_Gc_white_paper_v0.3.md`](paper/Model_Gc_white_paper_v0.3.md)**,
frozen as covering the research through EX-GC25 after two rounds of review by ChatGPT. Everything except MG195 is
EXPLORATORY and supports no claim on its own.

| Line of work | Studies | Where it stands |
|---|---|---|
| Gravitational sourcing | EX008–EX010, **MG195** (confirmatory, PI CLOSED), EX-GC09 | 1/r field and an inverse-square **drift** law (v ∝ 1/d²; P1, P2 SUPPORTED, P3 UNRESOLVED, C2 PASS). Survives as tokens on a 3D lattice. Not shown: acceleration (a ∝ 1/d²). |
| Equivalence principle | EX-GC01–04, EX-GC16 | Non-universal in Gc1 (type and binding). A recoil rule gives a universal response per unit of *imposed* inertia in the weak regime. |
| Spin-like topology | EX-GC05–15 | Rotation route failed. A degree-1 lattice knot is stable, keeps its degree, has the 360°/720° loop classes and the classical (Finkelstein–Rubinstein) exchange signature, and can carry a Gc particle that sources 1/r and drifts correctly. Not quantum spin. |
| Inertia | EX-GC17–24A | Absorbing G from a medium at rest gives drag, not inertia. With imposed inertia and a transport term the particle coasts (EX-GC23); the medium's own inertia gives coasting times of 0.043 / 0.0022 under the assumed drag law (EX-GC24A, conditional negative). |
| Frame-free direction | EX-GC25 | Net-absorption gravity defines a rest frame, which causes the drag. In 1+1 D, energy-coupled kinks fall alike in an external potential and coast. A self-gravitating frame-free particle is untested. |

**Next (PI decision):** the seven-step gated sequence in Section 10 of the white paper: a conservative 3D topological
particle, its rest energy and independently measured inertia, its boosted energy–momentum relation, then a reciprocal
gravitational field, coasting with that field, gravitational vs inertial mass, and two-particle dynamics. The
alternative is to keep Model Gc's consumption picture and pursue internal inertia.

Studies EX008–EX010 and MG195 were run in the ModelG repository and copied here with their code and records (see
`PROVENANCE.md`). EX-GC22–24A were designed and run by ChatGPT and reviewed by Claude; the rest by Claude.

## Layout
- `paper/`: the white paper (v0.3, frozen; earlier drafts v0.1 and v0.2 kept). `paper/figures/`: all PNG figures.
- `spec/`: model definition, study plans and frozen specifications. `GATES.md`: registered criteria and outcomes.
- `STAGE_LOG.md`, `explore/LOG.md`: stage history and exploratory log. `HANDOFF.md`: current state.
- `src/`: solvers and study scripts. `tests/`: regression tests. `results/`: committed records and a
  note per study. `review/`: reviews of studies run by the other AI.
- `docs/`: physics mapping, design notes (route 3, soliton particle), the 2026-11-19 status summary (superseded by
  the white paper), and `docs/lineage/`, notes from the ModelG studies that motivated Gc.
