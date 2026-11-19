# Handoff — current state (2026-11-18)

**Owner of the active work:** Claude (set-up). **No study is running.**

## State
- The repository was set up on 2026-11-18. It holds the Gc1/Gc2 code, records and notes copied from ModelG (source
  commit edab7d47; see `PROVENANCE.md`), the model definition (`spec/MODEL_Gc.md`), and the protocol.
- Inherited results:
  - EX008–EX010 (exploratory).
  - MG195 (confirmatory, PI CLOSED 2026-11-18): P1 and P2 SUPPORTED, P3 UNRESOLVED, C2 PASS.

## Candidate next steps (none approved yet)
1. **Physics mapping note:** done, `docs/PHYSICS_MAPPING.md`. Suggested tests from it: weighted action-reaction with an unequal pair; drift bifurcation (inertia); telegraph G carrier (finite speed); G budget.
2. **Stochastic survival (EX-GC01):** noisy 3D continuum (tau-leap or Langevin), to find the minimum Ω at which a Gc
   particle survives. This is a prerequisite for a token version.
3. **Token / hypergraph Gc (2D, then 3D):** G tokens only hop and are destroyed only at W, so conservation is exact.
4. **Existence window** in (A, B, k, c), and the P3 retest with the short-range part separated (registered).
5. **Gc2 with built-in charge symmetry**, and the gravity/charge strength ratio.

## Update 2026-11-18: EX-GC01 done
Weighted action and reaction holds for particles with equal mobility μ; free fall is not universal when μ differs
(doubled k: 1.8 times faster). New design target: a universal μ across particle types (results/EX-GC01/EX-GC01_note.md).

## Update 2026-11-18: EX-GC02 done
Inertia appears when the core species W responds slowly (drag → 0 at τ_W ≈ 28; effective mass confirmed by a direct
time response), but the threshold is type-dependent, so free fall becomes non-universal. Revised design target: one
quantity setting both G consumption and inertia (results/EX-GC02/EX-GC02_note.md).

## Update 2026-11-18: EX-GC03 done
Push ∝ M (active = passive gravitational mass) but drag is set by the shared X/Y core, so free fall differs by
type (35%). Suggested next: a bound pair of identical particles, to see whether binding breaks universality
(results/EX-GC03/EX-GC03_note.md).

## Update 2026-11-18: EX-GC04 done
Identical particles bind at d = 7.145 and 10.11 (short-range Turing force), and a bound pair falls 6% (3%) more slowly
than a single particle. Cause: consumption grows with local G, so the upper partner pulls harder; the deficit fades as
(dM/dG)/(8πd). In Gc1 the same G-dependence that makes particles fall breaks the equivalence principle for composites
(results/EX-GC04/EX-GC04_note.md). Next idea: a variant where drift does not come from G-dependent consumption.

## Update 2026-11-18: EX-GC05 done (Stage A)
No rotating Gc1 particle: non-round core modes never oscillate; the X/Y chemistry oscillates only as a whole medium
(ether Hopf τ_Y = 8.64). Route 1 to a 720° particle needs a core with its own clock (e.g. a slow inhibitor of W);
routes 2 (orientation field) and 3 (hypergraph ribbons) need new ingredients too (results/EX-GC05/EX-GC05_note.md).

## Update 2026-11-18: EX-GC06 done
A slow local core inhibitor V ("Gc3") gives the particle a breathing clock (period ≈ 18–20) near the existence
edge (κ ≈ 0.9–1). Rotation does not come first: l = 1 oscillates only later, and l = 2, 3 never. Next would be a 3D
time simulation (does breathing saturate, does l = 1 become circling, any 2:1 lock) (results/EX-GC06/EX-GC06_note.md).

## Update 2026-11-18: EX-GC07 done
3D (spectral-in-angle) time runs of Gc3: breathing grows until the core collapses; no steady pulse seen, no circling.
Route 1 would need a saturating brake (another tuned ingredient); routes 2/3 get 720° from structure instead
(results/EX-GC07/EX-GC07_note.md). Solver src/exgc07_sim3d.py (8 ms/step) is reusable for other Gc variants.

## Update 2026-11-19: route 3 design note and EX-GC08 done
docs/ROUTE3_RIBBON_DESIGN.md: 720° from a ribbon-structured Gc hypergraph tethered by its G supply. EX-GC08 (Stage 1):
a local three-node rule (no edge may swing through the half-turn ambiguity) on cube-rotation frames reproduces the
belt trick (1 turn stuck, 2 turns untangle); pairwise rules fail. Next: Stage 2, Gc tokens on this ribbon graph
(results/EX-GC08/EX-GC08_note.md).

## Update 2026-11-19: EX-GC09 done (route 3, Stage 2)
Gc1 as tokens on a 3D cubic lattice (h = 0.5): the particle survives at Ω = 1000 (M −6%) and 300 (−20%) over T = 30,
fades at 100, dies at 30 (about 10⁴ core tokens needed). Far G deficit matches the lattice point-sink field (1–4%).
A splitting bias in the first token update was found and fixed (results/EX-GC09/EX-GC09_note.md). Next: Stage 3,
tethered particle (needs the PI's go-ahead).

## Update 2026-11-19: EX-GC10 done (route 3, Stage 3)
Particle tethered by its G supply (frames on the EX-GC09 lattice, no-cut lift rule, twisted edges pass G with
probability 1 − ε). One turn: twist stuck next to the core, index −1, M −7.9% (deterministic, 3 seeds); two turns
untangle, M within 0.6%. Controls (ε = 0, pairwise rule) as predicted. Token runs (Ω = 3000, partly on the PC) agree on
topology but the mass signal is noisy (shared drift artefact). Open: exchange sign, a real rotating particle.

Overall status summary as of 2026-11-19: docs/STATUS_2026-11-19.md.
