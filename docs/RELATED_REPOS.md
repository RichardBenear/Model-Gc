# Related repositories and reusable tools (surveyed 2026-11-18)

Nothing below has been copied into Model Gc yet. Cite the source commit and file hash when adapting anything, and
revalidate it here: historical results do not transfer to modified code.

## RichardBenear/ModelG (parent programme)
- Already copied: the Gc solvers and records (see `PROVENANCE.md`).
- Candidates for the token / hypergraph Gc:
  - the tick-based token hypergraph engine (MG104+), the parent's causally invariant engine;
  - the Phase 2 tau-leap kernel `src/lattice.py` (`tau_leap_euler_1d`, the LNA);
  - tracking (`src/measure.py`) and deterministic confirmation (`src/confirm.py`);
  - seeding (`runner.seed_for`, entropy 20261027).
  - All of these are 1D and would need a 2D/3D version.

## RichardBenear/StochasticVacuum (commit 2161211, archival)
Classical stochastic graphs with finite, conserved token stores. The repository's own evidence boundary applies: it
established no native particles, and its controls validate apparatus only. Most relevant here:
1. `materials/pc_rd1b_l0b_narrowed_local_kit/.../pc_rd1b_l0b_v35_seedguard.cpp` handles explicit reservoir occupancy,
   channel transfers, depletion, interventions and lineage tracking. A reservoir exchanging tokens with the lattice is
   the discrete form of Gc's far-field G∞ supply.
2. The same kit's `analysis_core_v31.py` and `narrow_analysis.py`: lifetime, event-accounting and control analyses
   (structure survival under noise).
3. `materials/dn01d_v4_local_stage_kit/`: a frozen tracker (finite-horizon and proximity limits recorded).
4. E2/E3 kits (in `archives/`): separate spatial and exchange graphs, conserved finite tokens, graph-acceptance
   feedback. These are relevant if Gc moves from a fixed lattice to a rewriting hypergraph.

## RichardBenear/CRR-Investigation
On 2026-11-18 the repository was reachable but had no commits (empty), so there was nothing to survey. The CRR handoff
document is preserved in StochasticVacuum (`materials/STOCHASTIC_VACUUM_PROGRAM_HANDOFF.md`).
