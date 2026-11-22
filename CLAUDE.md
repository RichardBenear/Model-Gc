# Model Gc — standing instructions for Claude

Read automatically at the start of every Claude Code session here. Current state and next steps are in `HANDOFF.md`.

## What the project is
A research programme on **Model Gc**, our conserved-G reaction–diffusion model (`spec/MODEL_Gc.md`). It branched on
2026-11-18 from the Model G hypergraph study (RichardBenear/ModelG). The aim is to test whether particles with a
conserved, core-consumed carrier give long-range (1/r, inverse-square) gravity and a charge analogue, and eventually a
discrete, causally invariant rule version (the parent programme's goal of a link to quantum mechanics and gravity).
Model Gc is our own construction: always say so, and never present it as LaViolette's.

## Roles and protocol
- **PI: Richard Benear.** Direction, registration, budgets, closure sign-off.
- **AI collaborators:** Claude and ChatGPT / Codex take any role the PI assigns. AIs never sign off closures or
  approve compute.
- Follow `RESEARCH_PROTOCOL.md` (version 2, adopted from ModelG).
  - **Exploratory track:** the PI's go-ahead in words, an `explore/LOG.md` entry before running, a 2 core-hour cap,
    and the label EXPLORATORY.
  - **Confirmatory track:** a frozen spec plus a GATES entry committed before production, a separate recomputation
    path, and PI closure.
- Study numbers: GC001, GC002, … for registered stages; EX-GC01, EX-GC02, … for explorations. Inherited ModelG labels
  (EX008–EX010, MG195) keep their names.
- A registered result is never re-scored after the data are seen. Deviations are recorded, not smoothed over.

## Working conventions
- Commit and push regularly with stage-prefixed messages (e.g. `GC001: ...`). **Pull before you push.** The PI also
  commits from VS Code.
- Every number in a record must be reproducible from committed files. Keep the repo small: commit records, small
  profiles and summaries. Bulky fields and pair states stay out of git with a committed SHA-256 manifest.
- All PNG figures go in `paper/figures/`.
- Write for the PI in plain language: the conclusion first, then the evidence. Say plainly when something failed.

## Compute
- Worker count comes from the environment variable `MODELG_WORKERS` (name kept from ModelG; default 2 to 4). On the
  PI's 10-core PC use at most 8. Set numerical-library threads to 1 per worker.
- The PI's PC: Windows, Python 3.13.2, numpy 2.2.4, scipy 1.16.1. Cloud runs use newer versions, so record the
  platform with results.

## Code map (src/)
- `ex008_longrange.py`: 3D radial Gc1 solver (`Model`, `particle`, `drift` = l = 1 mobility, `leading_eigs`).
- `ex009_pair2d.py`: axisymmetric (ρ, z) two-particle co-moving bordered-Newton solver (`Axi`); EX009 scripts.
- `ex010_network.py`: generic mass-action network solver (`Net`) for Gc2 (charge plus gravity).
- `mg195_*.py`: the MG195 confirmatory study (inverse-square), with its recomputation and figure.
- `fig_particle_*.py`: illustrative figures (Model G vs Model Gc particles).
- `exgcNN_*.py`: one script family per exploration EX-GC01 to EX-GC26 (run, summary, figure). Reusable pieces:
  `exgc07_sim3d.py` (3D spectral-in-angle solver), `exgc08_belt.py` (cube-rotation frames, lifts, three-node move
  table), `exgc09_lattice.py` (3D lattice Gc1, tokens), `exgc12_soliton.py` (lattice knot, degree),
  `exgc14_knotgc.py` (knot-gated Gc particle), `exgc16_recoil.py` (token recoil absorbers), `exgc20_wave.py` (wave-like
  G medium), `exgc23_transport.py` (+ C++ kernel; ChatGPT), `exgc24_impedance.py` (ChatGPT), `exgc25_framefree.py`, `exgc26e1.py` (+ `exgc26e1_stencil.cpp`: Nordström-weighted 26B2 stencil, B = 2 rational map, baryon centroid), `exgc26c.py` (knot + dynamical Nordström field on a larger box, common action), `exgc26d.py` (boosted knot with its own field), `exgc26f.py` (two knots attracting through their own field).
- `whitepaper_figures.py`: white-paper figures computed from committed records only.
- White paper: `paper/Model_Gc_white_paper_v0.4.md` (frozen 2026-11-22, through EX-GC26); v0.3 frozen through EX-GC25.
- Test: `python tests/test_regression.py` (seconds; add `--pair` for the 2D pair regression, about 1 minute).
