# EX-GC10 (EXPLORATORY): route 3, stage 3 – a Gc particle tethered by its own G supply

Claude, 2026-11-19, at the PI's request ("Go ahead with the next test. If you exceed the 2 core-hours, you can use
my local cores for 40 core-hours."). Exploratory track; Gc1, the tether rules and the coupling are our construction,
not LaViolette's.

## Answer
**Yes, in the deterministic-chemistry runs: a particle turned once relative to its G-supply tether keeps a twist that
never relaxes, and it pays for it with about 8% less G consumption. Turned twice, the tether untangles and the
consumption comes back.** The Z₂ index (−1 for one turn, +1 for two) never changed in any run with the no-cut rule.
With token chemistry (Ω = 3000) the topology is the same, but the mass signal is weaker and noisier: one turn is
lowest, while two turns stay partly tangled and fall in between.

## Set-up (design choices, ours)
- Frames (24 cube rotations) on the EX-GC09 lattice nodes (h = 0.5, box ±6, 25³ nodes). EX-GC08's three-node lift
  rule, with steps of 90° and 120°, is the no-cut rule.
- **Tether = G supply.** Each node attempts frame moves at a rate proportional to the G traffic into it: at most one
  attempt per sweep, and a sweep every 5 chemistry steps, about 133 per time unit in the far field. Where G is
  depleted (the core region), the tether moves more slowly.
- **Twist impedes supply.** A G hop across an edge whose two frames differ succeeds with probability 1 − ε (ε = 0.5).
  That is the only way frames act back on the chemistry.
- Frame moves are Metropolis on the number of mismatched edges at a fixed temperature T_f = 0.6. This is a dynamics,
  not annealing. T_f was picked in a frame-only calibration: at 0.6, two turns untangle and one turn sticks; at 1.0
  both stay disordered; at 1.5 everything melts.
- **"Rotating the particle by hand":** start from the settled particle and spread n turns of twist about z radially
  between the core (|x|∞ ≤ 1, frames fixed) and the reservoir face (fixed), as in EX-GC08.

## Results
| Chemistry | Turns | Runs | Twist energy, late (from) | Index | M late / untwisted |
|---|---|---|---|---|---|
| deterministic | 0 | 1 | ≈ 50–100 (thermal) | +1 constant | 1 (5.878) |
| deterministic | 1 | 3 | 1036–1203 (5880) | −1 constant | **0.921** (5.413, 5.414, 5.421) |
| deterministic | 2 | 3 | 108, 476, 483 (12240) | +1 constant | 1.000, 0.994, 0.994 |
| deterministic, ε = 0 (control) | 1 | 1 | 962 | −1 constant | 1.000: no coupling, no mass shift |
| deterministic, pairwise rule (control) | 1 | 1 | 78 | flipped to +1 | 1.000: the twist leaks away |
| tokens Ω = 3000, T = 60 | 0 / 1 / 2 | 1 each | 108 / 1043 / 943 | constant | 5.20 / 5.14 / 5.69 (late means) |
| tokens Ω = 3000, T = 150 (PC) | 0 / 1 / 2 / 2 | 1 / 1 / 2 | 108 / 1228 / 468, 468 | constant | 4.87 / 4.55 (−6.7%) / 4.62, 4.70 (−5%, −4%) |

- **Where the twist sits.** The leftover one-turn twist collects in a shell around the core: about 90% of its
  mismatched edges lie in the first three lattice layers outside it. This is the tightest place it can be squeezed,
  and it is where G flows inward, so it chokes supply most there. The two-turn leftover (when it does not fully clear)
  is a looser cluster in layers 3–7, where low G traffic makes the frames slow.
- **A two-turn tether does not always clear completely.** It cleared to the thermal level in 1 of 3 deterministic
  runs and in 0 of 3 token runs by the end; the rest held a residue of about 470 edges for 50–100 time units. It
  carries index +1, so it is not topologically protected. It is a kinetic trap, probably made worse by the coupling
  itself (twisted edges carry less G, so they get fewer move attempts). In the deterministic runs it costs only
  0.6% in M.
- **The token runs have a shared drift.** In every token run, including the untwisted one, M falls by about 15–18%
  after t ≈ 40–50. At the same time the core's centre moves about 0.3–0.45 off the lattice node it started on.
  Coarse-lattice M depends on where the core sits relative to the lattice (EX-GC09: 10% between h = 0.5 and 0.4), so
  this looks like a lattice artefact common to all runs. It makes the late token comparison noisy. Also, the core
  frames stay fixed at the box centre and do not follow the drifting core.

## What it means for route 3
- The full chain works in a model where the 720° property comes from structure, not tuned dynamics. A particle's own
  G supply acts as the tether, a local no-cut rule protects a one-turn twist, two turns undo it, and the twist has a
  physical cost (less G consumed, so less gravitational mass while twisted).
- The "mass" of the twisted state depends on ε, which we chose. The finding is that a coupling of this kind makes the
  360° and 720° states physically different. The size of the difference (8% at ε = 0.5) is not a prediction.
- **Not yet shown:** the exchange sign (two particles swapping places compared with a 360° turn of one of them), and
  any link to a quantum sign. A real rotating particle, rather than a twist imposed by hand, is also still open.

## Method and records
- `src/exgc10_tether.py`: `Frames` (checkerboard Metropolis under the move table from `exgc08_belt.move_table`;
  `by_layer`; `z2` = `exgc08_belt.z2_lift`), `wlap` (weighted G Laplacian), `traffic`, `hop_w` (token G hops with
  edge pass probability), `tok_step` (EX-GC09's Euler-consistent token update), `run`, `frames_only`. The settled
  particle is cached as `settled_h0.5_L6.npy` (not committed; it is regenerated deterministically by `settled`).
- `src/exgc10_summary.py` → `summary.json`, `paper/figures/EX-GC10_tether.png` (reads the records only).
- Records: `frames_*.json`, `det_*.json`, `tok_*.json`.
- Platforms:
  - Cloud: Linux, 1 core, Python 3.13.18, numpy 2.5.3. Frame calibration, all deterministic runs, and the seed-0
    token runs (T = 60).
  - Richard's PC: Windows 11, Python 3.13.2, numpy 2.2.4, 1 thread per run. Token runs at T = 150, seeds 1 and 2, at
    commit 02380b1; pushed on branch `exgc10-pc-runs` (7cf467c) and copied into main.
- Compute: about 1.85 core-hours in the cloud and about 5.2 on the PC. The PC figure includes 1.9 core-hours for a
  first set of runs that was stopped and restarted to add the core-position diagnostic. Within the 2 + 40 core-hour
  allowance.

## Deviations and limits
1. **Frame dynamics speed and temperature are free choices** (ν = 1, a sweep every 5 steps, T_f = 0.6), set by the
   frame-only calibration before the chemistry runs. Lower T_f freezes the tether and higher T_f melts it; only a
   window near 0.6 separates one turn from two within the run times.
2. **The twist is imposed by hand**, as a radial spread, not produced by rotating the particle. The core frames are
   held fixed at the box centre.
3. **Seeds:** 3 deterministic runs for 1 and 2 turns; tokens 2 runs for 0 and 1 turns, 3 for 2 turns, of different
   lengths (T = 60 in the cloud, 150 on the PC). The seed-0 token runs lack the core-position diagnostic, which was
   added afterwards.
4. **The token mass comparison is confounded** by the shared drift described above. Only the deterministic runs
   support the 8% figure.
5. The logged prediction "2 turns untangle to 0" held fully in 1 of 6 runs with chemistry. In the rest a
   non-topological residue remained (see Results).
