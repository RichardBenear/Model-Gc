# EX-GC36C — does rule 66427's order-free particle persist in its own debris?

2026-11-25 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC36C_plan.md), with D1 (T cut to 600 after a smoke run). PI closure pending. Cloud container, Python 3.13. About 0.9 core-h, including the stopped first run and the smoke run.

**No. 66427's order-free particle is not persistent (Q1 fails).** Most runs die within a few events: the median lifetime is 3 from either phase. Only 5.5% (phase 0) and 15.5% (phase 1) of schedules are alive at 500 events. Every one of those survivors is no longer a single particle. Particle-shaped supports multiply into a swarm of 53–74 at once, and about 72% of events happen outside any live support. The ordered 55744 control lives longer (Q2 holds) and shows the same split between quick death and a proliferating swarm.

## Results (200 random schedules per start, T = 600)

| Rule, version, start | Alive at 500 | Median lifetime | Longest death | Survivors' peak live supports | Survivors' events inside a support |
|---|---:|---:|---:|---:|---:|
| 66427-C, phase 0 | 11 | 3 | 15 | 56–74 | 28% |
| 66427-C, phase 1 | 31 | 3 | 35 | 53–74 | 29% |
| 55744-L, phase A | 90 | 13 | 46 | 50–75 | 62% |
| 55744-L, phase B | 89 | 12 | 18 | 50–76 | 62% |
| 55744-L, phase C | 64 | 7 | 19 | 52–77 | 60% |
| 66427-L, phase 0 | 28 (27 halted) | 3 | 9 | 1–65 | 72% |
| 66427-L, phase 1 | 30 (26 halted) | 2 | 6 | 1–78 | 58% |

- **The outcome is bimodal.** Every run alive at 500 was also alive at 600. Every death happened by event 46. There is no gradual decay.
- **"Survival" means proliferation.** The tracker follows a lineage of renewals. In every surviving run of 66427-C and 55744-L, that lineage branches into dozens of particle-shaped supports. It is a growing medium whose motifs are particle-shaped, not one particle moving through debris.
- **The ordered 66427-L control** mostly "survives" by halting. In 53 of its 58 survivors no legal move is left, so the particle is frozen.

## Recomputation

A separate tracker replays 5 schedules per start for up to 100 events. It scans every candidate 5-set: all unions of two live relations' vertex sets that share exactly one vertex. Every phase contains such a pair. It agrees with the incremental tracker on all 35 schedules ([recompute.json](recompute.json)).

The first version of this check scanned only the 12 newest vertices. It missed supports on older vertices and disagreed on the long-lived runs. I replaced it before recording; the replacement is described in the [script](../../src/exgc36c_recompute.py).

## What it means (inference)

- The only order-free closed particle in CRR's catalogue does not survive its own debris as a single object. Isolated closure (36B) was necessary but far from sufficient.
- In both 66427-C and 55744-L, the long-lived outcome is a proliferating swarm of particle-shaped motifs. This fits CRR's own reading that particle and vacuum are the same machinery. It also means "one particle you can track, weigh or accelerate" is not what these rules produce from an isolated start.
- Together with 36A/36B, the 80-rule CRR catalogue now gives:
  - no order-free particle that persists;
  - no particle whose internal forks rejoin within the tested windows.

  For the programme's goal, the catalogue is exhausted. The remaining route is a search designed for the target from the start: order-free closure, persistence in its own debris, and fork reconvergence as selection criteria together.

## Files

[Code](../../src/exgc36c.py) · [Recompute](../../src/exgc36c_recompute.py) · [results.json](results.json) · [recompute.json](recompute.json)
