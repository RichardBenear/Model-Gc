# EX-GC35W — recoil with a persistent particle vertex

2026-11-25 | EXPLORATORY | ChatGPT | [Plan](../../spec/EX-GC35W_plan.md). PI closure pending.

**Recoil reintroduces order dependence even when the particle vertex persists.** A carrier on either adjacent site can be consumed and leave a signed deposit on P. The position relation is read for absorption. A recoil with enough net signed deposits consumes the current position and a threshold-sized group of same-sign deposits and moves one site toward that sign. Delayed initial carriers have one optional arrival event. The line has five sites and P starts at site 2. Each rule conserves the carrier accounting when deposits and spent deposits are included.

| Threshold | Fixtures | Complete histories | Fixtures with a recoil disabling an eligible absorption | Aggregate terminal and type-labelled DAG failures | Anchored graph failures |
|---:|---:|---:|---:|---:|---:|
| 1 | 1,470 | 1,063,497 | 948 | 948 | 948 |
| 2 | 1,470 | 270,957 | 304 | 304 | 552 |

The aggregate ledger counts free carriers by position and deposits by sign, irrespective of carrier identity. The event comparison includes positive dependencies on position, available carriers and consumed deposits. It does not encode negative eligibility checks. At threshold one, put two carriers on the left. Absorb both before moving and two right steps are possible. Move as soon as the first deposit arrives and the second carrier is left behind. At threshold two, three left carriers suffice: absorb all three before moving and one deposit remains, or move after two absorptions and leave the third carrier free. Thus failure includes aggregate outcomes and persists after imposing a larger threshold. The conflict is a read/write conflict on position, even though absorption does not consume P.

[Code](../../src/exgc35w.py) · [Data and minimal witnesses](results.json)
