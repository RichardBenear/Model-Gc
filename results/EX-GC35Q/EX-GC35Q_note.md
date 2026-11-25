# EX-GC35Q — whole-element carrier rules and critical pairs

2026-11-24 | EXPLORATORY | ChatGPT | [Registered plan](../../spec/EX-GC35Q_plan.md) · [Claude review](../../review/EX-GC35P_claude_review.md). PI closure pending.

**Only one-species exclusion passed the full finite screen.** The pass-through and nonrecoiling absorber examples in the review pass as plain causal graphs on their highlighted overlaps. Both fail in larger contexts.

Every matched element is consumed, and each output is a new version. An edge links a producer event to a later event that consumes at least one of its outputs. The labelled view retains the current match index and rule. The plain view keeps only event dependencies. Rule B also allows a single whole-graph reflection with R/L exchange if it preserves the initial word. All legal schedules were completed on every open word of lengths 1–6. A–C terminate by decreasing inversion count. In D/E absorption reduces carrier count and intervening hops strictly advance carriers to the right. These are termination arguments, separate from the graph comparisons.

| Rule | Starts | Terminal histories | Terminal-word failures | Labelled DAG failures | Plain DAG failures | Minimal overlap words |
|---|---:|---:|---:|---:|---:|---|
| A · One-species exclusion | 126 | 414 | 0 | 0 | 0 | None |
| B · Pass-through | 1,092 | 27,988 | 0 | 294 | 278 | R0L |
| C · Blocked R/L | 1,092 | 2,363 | 183 | 183 | 112 | R0L |
| D · Absorber | 1,092 | 4,022 | 122 | 446 | 261 | GPG, PG0, PGP |
| E · Recoil absorber | 1,092 | 2,384 | 201 | 274 | 155 | PG0, PGP0 |

For B the reflection allowance reduces the labelled failure count from 294 to 290. The paired-pool negative control F reproduced the labelled and unlabelled failure at (2,2), with canonical same-site and neighboring pool fixtures in the raw data. The representative pure three-token same-site pool has no distinguishable producer edge before any arrival. The pool control uses GC35L/35P's carrier-producer graph convention and is separate from the word-overlap enumeration.

## Minimal critical pairs and larger witnesses

All minimal overlap words and both completions were scored. Each of these overlap branches had one terminal completion.

- **B, R0L.** Both branches reach L0R with the same three-event triangle. The labelled mismatch is removed by reflection and R/L exchange. This local success does not extend to the extra leading vacancy in **0R0L**. Both branches end at L00R. One has edges 0→1, 0→2, 1→2, 2→3, and the other has edges 0→1, 0→3, 1→2, 1→3. Their longest directed paths contain four and three events. This is the smallest plain-graph counterexample.
- **C, R0L.** The two one-event graphs agree when labels are erased, but they jam at different words, 0RL and RL0. The smallest plain-graph difference is R00L, where one completion has a dependency edge and the other has two independent events.
- **D, GPG and PGP.** Each gives a plain-equivalent local comparison. PG0 fails terminal agreement because the carrier can be absorbed to leave P0 or hop away to leave P0G. The first plain graph failure is G0PG. After the first hop, the two absorptions may form a three-event chain or a pair of independent events feeding the last absorption.
- **E, PG0 and PGP0.** Both fail terminal agreement and have plain-equivalent local graphs. The first plain difference is GPG0. Both branches end at 00P, with two events in one completion and three in the other.

D uses the literal length-reducing absorptions from the review and one-way vacancy hops. E implements a left-side recoil only when the right destination is empty, with the explicit rule GP0→00P, and retains the right-side contraction PG→P. These are two specific absorber variants. Neither result establishes a general obstruction to absorber rules.

The terminating string overlaps are sufficient here to detect particular local failures and to support ordinary state-confluence reasoning for B. They are insufficient for graph equivalence, as the 0R0L example directly demonstrates. The independent A control consumes whole elements and has no competing match. No tested interacting rule passed the full graph screen. This study addresses Step 1 only.

Recorded CPU time was 3.64 seconds. [Code](../../src/exgc35q.py) · [Counts, terminal histories for the controls, and witness graphs](results.json)
