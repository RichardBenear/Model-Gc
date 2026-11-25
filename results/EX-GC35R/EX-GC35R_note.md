# EX-GC35R — systematic bordered search on short strings

2026-11-24 | EXPLORATORY | ChatGPT | [Plan](../../spec/EX-GC35R_plan.md). PI closure pending.

**Four competing rule sets pass the output-border and length-seven plain-graph screens. None pass the anchored input/output border.** This distinguishes an output-only comparison from a local certificate that can be glued into a fixed surrounding history.

The search enumerated every singleton and unordered pair in the registered slice of carrier-conserving nonexpanding rules. Binary sides have lengths 1–3 and ternary sides 1–2. Each eligible set has a minimal overlap containing the counted symbol 1. All completions of both sides of each minimal overlap are compared. The output border fixes each final element's position and producer dependency. The strict border also fixes each initial element's position and consumption dependency. A reachable repeated word is scored as a cycle. Only the output-border survivors continue to the exhaustive finite plain-graph check, stopping each failed system at its first witness.

| Alphabet | Elementary rules | Eligible sets | Local cycles | Output-border survivors | Strict survivors | Plain through length 7 |
|---|---:|---:|---:|---:|---:|---:|
| Binary | 31 | 470 | 37 | 12 | 0 | 2 |
| Ternary | 38 | 264 | 4 | 10 | 0 | 2 |

No history budget was hit. The four survivors are the following two pairs and their left-right reflections.

- Binary: 01→1 and 001→01. Both matches consume the carrier in 001. The two branches have the same visible output producer graph, while the first event consumes different anchored inputs.
- Ternary: 01→1 and 12→1. At 012 either side of the central 1 goes first. Both completions are a two-event chain and terminate at 1. The first event touches different anchored borders.

The four are genuine overlapping matches in the screened finite domain. They delete other symbols while conserving the distinguished carrier. Neither the finite graph agreement nor the strict border failure establishes a theorem for all contexts. Eighteen other sets pass the output-border screen and fail the plain-graph extension. This result limits the tested slice and leaves expanding systems, larger rule sets, and a different representation open. [Code](../../src/exgc35r.py) · [Counts and witnesses](results.json)
