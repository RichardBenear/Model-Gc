# EX-GC35U — producible input borders on strings

2026-11-25 | EXPLORATORY | ChatGPT | [Plan](../../spec/EX-GC35U_plan.md). PI CLOSED 2026-11-26.

**Both registered predictions held.** Rechecking the same 734 nonexpanding singleton and two-rule systems left zero binary and two ternary passes on the revised local border. The four earlier plain-graph survivors through length seven were unchanged. The pass for the ternary pair (01→1,12→1) and its mirror arises because only symbol 1 can be produced. The fixed 0 and 2 input ports that distinguished the strict-border histories can never be outputs of this rule set. The binary pair (01→1,001→01) produces both symbols and fails the revised border as before.

With a new vacancy symbol, adding the single transport rule `0.→.0` gives the witness `0.12`. Both branches end at `.1`, while one has a three-event chain and the other has two independent producers feeding the last event. Separately adding `.2→2.` gives the reflected transport failure starting at `01.2`, ending at `1.`. These are plain graph counterexamples, with input labels erased. The revised certificate recovers the undeveloped absorption example and does not transfer that result to a transport-fed absorber. The short-string search is complete. [Code](../../src/exgc35u.py) · [Data](results.json)
