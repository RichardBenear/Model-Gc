# EX-GC35L — competing consumptions on short marked lines

2026-11-24 | EXPLORATORY | ChatGPT | [Registered predictions](../../spec/EX-GC35L_plan.md). PI closure pending.

**Ledger and firing counts agreed in all 494 initial configurations. Both causal graph definitions failed.** Enumeration traversed 470,058 terminal schedules and token-producer choices on one to four interior sites with zero to seven initial tokens. Equal producer copies at one site were treated as indistinguishable. Exact labelled multigraph isomorphism distinguishes carrier dependency edges from predecessor-at-site edges and preserves labelled sites and sinks.

The ledger agreement was predicted by the chip-firing rule. The one-token branching control instead sends its token to different labelled sinks, and fails the ledger comparison. The token graph differs in 294 initial configurations. The site-history graph, with FIFO interchangeable arrivals and predecessor links, differs in 176 configurations. A different coarse-graining of arrival dependence would define a different graph test.

The smallest token counterexample is the two-site start (2,2). Fire the left site, then let the right site consume either its two original tokens or its new arrival and one original token. The branches share one original token. Both end with one token at each site and each sink and fire each site once. Only the second history has a carrier edge from the first firing to the second firing. The targeted shared-input fork starting from (2,3,2) also gives unequal completed token graphs.

The smallest FIFO site-history counterexample is the three-site start (0,5,0). Fire the center twice. Each edge site can then fire to return one token to the center. If the center fires as soon as the first return arrives, its third firing depends on the left or right producer according to which edge fired first. Both branches give terminal populations (1,1,1,1,1), including the sinks, and firing counts (1,3,1). The labelled dependency graphs differ even when tokens at each site are interchangeable under FIFO.

The data establish separation of the three readouts on this finite domain. They do not support causal-history invariance for either tested graph. [Code](../../src/exgc35l.py) · [Raw counts and witnesses](results.json)
