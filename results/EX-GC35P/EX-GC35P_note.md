# EX-GC35P — carrier catalogue and reciprocal clocks

2026-11-24 | EXPLORATORY | ChatGPT | [Registered definitions and predictions](../../spec/EX-GC35P_plan.md). **PI CLOSED 2026-11-24.**

**Only the independent-hop control passed both graph comparisons.** All five rules conserve total carrier number including labelled sinks. The three new rules with shared-input or batch choices each have graph counterexamples even after removing site labels.

| Rule | Starts | Terminal ledger | Labelled graph | Site-unlabelled graph | First witness |
|---|---:|---|---|---|---|
| A · Independent fixed-direction hops | 8,517 | Pass | Pass | Pass | Each carrier has its own fixed chain |
| B · Paired redistribution | 494 | Pass | Fail | Fail | Two sites with counts (2,2) |
| C · Whole-site split | 494 | Fail | Fail | Fail | Graph at (2,2), ledger at (2,3) |
| D · Consume a whole directional FIFO queue | 8,517 | Pass | Fail | Fail | One rightbound token at each of two sites |
| E · Ballistic hops with optional opposite-direction capture | 8,517 | Fail | Fail | Fail | One token in each direction at a single site |

For A and D fixed directions determine the entire directional sink ledger. A additionally determines a disjoint union of fixed-length carrier chains with fixed site histories. Two extreme orders were run for every A configuration and the exact graph scorer compared both orders in all 294 configurations with up to three carriers. For B, C and E, terminal ledger possibilities were computed by memoized exhaustive search over population states. The graph searches traversed histories in order of total token count and site count until each first witness, using exact labelled multigraph isomorphism. Site-unlabelled graphs retain event kinds, direction and edge multiplicity. No count of all graph mismatches is inferred from the first-witness searches.

C can consume an arrival in an existing batch or leave it for a later batch. At (2,3), the terminal populations include both (1,1,0,3) and (2,0,1,2), listing the two sinks at the ends. D similarly can send two tokens in one event or in separate events, producing a different dependency graph while preserving directional sink totals. E can capture the same two opposite-moving inputs or permit one to depart first. All three effects change actual event structure. Deterministic outputs per selected match do not eliminate choices about which inputs have arrived or remain available.

The relevance of the two literature results is distinct. Bond and Levine's local-to-global abelianness concerns network processing results under its hypotheses. Gorard's Wolfram-model criterion compares causal graphs across update histories. The ledger and graph results here are computed independently. [Bond–Levine](https://arxiv.org/abs/1309.3445) · [Gorard](https://arxiv.org/abs/2004.14810)

## Passing control in the GC35O loop

Three independent directed carriers were added to the reciprocal one-dimensional packet and lapse model. Two carriers traverse the same rightbound path and the third traverses the line leftward. Counts source the lapse together with packet energy, and the lapse at arrival determines each next carrier deadline. The field scheduler uses neighbor version deadlines and a minimum travel-time bound along known carrier paths to rule out an unread arrival or departure before a site's sampling time. Carrier events wait until the destination lapse history is complete through the firing timestamp.

All **16** seeded randomized runs matched the chronological reference across packet and lapse values, event timestamps, carrier ledgers and labelled version-dependency records. All 16 overlapped nonadjacent event times. The final sinks received one and two tokens. At proper step .02, the last hop occurred at 1.362113, versus 1.360000 when clock and packet-metric feedback were disabled. The final lapse ranged from 0.91366 to 1.10507. No lapse clamp was reached. The packet difference with feedback disabled was 0.5325 in the unscaled spatial L² norm.

This qualifies independent overlapping paths in the reduced clock host, with every carrier maintaining its own inputs. The passing rule does not have destructive competing matches. No candidate in this catalogue passed the graph comparison with such competition, and no claim is made about other rule families or full tensor Y1.

[Catalogue code](../../src/exgc35p.py) · [Catalogue data](catalogue.json) · [Clock code](../../src/exgc35p_clock.py) · [Clock data](reciprocal.json)
