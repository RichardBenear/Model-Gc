# A causal-invariance test with competing carrier events

2026-11-24 | Design only | ChatGPT | [GC35G–I review](../review/EX-GC35GHI_claude_review.md)

## What to compare

Require exact equivalence of the completed event dependency graphs up to relabelling of indistinguishable tokens. Preserve the initial marked sites, each transport channel, any boundary labels, and which events produced the inputs consumed by each later event. Distinguish genuine graph equivalence from agreement of the final site populations. Equal populations or transport observables would be valuable as a weaker result, with their own label. They cannot establish equality of causal histories.

## A branching control and a candidate

Start with a finite one-dimensional chain and two absorbing endpoints. Keep sink arrivals in the carrier ledger so total carrier number is conserved. Give a site one available token and allow either adjacent hop to consume it. A pair of simultaneously eligible hops genuinely conflicts and only one can use that token. A minimal branch with destinations leading to different labeled sinks should fail even the endpoint comparison. Keep the fixed-version GC35F graph as the scheduler bookkeeping control.

Then try paired redistribution. A site with at least two available tokens may consume any two and send one to each neighbor. There is no input reservation or canonical identity tie-break in the causal-invariance test. A token can appear in several legal matches and can be consumed only once. Include both site-order forks and same-site forks in which the two candidate pairs share one token but come from different producer histories.

On a finite line with sinks, the weighted interior token count with weights wᵢ=i(L+1−i) decreases by two per firing. This supplies a finite search bound. The rule resembles an Abelian chip-firing system at the level of site counts. The question here concerns the finer graph obtained by retaining producer-event links. Matching site counts or matching numbers of firings at each site does not decide it.

## Proposed readouts

Enumerate reachable branch completions on short marked lines with small token counts. Compare terminal site and sink ledgers, the firing count per site, and canonical directed dependency graphs preserving site and channel labels. Record the smallest inequivalent completed fork if one exists, including an explicit pair of legal schedules. If all enumerated cases agree, record the finite domain and identify what would be needed for a proof beyond it. Only after the small conflict controls work would a larger carrier model or metric-dependent rates make the test interpretable.

This is a design for a future experiment. No search or causal equivalence result is claimed.
