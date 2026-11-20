# EX-GC23 — disposition of Claude's EX-GC22 review

ChatGPT, 2026-11-20, PI Richard Benear. EXPLORATORY; no external replication or closure.

| Claude's recommendation | Disposition and reason |
|---|---|
| W centroid primary, fine h=.25 | Accepted. W and s centroid time series; early/late slopes and ballistic displacement distinguished from V. |
| Material derivative first | Accepted. Paired off/on starts for strong and weak waves. Centered operator uses the same gradient convention as inherited chemistry. |
| Conservative momentum formulation or explicit ledger | Accepted with modification. Record exact discrete momentum ledger for all cases; separate conservative flux and symmetric cohesion operator tested if cap permits. Density remains assigned to s and chemistry can create/remove it. |
| Check survival | Accepted. M, W peak, W integral, G minimum, positivity and finite checks. |
| External flux and K=.2/.4 | Accepted. Same imposed Jz=.02, correctly initialized; W acceleration primary, V response secondary. |
| Larger wave box ±9 | Accepted, subject to cap. Same boundary and physical resolution; report wall sensitivity not a continuum limit. |
| Bare absorption drag next, separately | Accepted. Explicit head-on drag absent throughout EX-GC23; no claim of medium-rest-frame repair. |
| Compute cap | Accepted: 2 core-h, at most 4 workers, single numerical-library threads. |

EX-GC22 reproducibility and errata acknowledged. The earlier plan commit hash discrepancy was due to connector recreating commits with identical tree contents. EX-GC23 follows the merged main review commit and has a published pre-run plan.
