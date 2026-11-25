# EX-GC35P — carrier rules with competing inputs

2026-11-24 | EXPLORATORY | ChatGPT | PI authorized implementation of the GC35L–O review | Cap .1 core-hours for the catalogue.

Use labelled absorbing lines of 1–4 interior sites and 0–7 initial tokens. Reuse the GC35L enumeration convention and exact labelled and site-unlabelled dependency graph comparisons. Keep node/edge kinds, direction and multiplicity in the unlabelled view. Measure terminal carrier ledgers and total carrier conservation. Predictions, in order ledger, labelled graph, unlabelled graph:

- A, independent single-token hops with immutable direction: pass/pass/pass. Every token advances monotonically toward its sink.
- B, GC35L paired redistribution with subset selection: pass/fail/fail.
- C, whole-site threshold-two firing with lower half sent left and remainder right: fail/fail/fail. Arrivals before firing can change batch size and odd-token allocation.
- D, forward per-link FIFO queues, consuming the entire nonempty queue per event: pass/fail/fail. Final exits agree, but batches can merge depending on readiness.
- E, single-token ballistic hops plus a legal local opposite-direction capture that consumes one token of each kind and deposits two inert tokens: fail/fail/fail. Departures can bypass capture.

Flags are finite-domain screens, with smallest witnesses by token count then site count. Directional rules test all left/right starting multisets up to the cap. Attach a passing rule to the reduced GC35O host and compare packet/lapse fields, timing, ledgers and dependency graphs under randomized legal schedules. For the no-overlap control use stable token IDs and local deadlines with fields progressed through each read time. Boundary labels and initial roots remain fixed; distinguish proof from finite testing. No full tensor gravity claim.