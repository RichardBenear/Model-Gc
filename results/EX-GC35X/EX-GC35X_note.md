# EX-GC35X — background-gated probe in exclusion

2026-11-25 | EXPLORATORY | ChatGPT | [Plan](../../spec/EX-GC35X_plan.md). PI CLOSED 2026-11-26.

**Equal speed in the probe's own gated ticks is built into the rule. The probe does not slow in proportion to the fixed-site background rate.** The independent probe is a read-only passenger on a 32- or 64-site exclusion ring. During each causal layer it may move one site right exactly when there is a background event at its current site. Eligibility is assessed at the start of each layer, so coincident events cannot drive two successive hops. Average over all initial probe sites on the common measurement window. Compare to a marked ordinary carrier, and repeat with the reverse-depth foliation.

All eight randomized schedules agree for each of eight density, size and arrangement fixtures. Balanced rings have fixed-site rates .25 and .50 at occupancies 1/8 and 1/4. The independent probe's speeds per depth layer are .9453 and .9844 at those densities on both ring sizes and both tested foliations. The ratios are about 1.04 for the probe and 2 for the fixed-site clocks. The marked ordinary carrier moves at one site per depth layer. Shuffled starts give speeds between .886 and .962 and direction-of-depth comparisons still agree in these uniform fixtures.

The gated probe follows moving regions of activity. A rate measured at fixed sites does not describe the sequence it encounters while moving. Its one-site-per-gated-tick result therefore neither establishes a lapse-dependent slowdown nor resolves GC35S's inhomogeneous tests. This was a uniform dilute diagnostic, with no absorber, spatial curvature or imposed deadline.

[Code](../../src/exgc35x.py) · [Data](results.json)
