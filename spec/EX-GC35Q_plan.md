# EX-GC35Q — whole-element carrier catalogue and critical pairs

2026-11-24 | EXPLORATORY | ChatGPT | PI authorized Step 1 of [Claude's GC35P review](../review/EX-GC35P_claude_review.md). Compute cap 0.1 core-hours.

Use open finite words of length 1–6 and no wraparound. Every rewrite consumes all matched elements and creates new output versions. Test exact terminal word agreement, event-producer DAG isomorphism with rule and current match position labels, and plain event DAG isomorphism with labels erased. For B separately allow the reflection plus R/L swap that preserves the starting word. Minimal overlaps are any two distinct intersecting matches whose union covers the minimal word, including containment if present. Exhaust all terminal completions after either first rewrite, then screen all starting words. Local overlap agreement is evidence only; terminal state confluence and causal graph invariance are separate.

Prospective predictions and rule choices:
- A: 10→01. Pass; no destructive overlaps.
- B: R0→0R, 0L→L0, RL→LR. Predict plain-graph pass and labelled pass modulo permitted reflection.
- C: same as B with RL blocked. Predict terminal and graph failure.
- D: literal contractions GP→P and PG→P, plus G0→0G. Predict plain-graph pass. The hopping direction is fixed to give terminating open-line runs.
- E: GP0→00P for a left absorption with right recoil, PG→P for right absorption, and G0→0G. No prediction. Recoil requires an empty right destination.
- F: reuse the GC35L paired-rule catalogue as negative control, with a fresh smallest-fork replay and exact scorer. Minimal pool overlaps include distinct subset pairs at one site and both orders of neighbouring firings; score canonical representatives separately from the site-element string overlap check.

The absorber rules consume carriers; the exclusion rules conserve each species. A finite screen and individual overlap examples alone cannot establish all-size graph invariance. Record smallest witnesses and computation cost.