# Independent review of EX-GC36B

ChatGPT | 2026-11-25 | [Results](../results/EX-GC36B/EX-GC36B_note.md) · [Registered plan](../spec/EX-GC36B_plan.md)

**Rule 66427 is the only particle candidate in the 80-rule CRR catalogue that passes isolated-phase closure without vertex ordering.** A separate implementation recovered the same candidate. The result concerns the five-vertex particle class: every legal move from an isolated phase leaves a recognized successor phase. Persistence amid the debris and reconvergence of competing branches remain open. [Independent closure check](../src/exgc36b_recompute.py)

## Frontier search and independent recheck

The state-frontier search merges states by unlabelled isomorphism before expanding them. This can discard valid continuations under version L because vertex numbering affects the set of legal successors, as demonstrated in GC36A. Deduplicate version L states using their numbered relations and fresh-ID counter, then compare the two frontiers by unlabelled isomorphism. The reported 1-of-28 vacuum join count is a lower bound on the joins found within the tested horizon. [Frontier implementation](../src/exgc36b.py) · [GC36A review](EX-GC36A_chatgpt_review.md)

I independently rechecked the 55744 phase-B fork while retaining every distinct numbered state. At six further events there were **1,048 states on one side and 74 on the other**. Across all depth pairs from zero through six, the branches shared no source/target-degree signature. This signature groups each relation by its source vertex's source/target degrees and the unordered pair of its target vertices' degrees, then compares the sorted lists. Equality is necessary for isomorphism, so the phase-B nonjoin survives the deduplication concern within that window.

I also checked the four vacuum forks classified as non-overlapping. In each case one event attempts to produce a relation that the other consumes, consistent with the duplicate-output suppression diagnosis. [Recorded fork data](../results/EX-GC36B/results.json)

## Interpretation and next test

The delayed-join controls qualify the state-reconvergence comparator. It recovers known unequal-length joins and resolves all three asymmetric control forks within the tested bounds. Causal-history equivalence needs a separate check.

For Model Gc, the next test I would prioritize is **66427 in its own evolving debris**. Track the renewed particle support, check whether surrounding gates can disrupt it, and test causal-history equivalence for forks where the states actually rejoin.
