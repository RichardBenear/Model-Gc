# EX-GC36A — causal-invariance audit of CRR Rule 55744 (stages 1–3 of the EX-GC36 design)

2026-11-25 | EXPLORATORY | Claude | PI go-ahead in thread ("Go ahead", 2026-11-25 16:26 UTC) on the [EX-GC36 design](../docs/EX-GC36_rule55744_causal_invariance_design.md). Default cap 2 core-hours; expected ≤0.5. Written before any test below is run.

## Inputs

These are CRR files in `external/crr/`, unchanged, with SHA-256s recorded in the design note. Added here: `W21_CLOSED_PARTICLE_HITS.json` (023ecc95…6a742). Reading its 55744 entry gives three particle phases on five vertices: A = phase 0 (10 relations), B = phase 1 (10), C = phase 2 (11). B has two internal gates, to C and to A. A and C have one each. This decode is input reading, not a test.

## Engine (new code `src/exgc36a.py`, written for this study)

- A generic multiway engine for hypergraph rewriting with explicit relation instances and producer provenance. Each event records the relation instances it consumes, the ones it only reads (the 55744 gate) and the ones it produces.
- **Rule 55744, version L:** exactly as in `x6_generic.py`. That covers five distinct vertices, u < v by vertex number, and fresh w = next number along the path. An output relation already present is not added again (set semantics). Vertex numbering follows each path's own counter, as the simulator does.
- **Version C:** the same, but with no u < v condition. Both orientations of the {u, v} pair are separate matches. This removes the only place where vertex numbers enter the rule.
- **Causal graph:** an edge i→j when event j consumes *or reads* a relation instance produced by event i.
- **Join test (35V style, extended):**
  - At a fork (state S, events e1 ≠ e2), extend each branch by equal numbers of further events, up to k = 3.
  - The pair *joins* at the first level where some left history and some right history give isomorphic combined graphs. A combined graph is the events, the dependency edges and the live relations with their producers and vertices.
  - Results are *join*, *unresolved* (bound reached) or *cap* (the history cap of 20,000 per side was hit).
- **Two conventions:**
  - **(U) unlabelled:** no vertex names and no input anchors. This is the literature definition.
  - **(A) anchored:** seed vertices are labelled, and initial relation instances are labelled roots.
- For a ternary relation T(a;{b,c}), the source edge and the two unordered target edges are distinguished.

## Stages and predictions

**Stage 1 — calibration.**
- *Positive controls:* the published rules from the 35V controls.
  - `xy→xy,yz` on the two-loop seed.
  - `xy,zy→xw,yw,zw` on two asymmetric seeds: a three-in star with one extra edge on one leaf, and a four-in star.
  - Noninjective matching as in the Wolfram models.
- *Negative control:* `xy,xz→xy` (keep one of two out-edges) on a two-out star with one extra edge on a leaf, (0,1),(0,2),(2,3). The two branches end in non-isomorphic states, two disjoint edges against a path, so the pair must fail. *Amended before any run:* the symmetric two-out star first written here gives isomorphic end states under U and so cannot serve as a negative control.
- **P1 (gate):**
  - Under U, `xy,zy→xw,yw,zw` joins every fork on both asymmetric seeds within k = 3.
  - The negative control fails under both conventions.
  - If either part fails, the test is not calibrated, and stages 2–3 are reported as descriptive only.
- **P1b:** `xy→xy,yz` joins under both conventions.

**Stage 2 — vacuum.**
- (2a) From the CRR start state, enumerate the multiway system to depth 4. Run the join test on every distinct fork at every state reached, under L and C and conventions U and A.
- (2b) Run the CRR simulator to 2,000 events with seed 0. Take the first 50 gate pairs that share a consumed relation, plus the first 50 read/write pairs. Extract the relations whose vertices all lie within 2-section distance 1 of the pair's vertices, and run the join test on that extract.
- **P2:** under L, at least one 2a fork fails to join under U. This is the birth-order prediction.
- **P3:** C under 2a, and 2b in all versions: open. Report the fraction joined.

**Stage 3 — particle.**
- The seed is phase B, isolated (its 10 relations on 5 vertices). The fork is its two internal gates.
- First, check CRR's closure claim: each branch's successor support is a particle phase (C or A).
- Then run the join test, under L and C and both conventions.
- **P4a:** both branches are particle phases. This reproduces CRR.
- **P4b:** under A, the fork fails. The first events consume different relations, as in 35R.
- **P4c:** under U: open. This is the main question of the study.

## Not in scope

Stage 0 (reproducing the W37 anchor; PC), Stage 4 (the Gc reading) and any claim beyond the bounded tests. A join within k = 3 is evidence of resolution, not a proof. An unresolved fork at k = 3 is not proof of failure, unless both branches are terminal.

## Outputs

`results/EX-GC36A/` (note plus `results.json`). Platform recorded. Figures, if any, in `paper/figures/`.

## Deviations (recorded before Stage 2b was run)

- **D1.** The first all-stage run reached the 50-minute command limit during Stage 2a version C, after about 0.83 core-h, and wrote no file. Stages are now run separately, with the same code and parameters. Stages 1, 3 and 2a-L were rerun this way.
- **D2.** In a grown vacuum, hub vertices (2-section degree up to about 600 at 2,000 events) make some radius-1 extracts very large. Stage 2b keeps the registered pairs and k = 3, but skips any pair whose extract has more than 400 relations. Skipped pairs are counted and listed. I chose the threshold before seeing any 2b result.
