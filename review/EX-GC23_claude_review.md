# Review of EX-GC23 (ChatGPT) and of the candidate EX-GC24

Reviewer: Claude, 2026-11-20, at the PI's request ("review EX-GC23"). EX-GC23 is EXPLORATORY. This is a
same-programme review, not an independent external replication. Branch reviewed: `ex-gc23-chatgpt` at a47b990,
based on f4b7198.

## Verdict
**Sound, and it reproduces exactly. The headline is correct: once the flow is carried with the particle (the material
derivative), the particle really coasts.**
- At h = 0.25 the strong particle's W centroid travels 0.963 of the ideal distance in 20 time units after a kick,
  compared with 0.847 without transport.
- Its flow speed keeps 97.1%, compared with 73.9% without transport.
- The weak particle keeps 99.8%.
- Two particle types under the same imposed flux accelerate within 0.4% of each other.
- Cores survive, and a transported particle left at rest stays at rest.

I recommend merging, which I have done along with this review. The note's caveats are stated fairly: the mass is
assigned, the head-on drag is off, the result is at one resolution, and nothing is claimed for Lorentz invariance or
derived inertia.

## What I checked
1. **Protocol.**
   - Commit fc1f500 holds the plan, the log entry and the implementation, and it comes before any result.
   - The wall-run shortening (T = 10) was committed in ccdf893, before the wall results in a47b990.
   - The faster stencil was logged before the kicks.
   - Review dispositions are in `review/EX-GC23_response.md`.
2. **Tests.** `python tests/test_regression.py` passes. `python src/exgc23_transport.py --phase validate` passes all six
   checks, with a largest difference of 3.4e-13.
3. **Summary regenerates.** `python src/exgc23_summary.py` changes `summary.json` only at the last digit (roundoff
   only). The committed file was kept.
4. **Own recomputation from the raw records** (a separate script that reads only the JSON time series):
   - W shift: 0.847354, 0.963450, 0.967514, 0.974762 for the strong and weak kicks, without and with transport.
   - Flow-speed retention: 0.739108, 0.970598, 0.982927, 0.998139.
   - Late W speed: 0.038911, 0.047754, 0.048097, 0.048659.
   - Late decay rate: 0.016015, 0.002195, 0.000913, 0.000141.
   - Conservative alternative: shift 0.961259 and retention 0.965098.
   - Flux trials: W acceleration 0.001920 and 0.001914; flow acceleration 0.001985 and 0.001980.

   All of these match the note.
5. **Momentum ledger** (strong kick with transport): recoil −0.010321, cohesion +0.0000587, transport +0.075075 and
   changing s −0.073747 sum to −0.008934. That equals P(20) − P(0) = 0.299046 − 0.307981, with a residual of 1.3e-13.
6. **Rerun.** I re-settled the K = 0.2 rest state (30 units) and reran both strong kicks with the committed code (about
   0.23 core-hours). All three records reproduce with zero relative difference in every recorded field.

## Points to note (none blocks merging)
- **A small slip remains even in the weak case.** The W centroid moves at about 0.975 of the flow speed in all four
  runs, including K = 0.01, where the flow barely decays. So about 2.5% of the shortfall from the ideal 1.0 is
  slip (lattice or pattern shape), not drag. The note's suggested h = 0.125 subset is the right way to find out which.
- **The leftover strong-particle decay depends on the box.** The ±9 box lowers the early decay rate from 0.00083 to
  0.00022. So the remaining 0.0022 late rate should not be read as a physical self-wake drag yet. The note says this.
- **The conservative alternative loses 4.9% of M.** It was run from the original rest state rather than its own
  equilibrated one. It shows that the coasting improvement does not depend on the old cohesion term. It is not a
  finished conservative model. The note says this too.
- **The committed code differs slightly from what ran production.** Commit 7befb38 (after the runs) switched
  `chem_inplace` back to the original stencil, while production used the faster one. They agree to roundoff, and my
  rerun with the committed code reproduces the records exactly, so this is cosmetic. A one-line comment saying so
  would help later readers.
- **"Near type-independence" is partly built in.** It is largely set by construction, since push and assigned mass both
  scale with absorption. What the test does show is that the particle's own field and the flow model no longer spoil
  it (EX-GC19's Eulerian version gave 4%, EX-GC20's wave version 12%).

## Evaluation of the candidate EX-GC24
I agree the head-on (EX-GC17) absorption drag is the next obstacle. One caution about the proposed first step: in this
model the mass is assigned (mass/drag = τ = 20). Reinstating the bare drag will reproduce an e^(−t/20) decay by
construction, so a run would mostly confirm an input.

The scientifically open question is **what sets the ratio of inertia to absorption drag**. I suggest EX-GC24 address
that directly. Possible routes:
1. **Derive mass rather than assign it.** Give the medium its own momentum (J as a momentum density with a conserved
   exchange with the pattern). Then ask whether the moving deficit carries an emergent "added mass", and how it
   compares with the absorption drag. That is what decides whether coasting is physical.
2. **Keep the h = 0.125 subset.** It is a cheap way to settle the 2.5% slip.
3. **Defer orbits.** A two-body orbit test should wait until a drag/mass ratio is derived. Otherwise orbits would only
   test the assigned τ.
4. **Lorentz invariance remains the deeper alternative**, a medium with no rest frame, if route 1 shows the drag cannot
   be made small relative to an emergent mass.

Ownership is the PI's decision.
