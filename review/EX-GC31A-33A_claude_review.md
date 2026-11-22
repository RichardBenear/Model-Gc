# Claude review: EX-GC31A to EX-GC33A, and the GC33B proposal

2026-11-22 | Claude (reviewer) | Review of ChatGPT's exploratory series on `main` up to a7d3657. No closures recorded; closures are the PI's.

## Bottom line

1. **The work is technically sound.** Every exact claim I checked is correct, and the two scripts I re-ran reproduce
   the published numbers. The negative results are honestly labelled and nothing is over-claimed.
2. **The direction has drifted away from the programme.** Eleven studies in one day went from "does a causally
   invariant rule force one light cone?" to causal sets, quantum graph Hamiltonians, and finally the IKKT matrix
   model from string theory. Each study answered the previous one's obstruction by changing the ontology, and none
   used Model G or Model Gc dynamics. The series mostly re-derived known obstacles in the emergent-gravity
   literature, which it says openly (Gorard, Konopka et al., Hamma et al., Steinacker).
3. **GC33B: I recommend not running it as proposed.** Its main pass condition, a shared light cone for scalar,
   spinor and gauge fields on a matrix brane, is already a known theorem of that model, so the test can't fail in an
   informative way. The hard parts (does energy really source the metric, is there a proper spin-2 mode) are open
   research problems in that literature and won't be settled by a bounded screen. And a pass would be a result
   about IKKT, not about Model G or Gc.
4. **Keep one lesson from it.** In matrix models the cone is shared because *one operator* builds every field's
   kinetic term. That principle can be tested inside our own framework (suggestion below).

## What each step established (my reading)

| Study | Result | Value for the goal |
|---|---|---|
| 31A | Exact counterexamples: a schedule-invariant rule can carry two different signal speeds (fronts 1 and ½; or speeds 0.866 and 0.5 with a common microscopic front). | High. Causal invariance is necessary for the programme but does not give a common metric. Worth keeping. |
| 31B1 | 80 CRR Stage W rules, 0/14,320 two-event diamonds; driven by state mismatch, not shown DAG mismatch. | Low. These are CRR rules, not Model G rules. Inconclusive by its own account. |
| 31B1B–F | Five hand-built rules, each causally invariant *by construction* (no conflicting matches allowed), each with a proven structural defect: hubs, causal forest, 3-port bottleneck, ancestry-locked boundary, saturating rewiring. | Low to medium. Correctly stopped. Avoiding conflicts trivialises causal invariance, which is the hard part. |
| 32A | Same quantum ground state, same vacuum energy, different speeds (1±b). Fixed quadratic vacuum can't respond to energy. | Medium. Clean, correct; known in spirit. |
| 32B | Response tomography recovers graphs; first moment blind to kinetic energy at fixed density. | Medium. Good method; it's a benchmark on supplied geometry. |
| 32C | Positive zero-point energy favours deleting links, cliques and disconnected blocks, not lattices. | Medium. Correct and decisive for that idea. |
| 33A | General theorem (positive links → energy can't select geometry); fermion/boson comparator; universal-source condition = lapse identity S = κ(H−E0). | High for the last item: it states exactly what any "energy sources geometry" mechanism must satisfy. |

## Checks I did

* GC33A discriminant B: recomputed half-filled band energies, C12 −7.464102t vs 2C6 −8t, C8 −4.828427t vs 2C4 −4t,
  and the boson zero-point table at μ = 0.1, 1, 5 (7.67056/7.58234, 10.06566/10.06450, 31.16621/31.16621). All match.
* GC32A: verified the b-independent Gaussian ground state (B = KAK holds since K = I⊗√L commutes with M⊗I).
* GC33A discriminant A (concavity and monotonicity of ground energy in positive link weights) is a standard and correct argument.
* Re-ran `src/exgc32b_response.py` (0.9 s): m = 6 bond response difference 0.297446 vs published 0.29745.
* Re-ran `src/exgc32c_relational.py`: N40 degree-4 necklace U0/N = 0.934996 at μ = 0.1, as published.
  (Same author code; this is a rerun, not independent recomputation.)

## Concerns

1. **Pending-review backlog.** 31A, 31B1, B1B–F, 32A–C and 33A (11 studies) are all "review/closure pending"; several
   say "review deferred". Each recommended next step was taken before the previous one was reviewed.
2. **Straw-man negatives.** Several studies refute sufficiency claims nobody holds (schedule invariance alone ⇒
   common cone; vacuum state alone ⇒ speeds). Correct, but low information per study.
3. **Lost connection to Model G/Gc.** The one result in the Gc line that actually worked (EX-GC25/26: one lossless
   knot field, common-action coupling, universal fall, 1/d², binding gravitates) isn't used. Its open problems were
   that the common coupling was *imposed* and that scalar gravity fails light bending. Those are the right targets.
4. **The parent question is still open.** MG189 left causal completion of the Model G hypergraph rules unresolved.
   31B1 tested CRR rules instead.

## Recommendation on GC33B

**Do not run as proposed.** If you want something from it, do a modified version that stays in our framework:

**Suggested alternative (EX-GC33B′, paper first, then a tiny lattice check):** a "single hopping operator" rule.
Make every token species in a Gc-style lattice/hypergraph move only through one shared link variable (one
operator, like ad X in the matrix model), so a common cone holds by structure rather than by tuning. Then apply
33A's lapse test: when a particle's energy changes the link variable, does the response equal κ(H−E0) for all
states, including binding and kinetic energy? Two clear outcomes:
* pass → a concrete, on-programme mechanism for "energy sources a common metric", which can then be made
  causally invariant and tested for tensor (light-bending) behaviour;
* fail → a sharp, Model-Gc-specific negative.

Either way, review and close or set aside the 31A–33A backlog first.
