# Independent review of EX-GC45

**GC45 demonstrates a low-spatial-charge limit for a single carrier, with a mobility trade-off. It also finds substantial binding-dependent spatial charge in the tested contact model.** The binding result survives at several representative points, but the claim that every point exceeds the registered threshold does not survive the four-band check. [Plan](../spec/EX-GC45_plan.md) · [Results](../results/EX-GC45/EX-GC45_note.md)

## GC45A — The previous exciton remains highly sensitive to strain

Using the best momentum from the wedge scan, the total spatial charge is 1.40–2.17 times the exciton energy for the unstrained contact law. The strain exponent needed to cancel the whole body's charge is positive and varies substantially with binding, confirming the distinction raised in the GC44 review. The zone-corner result at λ₀ = 0.1 has a sizeable energy separation from K = 0. At λ₀ = 0.2 the separation is as small as 0.0007. Evaluating the selected momenta on a finer lattice checks their energies but does not refine the momentum search itself. The charge conclusion is stronger than the claim that each exact minimum has been located. A finite crystal momentum at a dispersion minimum also does not imply a moving composite. Group velocity depends on the dispersion gradient. [Source](../src/exgc45.py) · [Results](../results/EX-GC45/EX-GC45_note.md)

## GC45B — Nearly dust-like at rest, with mobility that falls at the same rate

At g = 0.025 the three single carriers have Q_a/ΔE between 0.0075 and 0.0087. Both that ratio and the dispersion curvature decrease approximately as g² in the small-g comparison. This is a useful controlled limit: the on-site term stores rest energy while hybridization supplies mobility and a strain response. The response stays 28–74% of the bandwidth over the tested range, and the smallest-g bound composites were not measured. The strongest claim is therefore about a single carrier at its band minimum. The study has yet to show how small the full spatial charge can be at a specified nonzero composite mobility. [Results](../results/EX-GC45/EX-GC45_note.md)

## GC45C — Binding stress is real in this model, and its size depends on the comparison

The reported η is the binding-induced change in spatial charge divided by the positive binding energy. At fixed interaction parameters, this is also minus the logarithmic derivative of the binding energy with respect to the hopping scale. It measures sensitivity to strain under the chosen contact law. The experiment does not calculate passive free fall or a gravitational self-energy contribution, so η is not an observational Nordtvedt parameter. In the GC44 source equations it changes the active source through Q_a.

The four-band check gives a more precise reading of the threshold failure. Using each four-band energy to recompute the binding denominator against the same free threshold, I get η ≈ **−0.039** at the registered middle point. The strongly hybridized cases still give approximately **−2.10 to −1.01**, and the weakly hybridized cases give **+1.10 to +4.83**. Thus the broad spread and the opposite signs survive, while the all-points threshold does not. These are comparisons within the saved small-box four-upper-band model. It still freezes the valence sea. Also, m₀ and g change together across the three backgrounds, so the sign reversal cannot be attributed to g alone. [Projected outputs](../results/EX-GC45/C.json) · [Four-band outputs](../results/EX-GC45/C5_extended.json)

The normalization also matters for the gravity implications. In the weakest-hybridization pair tested, the binding contribution is only **0.55–0.74% of the total pair energy**, despite η exceeding one. Its **total** spatial charge is still **3.8–4.0% of the pair energy**. Report both measures when comparing with the approximately 0.8% single-carrier result. [Pair outputs](../results/EX-GC45/C.json)

## Next tests

1. **Promote the four-band calculation to the main solver.** The existing contact matrix is diagonal plus a rank-one term, even with four upper bands per carrier. Extend the scalar bound-state solver to that larger basis and validate it against the saved dense checks. At larger volumes, compare at fixed U for projection error and at matched E_b/W for comparisons between different carrier models. Score both η and total Q_a/ΔE.
2. **Test the small-g bound-pair limit.** Hold m₀ fixed, include g = 0.025, and tune the binding to the same E_b/W across runs. Track the composite curvature, full spatial charge, and binding contribution separately. Repeat at another m₀ to separate the effects of hybridization and the on-site scale.
3. **Compare one independently specified strain law with the contact control.** Use a finite-range interaction whose distance dependence is fixed before the run. Measure the total derivative of the bound energy and the binding-only derivative across several binding strengths. A useful result would be suppression that persists across states at finite mobility. A separate propagation test would still be needed to connect the resulting sourced fields to actual deflection and delay.

I reviewed the plan, implementation, results note, and saved projected and four-band outputs. The renormalized binding ratios are a reanalysis of those outputs. The simulations were not rerun.
