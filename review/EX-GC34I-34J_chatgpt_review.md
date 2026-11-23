# ChatGPT review of EX-GC34I and EX-GC34J

2026-11-23 | Independent review of the exploratory results

**Verdict.** GC34I and GC34J advanced the repair of the vacuum-wave instability. GC34I constructed a kinetic action with exact linear re-slicing symmetry. GC34J found two candidates with two near-light-speed waves and no growing mode in the sampled spectrum. Neither study produced a candidate that passed all registered criteria. The remaining deficit for J1 and J2 is the first-order error in the moving-source benchmark, alongside the registered time-domain failure and its encouraging post hoc follow-up.

## GC34I

GC34I imposed exact re-slicing as a null identity of the linear lattice equations. Short-range solutions exist for the forward-difference gradient. The tested centered gradients gave no solution at ranges 1 and 2. The seventh, non-metric edge mode remained at the lattice cutoff in the registered candidate, leaving two long-wavelength branches.

The registered If1 candidate failed G2 and G5. Its waves were about 1.1% slow and its moving-source error was 2.3% against a 1% limit. Fitting the normalization to GC34E was a poor continuum anchor because that kernel violates the identity even at leading order. Choosing the continuum normalization, μ = −1, recovered the wave-speed gate post hoc. The moving-source error remained. Registering that normalization in GC34J was a sound follow-up.

## GC34J

- **Vacuum spectrum.** J1 and J2 gave two low-frequency branches with phase speeds within 0.076% of the matter-wave reference in the tested long-wavelength range. Neither had a growing root in the sampled cases. J3 grew at short wavelengths. Exact re-slicing alone therefore does not guarantee stability for every kinetic action.
- **Time evolution.** All three failed the registered G4′ test, whose initialization seeded constraint-violating zero-cluster directions. The post hoc rerun restricted initial data to the constraint-compatible subspace. Residuals stayed below 8 × 10⁻¹², and the allowed zero-frequency seeds changed the tidal signal by less than 7 × 10⁻¹³ relative to the wave. This supports harmless residual coordinate freedom for the tested runs. The two directions at one spacing do not settle J3’s short-wavelength instability.
- **Moving sources.** All three failed both registered E4 gates while preserving the GC34D static kernel. The committed worst-case errors nearly double between kh = 0.025 and 0.05, reaching about 1.25%. I checked the ratios from the committed table. They give first-order convergence over the tested range.

**Claim to narrow.** The results strongly suggest a structural origin for the first-order error. Three substantially different kernels give almost the same error, and neither searched range could match the earlier kernel’s first-order coefficients. The E4 table reports a maximum over nine cases. It does not identify which curvature component sets that maximum or establish that the leading observable error is invariant under all permitted changes to the kernel. Failure to match the full set of coefficients still allows the possibility of matching the combination that matters for the observable. The staggered-shift explanation is plausible and remains a hypothesis.

## Priorities

1. Keep J1 as the linear vacuum-wave benchmark. Preserve the registered verdicts and label the constraint-filtered evolution as post hoc.
2. Before searching wider stencils, identify the moving-source case and curvature component that set the leading error. Test whether the coefficient varies within the space of kernels that obey the identity, with source coupling and curvature reconstruction consistent with the lattice transformations.
3. If the leading error is fixed in that space, test a staggered shift. Initialize the registered time evolution with constraint-compatible seeds from the outset.

The broader implication for Model Gc is that enforcing a symmetry succeeded where averaging strain rates stalled. The symmetry was a design condition in these studies. The emergence of shared geometry and universal matter coupling from the original Model Gc rules remains open.

## Review scope and sources

I reviewed the plans, source, and available results. GC34J was on Claude’s branch at the time of review. I checked the convergence ratios from the committed data. I did not rerun the spectral or evolution calculations.

- GC34I [plan](../spec/EX-GC34I_plan.md), [results](../results/EX-GC34I/EX-GC34I_note.md), and [source](../src/exgc34i.py)
- GC34J [plan](https://github.com/RichardBenear/Model-Gc/blob/claude/project-thread-2twuh8/spec/EX-GC34J_plan.md), [results](https://github.com/RichardBenear/Model-Gc/blob/claude/project-thread-2twuh8/results/EX-GC34J/EX-GC34J_note.md), [E4 data](https://github.com/RichardBenear/Model-Gc/blob/claude/project-thread-2twuh8/results/EX-GC34J/e4.json), and [post hoc evolution](https://github.com/RichardBenear/Model-Gc/blob/claude/project-thread-2twuh8/results/EX-GC34J/sync_posthoc.json)
