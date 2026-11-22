# EX-GC27A/B and EX-GC28A — Claude review (EXPLORATORY)

Reviewer: Claude, 2026-11-22. Reviewed ChatGPT's EX-GC27A/B (plan c402118, LOG e31257d, results through 86603aa) and
EX-GC28A (plan 5c4b5a3, LOG aae6dbb, regulator amendment d1f142d, results through 77f12fd). Plain Skyrme and φ⁴
reference models, not Model Gc. Closure is the PI's.

## Verdict

Both studies are sound and honestly reported. I found no mathematical or coding error, and every number I checked
reproduces. The conclusions as written stand. The issues below are a state mismatch in 27A that was disclosed but not
carried through to the binding number, one comparison that 28A leaves out and that I think is its most useful result,
and some wording. None of them needs a rerun.

## What I checked

1. **27A/B algebra.** The identity `k = (a+b)/2 + (a−b)v/2 + (c−2a+b)p` follows from E = E₂+E₄+E_V and the 3D
   virial v = (E₂−E₄+3E_V)/E. Every table entry in the note recombines from the committed baseline `parts`. For
   example, independent minus consistent B2−B1 = −4Δp = −4 × 0.000659 = −0.002636 (h = 0.2) and −4 × 0.000444 =
   −0.001776 (h = 0.15). The raw binding sensitivity .630 (h = 0.2) recomputes from the 26E1 energies.
2. **27B mass normalisation.** C = 2m² is right for this energy normalisation: (1−φ₀) ≈ π²/2 against a gradient term
   (∂π)² with no ½. The potential gradient sign (−C on φ₀) is correct, and the 1.77e-9 FD check covers it.
3. **Consistency with closed EX-GC26 results.** The (2,0) coupling is the Nordström weighting used in 26E. 27A's
   B2/B1 k ratio of 0.948418/0.935938 = 1.0133 at h = 0.2 reproduces the 26E1 lattice prediction of 1.013. At
   h = 0.15, 27A's 0.9638/0.9704 matches 26E2's 0.965/0.970. The B = 2 states match 26E1/26E2 (coarse is
   bit-identical).
4. **28A rerun.** I ran `src/exgc28a.py` and `src/exgc28a.py canonical` on numpy 2.5.3, scipy 1.18.1 and Python
   3.13. All 92 rows (12 kink, 40 radial, 40 canonical) reproduce raw, Born and subtracted values **bit-identically**,
   in about 8 CPU-s.
5. **28A independent recomputation** (`review/EX-GC28A_check/zchk.py`). I wrote my own hedgehog minimisation and built
   the canonical potential by numerical second variation of the energy density, not from the analytic formula in the
   code. R = 24:

   | case | h | mine | ChatGPT |
   |---|---:|---:|---:|
   | C = 0 | .1 / .05 / .025 | −.169259 / −.166822 / −.165989 | −.169261 / −.166827 / −.165998 |
   | C = .178396 | .1 / .05 / .025 | −.177390 / −.174589 / −.173622 | −.177393 / −.174593 / −.173631 |

   They agree to 1e-5, which confirms the effective potential (including the (√K)″/√K term and the matched l = 1
   vacuum 2/r² + C/2A) and the first-Born subtraction. My classical radial energy of 1.23155 is the standard B = 1
   Skyrmion value.
6. **Kink benchmark.** −0.666255 = μ(1/(2√3) − 3/π) is the DHN mass correction for meson mass 2μ. First-Born
   subtraction with the same lattice in defect and vacuum is the no-tadpole scheme. The 0.041% agreement is a genuine
   external check of the method.
7. **Protocol.** For both studies the plan and LOG row were committed before the results. Budgets were 0.40 and
   0.002 core-h, well under the 2 core-h caps. Both are labelled EXPLORATORY. The white paper and README are untouched.
   27 respected "no ZPF/dynamics". Commit times are consistent with the stated CPU (27: 1444 CPU-s between 06:46 and
   07:11; 28A: 8 CPU-s between 07:40 and 07:43). I can't independently verify the ordering, because all pushes come
   from ChatGPT's environment.

## Issues

**1. 27A fine-grid B = 1 state is not the 26E2 state, and the binding number depends on it (moderate, disclosed in
part).** The fine B1 array was reused from the 26B2 rest state. Its E₂ is 0.1% below 26E2's (0.586539 vs 0.587124),
with the same total to 1e-6, and its gradient (7.2e-6) is ~200× larger than the B2 state's (4.0e-8). The plan's
equivalence test (total energy within 0.1%) is too loose for subtractive quantities. Recomputed with the 26E2 B1
energies:

| quantity (h = 0.15) | 27A as reported | with 26E2 B1 state |
|---|---:|---:|
| k(2,0), B1 | .963826 | .964790 |
| raw binding sensitivity (2,0) | .812 | .835 |
| virial difference Δv (B2−B1) | .00593 | .00559 |

The note mentions the 0.1% E₂ difference but does not carry it into the binding sensitivity. That number shifts by 3%
on the state alone. The 27B virial-separated result (−4Δp) uses freshly relaxed baselines and is barely affected,
although the fine B1 B-baseline also stopped at its iteration limit with a gradient of 1.2e-6, the least converged of
the four. Suggested fix: one sentence in the note, or quote both values.

**2. The "scale-consistent" coupling is metric coupling (framing, not an error).** (a,b,c) = (2,0,4) is exactly
minimal coupling of the Skyrme action to the conformal (Nordström) metric g = e^{2η}η. √−g·g^{μν}g^{ρσ}… gives
e^{2η} on the two-derivative term, e⁰ on the quartic Skyrme term (conformally invariant in 3+1D) and e^{4η} on the
potential. So 27's main result is the known statement that metric coupling gives a universal static source (von Laue,
as already used in EX-GC25) and non-metric coupling (2,0,0) does not. The note's "common scaling principle across all
independent contributions" is the equivalence principle in its metric form. I'd say so explicitly, because it names
what a microscopic vacuum theory would have to produce. It also makes 27 a calibration of the lattice against a known
theorem, which is what the outline predicted.

**3. "Small" composition difference (wording).** The virial-separated species difference is −4Δp ≈ −7e-4 at
h = 0.15. A two-point h² extrapolation of Δp (.000659 → .000444) leaves ≈ 1.7e-4, so it does not obviously vanish.
The note correctly calls the continuum value unresolved. If real, though, it would be about ten orders of magnitude
above Eötvös limits. So "small" here means small relative to lattice error, not physically tolerable. That supports
the note's recommendation: the principle matters, and the decimal places don't.

**4. 28A leaves out the universality comparison that matters (substantive, from committed records).** For gravity, the
question is whether the zero-point part of a particle's mass sources in the same proportion as its classical part. The
note compares consistent with independent coupling instead. From `radial_summary.json` and `canonical_summary.json`
(R = 24):

| C, h | k classical, consistent / independent | k_Z = (dZ/dη)/Z, consistent / independent |
|---|---:|---:|
| .178396, .1 | .9998 / .9173 | 1.0365 / .8698 |
| .178396, .05 | 1.0000 / .9176 | 1.0132 / .8510 |
| 0, .05 | .9997 | 1.0116 |

Under non-metric (independent) coupling, the ZPF part sources about 7% less per unit energy than the classical part of
the same particle (≈ 0.84 after removing the ~1% fixed-grid bias, vs 0.918). That is a composition dependence inside a
single particle. Under metric coupling with a comoving cutoff, both are 1 by construction. The fixed-grid consistent
values (1.033 at h = 0.1, 1.012 at h = 0.05) show that a fixed, non-comoving cutoff itself breaks universality of the
ZPF term at O(h²). That error vanishes here only because the Born-subtracted Z is finite.

**5. The subtracted Born term is the open part (scope).** The divergent first-Born (tadpole) piece is subtracted and its
η-dependence is not reported. In a Sakharov-type argument, whether that cutoff-dependent piece, or its counterterm,
couples to η is the actual question. Its η-derivative is free to compute from the existing code. I'd record it, with
its cutoff dependence, before 28B.

**6. 28B expectation should be registered (advice).** The proposed connected-vacuum-energy test is a good control, but
the outcome for weak localised potential defects in a massless 3D scalar is known at leading order: the connected
one-loop energy falls off as 1/d³, not the 1/d needed for a Newtonian potential. In 1D it is different. The 28B plan
should state the predicted exponents in 1D and 3D before running, so that a Casimir-type result isn't read as a
gravity candidate.

**7. Minor.**
- The 28A table header "C0" means C = 0, while C0 elsewhere (27B, and the code's variable) means 0.178396.
- "Baseline radial energy shifts ~4.6% between cases" refers to Z (−.1660 → −.1736). The classical radial energy
  shifts 2.5% (1.2315 → 1.2625).
- In 27, `optimization.E2` in the B-baseline files is E₂ + E_V (`parts` is correct).
- The 28A `raw_manifest.md` hashes can't be verified by regeneration, because the raw JSON embeds CPU timestamps. The
  numerical content does regenerate exactly (item 4). Hashing a timestamp-free copy would make the manifest checkable.
- Raw 28A spectra exist only in ChatGPT's scratch. That's acceptable given the 8 s deterministic regeneration.

## Bottom line for the PI

- 27A/B correctly shows that a universal static source comes from the coupling being metric, not from topology. This is
  a calibration of a known theorem, and it is consistent with 26E1/26E2.
- 28A's numerics are right: kink benchmark, rerun, and my independent Z all agree.
- Its one-channel result already says something the note doesn't: with a non-metric coupling, the ZPF part of the mass
  would not gravitate in proportion to the classical part.
- Nothing here suggests ZPF as a source of 1/r gravity. Before 28B, I'd record the Born-term η-derivative and register
  the expected Casimir-type distance exponents.
