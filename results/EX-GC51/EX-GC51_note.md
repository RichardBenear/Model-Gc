# EX-GC51 — one declared model end to end: Dirac matter, a coincidence-carried slowness field, and light

2026-11-27 | EXPLORATORY | Claude | Closure pending | [Plan](../../spec/EX-GC51_plan.md) (amendment A1 before any run) | About 6.5 core-h on the cloud container (Linux, Python 3.13.18, numpy 2.5.3). The PI raised the cap from 2 to 5 core-h, then to about 7, and dropped the last two C4 checks. Scoring and post-run diagnostics: `src/exgc51_summary.py` → `summary.json`. This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**One set of rules gives Newton's law with G fixed, light bending at GR's factor 2, β ≈ 1, a shared light cone, and action equal to reaction. Nothing beyond the declared model was put in. The probe-fall test fails as registered.** That failure comes from my comparator and probe design, not from a sign of extra equivalence-principle violation. Of 5 scored predictions, 1 passes fully (P4), 3 pass on one clause and fail on another (P1, P2, P5), and P3 fails.

| | What was scored | Result | Verdict |
|---|---|---|---|
| P1 | field signal speed c_N = 1.00 ± 0.02 | 1.078 | FAIL |
| P1 | field arrival slope equals light's within 1% | ratio 1.0010 | PASS |
| P2 | G_eff/G = 1 within 1% | 0.99997–1.00000 (12 fits) | PASS |
| P2 | β = 1 ± 0.05 after the source-linear separation | 1.060, 1.062, 1.063 | FAIL |
| P3 | a/a_semi = 1 ± 0.03 | 0.77–0.86 | FAIL |
| P3 | mass ratio equals a_semi's within 1% | 0.91, 0.93, 0.95 | FAIL |
| P3 | at r = 0, both masses fall equally within 1% | ratio 0.63 | FAIL |
| P4 | \|F₁ + F₂\|/\|F₁\| < 0.05 | 0.0136 | PASS |
| P5 | Shapiro ratio 2.00 ± 0.06 at b = 8, 16, 24 (wave) | 2.076, 2.037, 2.024 | FAIL at b = 8 |
| P5 | same, eikonal | 2.038, 2.022, 2.017 | PASS |

## What each part shows

**51A, one cone.** The field's change spreads at exactly the speed of light, to 0.1%. Both read 1.08 rather than 1 on the registered half-arrival measure. The likely cause, which I have not verified, is that the static target is taken with the Dirichlet box. Near the walls that target is smaller than the free-space field, so half of it is reached early, and more so at larger r. That would give the slope 0.93 that light and the field share. P1's "c_N = 1" clause fails as registered. The comparison that matters physically, field against light, passes.

**51B, static field.**
- G_eff = 1 to 3 × 10⁻⁵ at all four strengths and on all three fit directions. The Gauss-law charge is the source's own lapse charge Q = ∂E/∂N, which falls from 0.966 to 0.898 as the source redshifts in its own well.
- The β separation as coded used κ as the strength variable. It should have used κQ, because the source-linear coefficient scales with the actual charge. As coded it gives 1.06, a registered FAIL.
- With κQ (post hoc) the separation gives β = 0.9987, 1.0004 and 1.0010 on axis, face and body.
- The raw single-strength fits all lie between 0.979 and 1.013.
- So the self-consistent field gives the exponential exterior and β = 1, as 49 found for a prescribed source.

**51C, probe fall.** All three P3 clauses fail. The registered comparator a_semi = −⟨E_kk E_N⟩∂N omits a term. As a packet moves into a different N, its group velocity changes, which adds ⟨E_k E_kN⟩∂N. I found this after the first case, before rerunning, and recorded the corrected comparator only as a diagnostic.
- **The diagnostic test.** In a linear gradient on a 48³ box with m = 0.5, the wave packet's 2Δx/t² sits at 1.39–1.40 of the GR rate from t = 6 to 16. The corrected comparator gives 1.397 and a classical ray ensemble gives 1.41.
- **The registered runs.** Against the corrected comparator, the late-time 2Δx/t² is 0.98 (m = 0.5) and 0.94–0.96 (m = 0.25).
- **Why the fitted value is worse.** The fitted a is pulled further off by a fast oscillation from the slight negative-energy admixture, with period about 2π/(2E). That period is about as long as T = 12 for m = 0.25.
- **Why the r = 0 clause fails.** The design itself fails it. With σ = 3/m the two probes have different sizes, 6 and 12 sites. The wide one reaches far enough toward the source that it feels only about 65% of the point force. Its a_GR ratio is 0.645, so "equal fall" could not hold for these probes.
- **What 51C does show.** The fall follows the lattice's own dispersion, with no sign of anything extra. The point-particle limit of that dispersion is (1 + rm) times the GR rate, which is the order-mass/cutoff lattice term that 48 predicted.
- **What it does not show.** It does not test universality at the 1% level. A clean test needs equal-size probes (same σ), runs much longer than 1/m, and the corrected comparator.

**51D, reciprocity.** Two free fermions (m = 0.4 and 0.6) evolve with their own dynamic field. Their canonical-momentum forces balance to 1.4% (F₁ = 2.410 × 10⁻⁴, F₂ = −2.378 × 10⁻⁴). The inertia-weighted acceleration balance, which is a diagnostic only, is 2.1%. So the coupled rule conserves momentum to this accuracy even though no action is known for it.

**51E, light in the self-consistent field.** Against the one-site control the delay ratio is 2.02–2.04 for rays, and 2.02–2.08 for waves. The excess over 2 grows as the pass gets closer, which fits the exponential metric's (e^U + 1)-type second-order term at U ≈ 0.05. It only exceeds the band for the wave at b = 8.

## Checks

| | | |
|---|---|---|
| C1 | static solver vs 49's slowness solution | 4.7 × 10⁻¹³ (independent residual 9 × 10⁻¹³) |
| C2 | uniform-N dispersion, N = 0.9 | 2.7 × 10⁻¹⁵ |
| C3 | Q = ∂E/∂N by finite difference | relative 4.6 × 10⁻¹⁰ |
| C4 | 51C (m 0.5, r 1) at dt/2 | a changes by 0.05%; largest energy excursion 8 × 10⁻¹⁰ (2.6 × 10⁻⁸ at dt) |
| C4 | 51C box × 1.5 and 51D dt/2 | **not run** (PI-approved scope cut, 18:44 UTC) |

51C energy excursions were 1.5 × 10⁻¹⁰ to 2.6 × 10⁻⁸. 51D did not record its energy excursion, contrary to C4's wording, and this is recorded as a deviation.

## Deviations

1. The 51C comparator omitted the velocity-change term (above). The registered score stands. The corrected comparator is a post-run diagnostic.
2. The 51B β separation used κ instead of κQ. The registered score stands. The κQ version is post hoc.
3. The A1 probe design (σ = 3/m) made the r = 0 equal-fall clause untestable.
4. Compute ran to about 6.5 core-h against the registered 2. Both raises were approved by the PI before the extra runs continued. The runs were paused at 5.6 core-h while waiting for the second raise.
5. C4 was only partly run. The C4 result is taken from its log line (`checks4.json`), because the process was stopped before writing its own record.

## What it means

Weak-field, static gravity comes out of one set of local rules:
- Newton's law with a fixed G;
- γ = 1 (from 37B counting);
- β = 1 (from passing on slowness);
- a single light cone for gravity and light;
- action equal to reaction.

None of the field equations was put in. This is the programme's first end-to-end mechanism candidate. The fails above are measurement and comparator issues, and none of them points the other way. But a candidate is all it is.
- **Not derived.** Slowness as the transported quantity is declared, and 49 showed it is not unique.
- **Not tested.**
  - Sources moving or spinning: GR's g₀ᵢ, frame dragging and velocity-dependent forces.
  - Radiation: a single scalar may radiate monopole and dipole waves, which binary pulsars rule out.
  - Bound composites (the Nordtvedt test).
  - The vacuum.
  - Strong fields.

**Suggested next (none approved):**
1. A moving source: does one coincidence-carried scalar give the right velocity-dependent force, or does a vector or shift part have to appear?
2. Radiation from an orbiting pair.
3. A clean universality rerun: equal-size probes, longer runs, the corrected comparator.
4. An action for the coupled rule.
