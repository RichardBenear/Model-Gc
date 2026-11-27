# EX-GC46 — heavy pairs with all four upper bands, and the small-g limit

2026-11-27 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC46_plan.md), pushed before any run. PI request 2026-11-27 04:06 UTC, go-ahead 04:09 UTC ("go ahead and close 43, 44, 45 and write the plan for 46 and then execute the runs"). Follows the [Codex review of GC45](../../review/EX-GC45_codex_review.md), tests 1 and 2. PI CLOSED 2026-11-27. Cloud container (1 core), Python 3.13, numpy 2.5. About 0.02 core-h. This is our construction. It is not LaViolette's Model G, and it is not the Chain Theory construction.

**With all four upper bands, the binding stress survives and its sign flips with the hybridisation g, not with the on-site scale. As g → 0 it approaches a finite, positive ratio.** All of P2 passes; P1 passes except its first clause, a 5% projection-error limit, which fails at the strongly hybridised point.
- **The projection was a good approximation except at strong hybridisation.** The four-band binding energy differs from the projected one by 0.3% at (0.7, 0.05), 2–3% at (0.5, 0.1) and 13–18% at (0.3, 0.2) (fixed U).
- **η at the middle point is −0.033 (matched) or −0.038 (fixed U).** That reproduces Codex's −0.039. The values for strong hybridisation are −2.34 to −1.06 and for weak +1.11 to +4.89.
- **The small-g limit is dust-like and mobile at the same rate.** At g = 0.0125 the pair's total spatial charge is 0.2–0.3% of its energy, its binding energy is 0.02–0.06% of its energy, and its mobility has fallen as g². The ratio η tends to a finite positive value (0.4–2.7).

## 46A — four-band solver at 45's nine points (L = 32)

η = δQ_a/E_b, with δQ_a the pair's t∂E/∂t minus its free-pair threshold's. "Fixed U" is 45's U; "matched" re-tunes U so that E_b is the same fraction of W_h in the four-band model.

| (m₀, g) | E_b/W_h | η, 45 projected | η, four-band, fixed U | η, four-band, matched | E_b change at fixed U | δQ_a change at fixed U | total Q_a/ΔE (matched) |
|---|---|---:|---:|---:|---:|---:|---:|
| (0.3, 0.2) | 0.05 | −2.68 | −2.17 | **−2.34** | +18% | −4% | 0.771 |
| | 0.15 | −1.69 | −1.35 | **−1.43** | +14% | −9% | 0.729 |
| | 0.30 | −1.31 | −1.01 | **−1.06** | +13% | −13% | 0.744 |
| (0.5, 0.1) | 0.05 | +0.31 | +0.38 | **+0.40** | +2.6% | +27% | 0.118 |
| | 0.15 | −0.104 | −0.038 | **−0.033** | +2.1% | −63% | 0.117 |
| | 0.30 | −0.233 | −0.164 | **−0.162** | +2.1% | −28% | 0.113 |
| (0.7, 0.05) | 0.05 | +4.87 | +4.88 | **+4.89** | +0.3% | +0.5% | 0.038 |
| | 0.15 | +2.06 | +2.07 | **+2.08** | +0.3% | +0.9% | 0.039 |
| | 0.30 | +1.09 | +1.10 | **+1.11** | +0.3% | +1.8% | 0.040 |

**Scoring P1.**
- **E_b changes by less than 5% at fixed U: FAIL** at (0.3, 0.2) (13–18%); passes elsewhere.
- **δQ_a changes by more than 5% at (0.5, 0.1): PASS** (+27%, −63%, −28%).
- **Signs: PASS.** η is negative at (0.3, 0.2) and positive at (0.7, 0.05) for all three fractions.
- **|η| < 0.1 at (0.5, 0.1), 15%: PASS** (0.033 matched, 0.038 fixed U).
- **Spread of η across the nine points above 20%: PASS.**

## 46B — the small-g limit (L = 32, four bands, U tuned for E_b/W_h)

| m₀ | g | E_b/W_h | total Q_a/ΔE | binding part δQ_a/ΔE | E_b/ΔE | **η** | pair mobility t_pair |
|---|---|---|---:|---:|---:|---:|---:|
| 0.5 | 0.0125 | 0.15 | 0.00209 | 2.2e-4 | 2.8e-4 | **+0.763** | 1.47e-5 |
| | | 0.30 | 0.00211 | 2.3e-4 | 5.6e-4 | **+0.413** | 1.27e-5 |
| 0.5 | 0.025 | 0.15 | 0.00829 | 8.1e-4 | 1.1e-3 | **+0.713** | 5.88e-5 |
| | | 0.30 | 0.00834 | 8.5e-4 | 2.3e-3 | **+0.376** | 5.07e-5 |
| 0.5 | 0.05 | 0.15 | 0.0322 | 2.4e-3 | 4.6e-3 | **+0.529** | 2.34e-4 |
| | | 0.30 | 0.0321 | 2.2e-3 | 9.2e-3 | **+0.240** | 2.02e-4 |
| 0.5 | 0.1 | 0.15 | 0.1166 | −6.5e-4 | 0.0196 | **−0.033** | 9.27e-4 |
| | | 0.30 | 0.1132 | −6.5e-3 | 0.0400 | **−0.162** | 8.00e-4 |
| 0.7 | 0.0125 | 0.15 | 0.00266 | 5.7e-4 | 2.1e-4 | **+2.73** | 1.14e-5 |
| | | 0.30 | 0.00275 | 6.6e-4 | 4.2e-4 | **+1.57** | 1.01e-5 |
| 0.7 | 0.025 | 0.15 | 0.01046 | 2.2e-3 | 8.4e-4 | **+2.59** | 4.58e-5 |
| | | 0.30 | 0.01076 | 2.5e-3 | 1.7e-3 | **+1.47** | 4.04e-5 |
| 0.7 | 0.05 | 0.15 | 0.0393 | 7.1e-3 | 3.4e-3 | **+2.08** | 1.84e-4 |
| | | 0.30 | 0.0399 | 7.5e-3 | 6.8e-3 | **+1.11** | 1.62e-4 |
| 0.7 | 0.1 | 0.15 | 0.1294 | 1.1e-2 | 0.0141 | **+0.783** | 7.53e-4 |
| | | 0.30 | 0.1265 | 6.5e-3 | 0.0286 | **+0.228** | 6.62e-4 |

**Scoring P2: PASS.**
- Total pair Q_a/ΔE at g = 0.0125 is 0.0021 and 0.0027 (below 5e-3).
- The binding share E_b/ΔE is 2–6e-4 (below 1e-3).
- t_pair rises by 4.00 (m₀ = 0.5) and 4.01 (0.7) from g = 0.0125 to 0.025: g² within 1%.
- η has a finite limit: η(0.0125)/η(0.025) − 1 is +7%, +10% (m₀ = 0.5) and +6%, +7% (0.7), all below 20%. It is still creeping up slowly at g = 0.0125.

**P3 (descriptive).** At small g, η is positive at both m₀, and about 3.5× larger at m₀ = 0.7 than at 0.5. At m₀ = 0.5 it changes sign between g = 0.05 and 0.1. At m₀ = 0.7 it stays positive up to g = 0.1. So the **sign flip belongs to g**, and the on-site scale m₀ sets the size. The negative values at (0.3, 0.2) in 46A are in the same strongly hybridised regime.

## What it means (inference)

- **The 45 result is robust to the band restriction.** Four bands change η at the percent level where hybridisation is weak, and by about 15% where it is strong. The binding stress is not a projection artifact.
- **In the dust limit, binding stress is a fixed fraction of the binding energy, not a vanishing one.** As g → 0, E_b and δQ_a both fall as g², so η approaches a constant of order 1 (0.4–2.7, set by m₀ and the binding fraction). Dust-likeness (total Q_a/ΔE → 0) and an order-one stress ratio therefore coexist. The stress is tiny in absolute terms: binding contributes a few parts in 10⁴ of the pair energy.
- **Where this leaves the equivalence-principle question.** The violation is order η·E_b/ΔE. For the weakest pair here that is below 10⁻³ of the energy. Whether that is acceptable depends on the nuclear binding fractions one wants to mimic (≈ 1% for nuclei), so it is not negligible in principle. The contact law has no free exponent: a finite-range interaction with a pre-fixed strain law (Codex test 3) is the next way to ask whether η can be made to vanish.

## Checks

- **C1: PASS.** Q_N = ΔE for the four-band pair to 5.6e-12 relative.
- **C2: PASS.** With the basis restricted to band 3, the four-band solver equals 45's projected solver exactly (difference 0.0). Hellmann–Feynman U ∂E/∂U = ⟨V⟩ to 2.7e-10 relative.
- **C3: PASS.** A dense four-band pair matrix at L = 6 equals the rank-1 solver to 1.5e-14. (The driver's dense routine is the 45 one, so this check shares code; the recomputation below does not.)
- **C4: PASS.** E_b > 0 and below the four-band threshold at all points.
- **Separate recomputation (`src/exgc46_recompute.py`): PASS.** An independently built 7-band H (opposite phase convention, own loops, own pair basis) with a dense four-band matrix at L = 6 gives η = −1.347049, −0.038511 and +2.074241 at the three 15% points (45's U). The driver's solver at the same L and U agrees to ≤ 4.3e-8. Shared: the model definition and the U values.

## Deviations and limits

- The valence sea is still frozen, as in 45. The four bands are those above μ only.
- The 46A "fixed U" results at L = 32 compare with 45's projected L = 32 values. The L = 6 comparison in the recomputation is the only dense one.
- Mobility is t_pair = [E₀(K) − E₀(0)]/|K|² at one K. It is a finite-difference curvature, not a fitted mass.
- Not done: the finite-range interaction (Codex test 3), K* refinement beyond 45A's, and any dynamics.
