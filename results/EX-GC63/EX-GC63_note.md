# EX-GC63 — the directed signal count, per local tick: the moving-mass shift with nothing tuned

2026-11-28 | EXPLORATORY | Claude | PI closure pending | [Plan](../../spec/EX-GC63_plan.md) (pushed before any run; PI go-ahead 17:47 UTC on the decision card) | 0.92 core-h on the cloud container (Linux, Python 3.13.18, numpy 2.5.3; 4 processes). Scoring: `src/exgc63_summary.py` → `summary.json`; recomputation: `src/exgc63_recompute.py` → `recompute.json`. This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**Counted per local tick, 54's directed signal-count rule keeps the co-moving pair Lorentz-invariant with nothing tuned: R_dyn(0.2) = 1.0028, against the predicted 1.002 ± 0.008.** That closes the gap 54 left (1.019). The gap was the clock-well weighting of the counted current, as 54's post-hoc analysis said.

The count does **not** carry a resting spin's source. Its magnetic moment for 59's polarised packet is 9.2% of the PS rule's, so P5 (< 5%) fails as registered. But the reason confirms the prediction's physics rather than contradicting it. The counted source's moment is 0.93 of the canonical momentum P's own moment, and P itself carries only 9.9% of PS on this lattice (the lattice-orbital share 57 found). The counted source has no spin term of its own. The ¼Cᵀs plaquette term (58) remains a separate ingredient.

So, in 3D, **the moving-mass part of λ = 4 is now a count**: two-site signal weight (2) × Doppler square at the two endpoints (2), counted in the endpoints' own ticks. The spin part is not.

| | What was scored | Result | Verdict |
|---|---|---|---|
| P1 | in-place strength in the pair 3.85 ± 0.02 | 3.852 (54's in-place 3.368) | PASS |
| P2 | R_dyn(0.2) = 1.002 ± 0.008 | 1.0028 (line through 52/53 at 3.852: 1.0018) | PASS |
| P3 | \|R_dyn − 1\| < 2% | 0.28% | PASS |
| P4 | F₁(0) within 1% of 53's 4.699 × 10⁻⁴ | −0.007% | PASS |
| P5 | resting spin: \|m(count)\| / \|m(PS)\| < 0.05 | 0.0916 (both spins) | FAIL |

## Checks

| | | | |
|---|---|---|---|
| C1 | zero packet → zero source | 0 | PASS |
| C2 | source total flips with V, relative < 1e-12 | 1.5e-11 | FAIL (rounding level; see deviations) |
| C3 | N ≡ 1: count = 54's D2-E | 0 | PASS |

**Recompute.**
- R from an independent least-squares design matrix: 1.0027651, the same as the scorer to 10⁻¹⁵.
- Fitting only t ≥ 2 gives 1.0022.
- The z-moments from an explicit x j_y − y j_x sum reproduce 0.0916.
- Shared dependencies: run records, packet builders, operators.

## Deviations

1. **C2's first construction was wrong.** `k_for_v` bisects on k > 0 only, so `packets(−0.2)` built a packet at rest. Before the scored runs I replaced it with the +V packet's k₀ negated. The rerun residual of 1.5 × 10⁻¹¹ misses the 10⁻¹² gate. The total is a sum over 65³ links of terms of order 10⁻³, so this is plausibly summation rounding, but I have not shown that.
2. **P5's threshold assumed the lattice canonical P carries no moment for a resting packet.** That was wrong on this lattice: it carries 9.9% of PS (57 found the same orbital share). Post hoc, not re-scored: count / P = 0.929, close to the strength ratio 3.85/4 = 0.963.

## What it means

- **The orbital shift is a count.** With the directed signal count measured in local ticks, the shift a moving mass sources is fixed with no free strength, to 0.3% of Lorentz invariance in the pair. 52 (no shift) was 11.5% off, and 54 (coordinate-time count) 1.9%.
- **Spin needs its own event.** The counted energy current is ½{Ĥ, v̂}, which in the continuum is the canonical momentum (no spin curl). GR's source is the Belinfante current, P + ½∇×S. The extra half lives on plaquettes (58), and this study does not derive it from counting.
- **The remaining 3.85 vs 4 is the lattice difference between E·v and p.** It moves R by about 0.004, inside the band.
- **Still put in:** the two-site weight (37B: 1 + γ = 2 in 3D only), the Doppler square (EX-GC62 1b, an argument, not a derived rule), the ¼Cᵀs spin term, and the field law. The D study is where the two-site weight would be tested.
- **Suggested next (not approved):**
  1. A counting origin for the spin term: whether plaquette circulation events counted per local tick give ¼Cᵀs.
  2. The D study.
  3. Score the per-local-tick count against 59's precession with the PS spin term kept.
