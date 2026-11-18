# Gates — Model Gc

Registered criteria, written and committed before the runs that test them, followed by outcomes and closures. Format
as in ModelG `GATES.md`.

## Inherited from ModelG (copied verbatim, source commit edab7d47)

## MG195 REGISTERED — inverse-square attraction in the conserved-G test model (3D axisymmetric)

PI approval: Richard Benear, 2026-11-18 10:30 America/Denver ("yes", adopting D1–D4). Specification `spec/MG195_inverse_square.md`, frozen with this commit before computation; implementation `src/mg195_inverse_square.py`. Run by Claude in the cloud, cap 2 core-hours, at most 4 processes. Model is our EX008 construction, not LaViolette's.
- Points P1 ("G feeds X", k = 0.2), P2 ("G feeds Y", k = 0.2), P3 ("G feeds X", k = 0.05); pairs at d = 18, 22, 27, 33, 40. Prediction v = μM/(4πd²) from the R = 120 radial solve.
- Per point SUPPORTED iff residual ≤ 1e-8, v < 0 above the null floor (10·max|v_C1|) at all d, |v/pred − 1| ≤ 0.05 for d ≥ 22, exponent over d = 22–40 within 2.00 ± 0.10. UNRESOLVED on non-convergence or relative grid/far change > 2% (P1, d = 22). Otherwise NOT SUPPORTED.
- C2 (λ(G∞ − G), λ = 0.01): v(40) within a factor 2 of the Yukawa prediction μM/(4π)e^(−κd)(1/d² + κ/d) and below 0.5 × the P1 inverse-square prediction. Study verdict INVERSE-SQUARE CONFIRMED iff P1–P3 SUPPORTED and C2 passes; otherwise per point.
- Pre-registration corrections to the draft (C2 exchange form and its 0.092 factor; floor split into absolute and relative parts) are listed in the spec header.

## MG195 COMPLETE — P1 SUPPORTED, P2 SUPPORTED, P3 UNRESOLVED, C2 PASS; study verdict per point (PI closure pending)

Report `results/MG195_inverse_square/MG195_results.md`. All 27 pair solves converged (residual ≤ 7.4e-12); grid and far-boundary changes 0.07% and 0.14%. P1 and P2: v/prediction 0.978–0.990 for d ≥ 22, exponents 1.984 and 1.981. P3: null floor 4.7e-6 (from C1 at d = 18; the short-range Turing interaction was larger than the spec assumed) exceeds all P3 predictions, so UNRESOLVED; P3 also ratio 0.927 at d = 22 and exponent 1.883. C2: v(40) / Yukawa prediction = 1.003, and 0.096 of P1. Post-hoc C1 subtraction (not scored) gives P3 ratios 0.995–0.997 and exponent 2.003. Separate recomputation agrees (speeds to 1e-6; M by flux to 1e-4). About 0.5 core-hours of 2.

## MG195 PI CLOSED — 2026-11-18

PI Richard Benear, 2026-11-18 11:06 America/Denver: "go ahead and sign off." Closed as reported: P1 SUPPORTED, P2 SUPPORTED, P3 UNRESOLVED (null floor set by the short-range Turing interaction), C2 PASS; study verdict per point. No result changed.
