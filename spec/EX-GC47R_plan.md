# EX-GC47R — corrections to EX-GC47 after the Codex review

2026-11-27 | EXPLORATORY | Claude | Plan only, pushed before any run. PI request 2026-11-27 13:32 UTC ("the 47 review is available. Take a look. You probably should make the corrections it suggests."), taken as the go-ahead for the two corrections and the convergence checks in the [Codex review of GC47](../review/EX-GC47_codex_review.md). The review's recommended next experiment (a massive band with 2 t_pair E_rest = c_signal²) is **not** part of this; it needs its own plan. 47's registered scoring is not changed; this adds a corrected record. Cap 2 core-h. Our construction; not LaViolette's Model G, not Chain Theory.

## Corrections

1. **Release packet seam.** 47C built the launch Gaussian from coordinates wrapped about the origin, so for L = 24, start 8, σ = 3, the site x = −12 got displacement −20 instead of +4 (about 61% of peak weight there). Fix: build the packet and track its centroid periodically about the packet's own centre. Repeat at L = 24 and 32, σ = 2 and 3, starts 6 and 8, fit windows T and T/2, dt and dt/2. Report energy drift ⟨H_eff⟩(t) − ⟨H_eff⟩(0).
2. **Strain (and lapse) form factor from the derivative of the full projected operator.** In the band-projected model the projector depends on the field, so the local density that the field couples to is F(q) = ⟨Ψ_q|P o(q) P|Ψ₀⟩ plus the projector-derivative terms ⟨Ψ_q|V Q R o(q)|Ψ₀⟩ + ⟨Ψ_q|o(q) R Q V|Ψ₀⟩ (per carrier, single-carrier denominators). Two variants: (i) Q over all other bands (exactly the derivative of the projected model), (ii) Q over the three other bands above μ only (Pauli-allowed with the frozen sea). Recompute V₁₂ and the release with (i).
3. **Force and box sweep.** Force F(r) = −ΔV₁₂/Δr along x for the microscopic pairs at L = 24, 32, 48, against the point-charge periodic force at the same L and the infinite-box 1/r² force. Report r²F/(G_eff·coefficient) against r.

## Predictions (scored)

- **Q1.** Variant (i) at q = 0 equals t∂E/∂t of the projected pair to 1e-7 relative, for B1 and B2; the lapse correction at q = 0 is zero to 1e-12.
- **Q2.** With corrected F, V₁₂ changes by less than 0.1% from 47 at L = 32, r = 3–12.
- **Q3.** Corrected release: acceleration equals −2 t_pair⟨∂V₁₂/∂x⟩ within 5% for every (L, σ, start, window); a/a_N stays within 10% of 47's 1.16–1.27 × 10⁻⁴; energy drift below 1e-10 relative to the field energy scale |⟨V₁₂⟩|.
- **Q4.** The B1/B2 acceleration ratio stays within 0.90–0.97 (47: 0.915–0.945).
- **Q5.** The microscopic force equals the periodic point-charge force within 0.5% at r = 4–L/3, at each L; r²F approaches the infinite-box value as L grows (descriptive).

## Outputs

`src/exgc47r.py`, `results/EX-GC47R/`, and a correction section appended to the 47 note.
