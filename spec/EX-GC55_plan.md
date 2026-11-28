# EX-GC55 — one interaction for the link shift, and a clock readout that respects gauge

2026-11-28 | EXPLORATORY | Claude | PI request and go-ahead 2026-11-28 00:15 UTC ("close 54 and write up 55 and run it"), after the Codex review of 54 (`review/EX-GC54_codex_review.md`). Cap: 2 core-h (expected about 1.5). This plan is pushed before any run. Our construction; not LaViolette's Model G and not Chain Theory.

## Why

The Codex review of 54 made two points that 55 takes up.
1. **Source and response were different operators.** In 54 the shift was sourced by the energy current but felt by matter through the momentum P (53's coupling). That is not one interaction, so energy need not be conserved, and the 54 force gap cannot be read cleanly.
2. **The clock control failed,** and choosing a clock site where β = 0 is not a fix.

## What I found while writing this plan (before any run)

**The 53 light rule is not the covariant wave equation.** The 53 loop code uses D_t² u = ∇²u with D_t = ∂_t − β∂_s (the "advected" form).
- **The covariant equation.** For the 1D metric ds² = −dt² + (ds + β dt)², the covariant wave equation at linear order is u_tt = u_ss + 2βu_st + (∂_sβ)u_t. The advected form lacks the last term.
- **Why that matters for gauge.** For a pure-gauge pattern β = ∂_sχ, that term is χ″u_t, and it is exactly what a clock resynchronisation t → t − χ(s) produces. So the 53 light rule is not gauge covariant, which explains 53A's and 54D's gauge residuals.
- **The covariant form as a Hamiltonian.** It is the flow of H = Σ[½Nπ² + ½N_sN_{s+1}(Δu)² + β π ∂_su]. The shift term is ADM's −β·(field momentum), the same structure as matter's −β·P.
- In first-order form, u_t = Nπ + β∂_su and π_t = (lattice Laplacian) + ∂_s(βπ). The 53 code has β∂_sπ in the last term. The two differ by (∂_sβ)π.

**What one interaction means for matter.** For each choice of current operator Ĵ, take the coupling H_β = −Σ β Ĵ in the matter Hamiltonian, and the source j = −∂⟨H⟩/∂β = ⟨Ĵ⟩. With the β field action L_β = −(1/2λκ)Σ[(∂_tβ)² − (∇β)²], the field law is the 53 law, and the total
E_tot = ⟨H_m⟩ − Σβ⟨Ĵ⟩ − (1/2λκ)Σ[(∂_tβ)² + (∇β)²]
is conserved (N held at 1).
- The sign of L_β is the "wrong" one, as GR's g₀ᵢ sector is: parallel mass currents must repel to cancel 52's excess.
- 53's law (Ĵ = P̂ for both source and coupling) is already of this kind. 54's mixed rule (source ½{Ĥ, v̂}, coupling P̂) is not.
- The new candidate is Ĵ = ½{Ĥ, v̂} for both: the reciprocal energy-current rule (RJE).

## Parts

**55A. The covariant clock (cheap).**
- **Setup.** The 53A/54D loop (L = 256, k₀ = 0.5, pulse width 16, one-clock Sagnac delay) with the covariant first-order rule.
- **Initial data.** The initial momentum is π = (u_t − βu_s)/N, with u_t taken from the right-moving pulse (the same pulse expressed in the gauge-transformed variables).
- **Centroid.** The pulse centroid uses the conserved energy density e = ½Nπ² + ½(gradient energy) + βπ∂_su.
- **Gauge cases.** Smooth gauge β = χ(s+1) − χ(s) with χ = A sin(2πs/Λ + φ):
  - Λ = 16, 32, 64 and 128 at A = 10⁻³ and φ = 0;
  - φ = 0, π/4, π/2 and 3π/4 at Λ = 64;
  - A = 10⁻⁴ and 10⁻² at Λ = 64.
- **Positive controls.** Uniform β = 10⁻³ (53A case a) and the one-sided pattern (53A case c).
- **Real circulating shift.** The 53C/D rotating ring, ℓ = ±2, n = 129: wave delay with the energy centroid on the enclosing loops of half-width 20, 28 and 36, against the readout 2Σβ/N⁴ × the lattice factor cos k₀/cos²(k₀/2) = 0.935.
- **Also reported, not scored:** both ingredients separately (covariant rule with the old initial data; advected rule with the new initial data).

**55B. Energy bookkeeping for three coupling choices (cheap).**
- **Setup.** One packet (m = 0.6, σ = 5, spin up) at V = 0.2. n = 49, N frozen at 1, β dynamic from its steady initial state, κ as in 53E, λ = 4, T = 8, dt = 0.1, source evaluated at every RK stage.
- **Rules.** P/P (53), J_E-source/P-coupling (54 mixed) and J_E/J_E (RJE).
- **Measure.** E_tot(t) as defined above, using the coupling current in the interaction term. The drift is max|E_tot(t) − E_tot(0)| over |Σβ j| at t = 0.

**55C. The reciprocal rule in full dynamics.** The 53E co-moving pair with N, β and both packets dynamic, at n = 49 to keep within the cap.
- **Geometry.** Packets at y = c ± 7, about 17 sites from the faces. Separation 14, T = 8, dt = 0.1.
- **Runs.** V = 0 and 0.2, each pair with its single-packet subtraction, for RJE and for the P/P reference at the same n and T.
- **Measure.** R_dyn(0.2) = γF₁(0.2)/F₁(0) for each rule, F₁ fitted as in 53E. The box and T are not 53's, so only the difference R_RJE − R_P is compared with predictions.

## Predictions (scored)

- **P1 (gauge).** In 55A, every smooth-gauge case (all Λ, φ and A) gives a delay below 10⁻³ of case (a), with the same Σ|β| normalisation as 54D.
- **P2 (controls).**
  - Uniform β: delay/readout = 0.935 within 0.5%, the lattice eikonal factor.
  - One-sided pattern: ¼ of (a) within 2%.
  - Ring: wave delay/readout = 0.935 within 1% on all three loops. 53 had 0.944–0.955 with the advected rule.
- **P3 (bookkeeping).**
  - P/P and RJE conserve E_tot: drift below 10⁻⁴ of |Σβj|.
  - The 54 mixed rule does not: drift above 10⁻³ of |Σβj|.
- **P4 (reciprocal energy current).**
  - R_RJE − R_P = +0.030 ± 0.012. The reason: in the packets' own N wells, J_E is about 0.84 (pair) and 0.88 (single) of its uniform-N value (54). With RJE the reduction enters both the source and the coupling, so the shift's effect on R falls to about 0.74 of P's.
  - So the reciprocal energy current is predicted to fail Lorentz invariance by about 3%, and the momentum P (53's law, already reciprocal) to remain the consistent operator.
- **P5 (reference).** R_P at n = 49, T = 8 lies within 1% of 53E's 0.9975. This is a size check on the smaller box.

## Checks

- **C1.** The 55A covariant code with β = 0 reproduces the 53 loop delay for β = 0 exactly (the two rules coincide there).
- **C2.** 55B P/P at dt = 0.05: the drift falls by at least 10×, confirming that the residual is time-step error.
- **C3.** The RJE coupling is Hermitian, ⟨φ|H_βψ⟩ = ⟨H_βφ|ψ⟩ to 1e-12 for random φ and ψ, and ⟨ψ|H_βψ⟩ = −Σβs to 1e-12.

## What the outcomes mean

- **P1 and P2 pass.**
  - The 53/54 clock residuals came from the advected light rule, not from the readout location.
  - The covariant rule is the one to use, and it puts light's and matter's couplings to the shift in the same ADM form (−β × momentum).
  - The 53 ring numbers move by the measured amount, and I would restate them.
- **P3 passes.** Reciprocity is now an explicit, tested property. Only P/P and RJE are admissible.
- **P4 as predicted (RJE fails by about 3%).**
  - The operator that both sources and feels the shift should be the momentum, as in GR's ADM form.
  - The counting argument then supplies the strength λ = 4, and 54's energy current is its uniform-field shadow.
  - Next would be to derive the momentum coupling from directed coincidence weights on the hops, so that the rule carries its own shift.
- **RJE comes out near 1 instead.** The well effects cancel in a way I did not foresee, and RJE becomes the leading candidate.

## Outputs

`src/exgc55.py`, `src/exgc55_summary.py`, `src/exgc55_recompute.py`, `results/EX-GC55/`, note, LOG row and HANDOFF section.
