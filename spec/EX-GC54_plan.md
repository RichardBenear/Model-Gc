# EX-GC54 — where does the link shift come from? Two event-counting source rules with no free strength

2026-11-27 | EXPLORATORY | Claude | **Plan only. Nothing is run until the PI gives a separate go-ahead.** PI request 2026-11-27 ("I agree with you to focus on finding the rule for the source of link shifting instead of the suggested orbital precession test"), after the Codex review of 53 (`review/EX-GC53_codex_review.md`). Cap requested: 2 core-h (expected about 1.9; the λ = 2 point is dropped first if the cap is reached). Our construction; not LaViolette's Model G and not Chain Theory.

## Why

53 showed that a signed timing offset β on each link, sourced by the matter momentum P with GR's strength λ = 4, removes 52's preferred frame (R_dyn(0.2) = 0.997; λ* = 3.91). But λ = 4 and the choice of P as the source were put in by hand. 54 asks whether an event-counting rule fixes both, with nothing to tune.

It also answers three points from the Codex review:
- **A third strength.** λ* was interpolated between two strengths (0 and 4). 54 adds a dynamic point at λ = 2, so linearity in λ is tested rather than assumed.
- **Smooth gauge patterns.** 53A's pure-gauge loop residual (7.7 × 10⁻³) came from a random pattern that also scattered the pulse. 54 repeats it with smooth gradients.
- **A mechanism should make its own predictions.** Each rule below is zero in a stationary spin-balanced vacuum and odd in the current, and its strength comes from the counting, not from GR.

The orbital precession test that Codex suggested is parked, as the PI asked.

## What a rule must do (from 52 and 53)

- The shift must be sourced by something odd in the motion of matter, of size 4 × (momentum) × κ in the continuum, on links.
- It cannot come from the clock field N alone. **No-go (54A):** any directed offset made from the N history at the two ends of a link is, at linear order, a gradient. For example, the time-ordered coincidence excess on a link with clock rates N_x(t), N_y(t) and window w is (w²/2)(N_xṄ_y − N_yṄ_x) ≈ (w²/2)N_xN_y ∂_t(ψ_y − ψ_x), with ψ = ln N. That is a lattice gradient of ∂_tψ, so its sum around any closed loop vanishes at linear order and it is pure gauge. So the source must read matter directly, not the field it has already made.

## The counting argument (worked before writing this plan)

**Setup.** In 51, light's one-way delay across a link is set by the slowness 1/N² ≈ 1 − 2ψ. The factor 2 comes from the two-site coincidence (rate N_xN_y), and it is the same factor 2 that gives GR's light bending. ψ is sourced by κ × (the coincidence charge of matter).

**Directed links.** Now let each direction of a link count its own coincidences, (x → y) and (y → x). A firing pattern that moves with velocity v along the link is seen Doppler-shifted, by (1 ∓ v) at first order in each one-site rate. A two-site coincidence is the product of two such rates, so it is shifted by (1 ∓ v)² ≈ 1 ∓ 2v.

**Strength.** The one-way delays then become −2ψ(1 ∓ 2v). The shift is the antisymmetric part, β = 4vψ, so it is sourced by 4κ × (charge × velocity).
- This is λ = 4: 2 from two-site counting (light's factor) times 2 from the Doppler square.
- One-site counting would give 2 (the counterfactual D1 below).

**What "charge" means.** The argument leaves open what a firing carries. There are two natural readings, which agree for a point particle at linear order but differ on the lattice:
- **D2-q.** Each firing carries its coincidence charge q = ∂H/∂N (the N source of 51). The shift source is 4κ J_q, where J_q = ½{q̂, v̂} is the charge current.
- **D2-E.** Each firing is one tick of the pattern's own clock, at its frequency E. The shift source is 4κ J_E, where J_E = ½{Ĥ, v̂} is the energy current, the conserved flux of firing rate.

In both, v̂_a = i[Ĥ, X̂_a] is the velocity operator. For a continuum point particle J_E = p exactly (GR's T⁰ⁱ). J_q = p(1 + v²) agrees with it only at linear order.

**Counterfactual D1.** One-site counting gives 2κ J_q. This is what the rule would be if light were a one-site process.

## Predicted strengths (k-space, computed for this plan; no dynamic runs)

The effective strength of each rule against 53's law is λ_eff = (rule source)/(κP), averaged over the 53E packet (r = 0, σ = 5, spin up). The script is in this plan's appendix.

| m | V | D2-q: 4⟨qv⟩/⟨P⟩ | D2-E: 4⟨Ev⟩/⟨P⟩ | D1 |
|---|---|---|---|---|
| 0.3 | 0.2 | 5.73 | 3.87 | 2.86 |
| 0.6 | 0.1 | 4.66 | 3.87 | 2.33 |
| 0.6 | 0.2 | 4.70 | 3.85 | 2.35 |
| 0.6 | 0.3 | 4.75 | 3.81 | 2.38 |
| 1.2 | 0.2 | 4.11 | 3.75 | 2.05 |

**What the numbers mean.**
- **D2-q overshoots, and its overshoot depends on mass.** The excess is the packet's internal motion (Q/E = 1.14–1.38). So D2-q would make frame dragging depend on what the body is made of, which GR forbids.
- **D2-E sits close to the strength Lorentz invariance needs** (λ* = 3.91 at m = 0.6), and changes much less with mass. The rest of its gap from 4 is lattice dispersion (Ev = sin k cos k against P = sin k).
- This is a prediction made before any run. The dynamic runs below test it.

**Prediction for R.** 52 and 53 give R_dyn(0.2) = 1.115 at λ = 0 and 0.9975 at λ = 4. If R is linear in λ (tested by 54C), each rule predicts R_dyn(0.2) = 0.9975 − 0.0294 (λ_eff − 4):
- D2-q: 0.977;
- D2-E: 1.002;
- λ = 2 with 53's law: 1.056.

## Parts

**54A. The no-go, checked numerically (cheap).**
- **Field.** The 52B steady comoving field of the V = 0.2 source (regenerated, about 2 minutes).
- **Measure.** The time-ordered coincidence offset b_{x,a} = N_xN_{x+a}(∂_tψ_{x+a} − ∂_tψ_x), with ∂_t = −V∂_x. Its curl (sum around each plaquette) is compared with the curl of 53's β for the same source, over the window 8 ≤ ρ ≤ 24.

**54B. The two counting rules in full dynamics** (53E unchanged except the source).
- **Setup.** In the β law, λκP is replaced by 4κ s, where s is the link-local current of the rule. For an operator Ô, the current on link (x, x + a) is s_{x,a} = ½ Re[ψ_x†(Ôψ)_x + ψ_{x+a}†(Ôψ)_{x+a}], with Ô = ½{q̂, v̂_a} (D2-q) or ½{Ĥ, v̂_a} (D2-E). q̂, Ĥ and v̂ are built from the 51 matter stencil at the current N.
- **Runs.** m = 0.6, n = 65, separation 14, T = 16, dt = 0.1. V = 0 and 0.2, each pair with its single-packet subtraction, for each rule.
- **Measure.** R_dyn(0.2) = γF₁(0.2)/F₁(0), with each rule's own F₁(0).
- **Also reported (static, cheap):** the 53C rotating ring's source J under each rule, against P's J = 1.976.

**54C. The third strength.** 53E at V = 0.2 with 53's law at λ = 2.
- F₁(0) at λ = 2 is taken as the mean of 52's and 53's V = 0 values. This is a declared assumption: the two differ by 0.3%, so the interpolation error is at most 0.15%.
- R_dyn(0.2; λ = 2) is compared with the linear prediction 1.056.

**54D. Smooth gauge patterns for the one-clock loop** (53A code, cheap).
- **Pattern.** β = χ_{s+1} − χ_s with χ_s = A sin(2πs/Λ) on the L = 256 loop.
- **Cases.** Λ = 16, 32, 64 and 128, at A = 10⁻³; plus A = 10⁻⁴ and 10⁻² at Λ = 64. The pulse width (16) and k₀ = 0.5 are as in 53A.
- **Measure.** The two-way pulse delay difference against the uniform-shift case (a) with the same Σ|β|.

## Predictions (scored)

- **P1 (no-go).** In 54A, the plaquette curl of the time-ordered offset is below 2% of the curl of 53's β, at every ρ in the window. Allowing 2% covers terms of order ψ ≈ 0.01.
- **P2 (operators).** For each rule, Σ s over all links equals ⟨Ô⟩ to 1e-10. For the 53E packet at V = 0.2, Σ s equals the k-space value in the table within 1%.
- **P3 (D2-q fails universality).**
  - R_dyn(0.2) under D2-q is 0.977 ± 0.008.
  - Its Lorentz-invariance clause, |R_dyn − 1| < 2%, is predicted to fail.
- **P4 (D2-E).**
  - R_dyn(0.2) under D2-E is 1.002 ± 0.008.
  - |R_dyn − 1| < 2%.
- **P5 (linearity).** R_dyn(0.2) at λ = 2 is 1.056 ± 0.006, the linear interpolation of 52 and 53.
- **P6 (smooth gauge).**
  - For Λ ≥ 32, the smooth-gauge delay is below 10⁻³ of case (a).
  - At Λ = 64 it scales linearly with A, within 10%.

The ± bands are the spread of 53's own lattice effects (spin up against spin down was 0.14%, and C2 was 7 × 10⁻⁷) plus the shape difference between each rule's source and P. I estimate the shape difference at about (σ/sep)² × 5% ≈ 0.6%. That estimate is not verified.

## Checks

- **C1.** With Ô replaced by the link momentum operator, the 54B code reproduces 53E's F₁(0.2) to 1e-10.
- **C2.** D2-E at V = 0.2 with dt/2, P₁(t) for t ≥ 2, within 1e-5 (as 53's C2).
- **C3.** Each rule's source vanishes for the spin-balanced ring at ℓ = 0 (below 1e-15) and flips sign with ℓ.

## What the outcomes mean

- **P3 and P4 pass (D2-q fails, D2-E works).**
  - The shift has a counting origin. Directed two-site coincidences of the pattern's own ticks, Doppler-shifted twice, give GR's frame dragging with λ = 4 and nothing tuned.
  - The firings that drag the frame are counted by their rate (energy), not by their N charge. That is the same split GR makes between the energy current (T⁰ⁱ) and the active mass (ρ + 3p).
  - Next would be to build the directed weights into the coincidence rule itself, rather than as a source operator, and to rerun the ring benchmark with the derived law.
- **Both fail.** Two-site Doppler counting is not the mechanism, or the lattice dispersion of the current matters more than 2%. I would then look at the non-symmetric rewriting rule the PI suggested.
- **P5 fails.** R is not linear in λ, so 53's λ* = 3.91 has to be re-read.
- **P1 fails.** The N history carries a physical directed offset at the order we can see, and the no-go argument is wrong.

## Outputs

`src/exgc54.py`, `src/exgc54_recompute.py`, `src/exgc54_summary.py`, `results/EX-GC54/`, note, LOG row and HANDOFF section.

## Appendix: the k-space table

Run from `src/` with the 52 helpers (r = 0, σ = 5, n = 96). The rows are the packet's k distribution w(k). Per k:
- q = ∂E/∂N at N = 1;
- v = ∂E/∂k_x;
- P = sin k_x.

D2-q = 4Σw q v / Σw P, D2-E = 4Σw E v / Σw P, and D1 is half of D2-q.

```python
import numpy as np, exgc52 as E52; G = E52.G
for m in (0.3, 0.6, 1.2):
    for V in (0.1, 0.2, 0.3):
        k0 = E52.k_for_v(V, m, 0.0); psi, K, w = E52.packet_k(96, m, 0.0, (48, 48, 48), E52.SIG, (k0, 0, 0))
        h = 1e-5; e1 = np.array([h, 0, 0])
        q = (G.disp(K, 1 + h, m, 0) - G.disp(K, 1 - h, m, 0)) / (2 * h)
        v = (G.disp(K + e1, 1, m, 0) - G.disp(K - e1, 1, m, 0)) / (2 * h); E = G.disp(K, 1, m, 0)
        P = np.sum(w * np.sin(K[..., 0])); print(m, V, 4 * np.sum(w * q * v) / P, 4 * np.sum(w * E * v) / P)
```

## Amendment A1 (2026-11-27, written with the code, before any scored run)

The operator check (P2) ran first, with the final code, and is scored as registered (`results/EX-GC54/EX-GC54_ops.json`). It showed something I did not foresee when writing the plan:
- **The registered D2-q operator is exactly twice the energy current.** At uniform N, Q̂ = Σ_x ∂Ĥ/∂N_x = 2Ĥ − mβ̂, and the Dirac mass matrix β̂ anticommutes with the velocity v̂ = i[Ĥ, X̂], so ½{Q̂, v̂} = 2 · ½{Ĥ, v̂}. Measured: Σs(D2-q)/Σs(D2-E) = 2.0000, so the D2-q operator has λ_eff = 7.70, not 4.70.
- **Why the plan's table differs.** The table's D2-q column used the band-diagonal product q(k)v(k). The operator also has interband (Zitterbewegung) parts, which the table left out. D2-E has none, because Ĥ commutes with the band projector, and it matches the table to 1.5 × 10⁻¹⁰.
- **Consequences, recorded and not re-scored.** P2's 1% k-space clause fails for D2-q (0.64). P3's number (0.977 ± 0.008) was computed from the wrong λ_eff and will be scored as registered. The registered operator is the one that runs. Because D2-q = D2-E at strength 8 (up to O(ψ) differences in the bond weights), its run is also a far point for the linearity test (P5). The linear prediction for it, written now, is R_dyn(0.2) = 0.9975 − 0.0294 × 3.70 = 0.889 (reported, not scored).

Implementation choices fixed now:
1. **Source refresh.** For the two rules the link current is recomputed once per time step from the state at the start of the step and held over the four RK stages. This keeps the cost near the cap. C2 (dt/2) measures the effect. For 53's law (P) the source is evaluated at every stage, as in 53E.
2. **C1.** 53's law run through the 54 code at V = 0.2 to T = 8 is compared with 53E's stored P₁(t) for t ≤ 8, to 1e-10.
3. **C2.** D2-E at V = 0.2 with dt = 0.05 to T = 8, against the dt = 0.1 run, P₁(t) for 2 ≤ t ≤ 8, within 1e-5 relative to max |P₁|.
4. **Ring (C3 and the reported J)** at n = 129, uniform N, ℓ = +2, 0, −2, spin-balanced as in 53C.
