# EX-GC34O — a Y1 scalar stress that conserves momentum exactly (EXPLORATORY)

Claude, 2026-11-24. PI go-ahead (Richard, 2026-11-23): "Yes, go ahead with 34O".
Source: `review/EX-GC34N_chatgpt_review.md` (resolve the scalar momentum mismatch before any coupled run).

This plan and the `explore/LOG.md` entry are pushed before any 34O computation. Cap: 2 core-hours on the cloud; I
expect a few CPU-minutes. Our construction, not LaViolette's Model G.

## What is held fixed

- The Y1 lattice and vacuum action of 34N, unchanged.
- The flat scalar evolution: φ̈ = Δ_lat φ − m² φ, with the 7-point staggered Laplacian.
- The scalar's lapse and shift couplings (so the energy source s_n and momentum density s_β are 34N's).

Only the scalar's couplings to the spatial metric (h_aa at sites, h_ab on plaquettes) are rebuilt.

## Derivation (done before running; the runs test it)

Write G_a(v) = (φ(v + e_a) − φ(v))/h on link v + e_a/2, and Δ_a f = f(v + e_a) − f(v).

The time derivative of the 34N momentum density s_β,a = −(V/2)(φ̇_v + φ̇_{v+a}) G_a(v), with the flat evolution,
splits into exact lattice differences:
- the parallel gradient part gives −(V/2h) Δ_a P_a, with **P_a(v) = G_a(v) G_a(v − e_a)** at site v;
- each transverse part (b ≠ a) uses the exact lattice identity Δ_b G_a = Δ_a G_b on every plaquette. It gives
  (V/2h)[Δ_a X_b − (Y_ab(v) − Y_ab(v − e_b))], with X_b(v) = ½(G_b(v)² + G_b(v − e_b)²) at sites and
  Y_ab = 2 Ḡ_a Ḡ_b on plaquettes (bars are the averages of the two bordering links).

Exact local balance then requires:
- normal stress s_aa = (V/4)[φ̇² − m²φ² + P_a − Σ_{b≠a} X_b];
- shear stress s_ab = V Ḡ_a Ḡ_b (already 34N's).

34N used X_a in place of P_a in the parallel term. The difference, X_a − P_a = ½(G_a(v) − G_a(v − e_a))², is the
second-order residue 34N measured, and is ChatGPT's 1D difference.

**The action (Y1-S2).** It is 34N's action with one change. In the link gradient term the parallel metric coupling
moves to the site in product form:
- links: −½(1 + n̄ + ½ Σ_{c≠a} h̄_cc) G_a² − ½(φ̇_v + φ̇_{v+a}) β_a G_a;
- sites: ½(1 − n + ½ tr h) φ̇² − ½(1 + n + ½ tr h) m² φ² + **¼ Σ_a h_aa(v) G_a(v) G_a(v − e_a)**;
- plaquettes: h_ab Ḡ_a Ḡ_b.

The continuum limit is the same: ¼ h_aa G_a² = ½ h_aa G_a² − ¼ h_aa G_a². The lapse and shift terms are
unchanged, so the energy identity should stay exact.

## Gates

- **O0, 1D algebraic control.** A 1D periodic chain, random φ and φ̇ (N = 64, h = 1 and 0.5). The local balance
  ∂_t p + Δ T / h with T = (φ̇² − m²φ² + G_v G_{v−1})/2 (the 1D reduction):
  - product form ≤ 1e-12 relative;
  - average-of-squares form (34N) visibly non-zero (≥ 1e-3), as the null control.
- **O1, derivation check (M0).** Finite-difference variation of the directly coded Lagrangian against Σ s·x,
  Richardson extrapolated, ≤ 1e-8, at h = 1 and 0.5.
- **O2, energy identity (M1).** ≤ 1e-10 relative, on the 34M packet (axis and generic; h = 1, 0.5, 0.25;
  t = 0, 6, 12), and on random 3D fields.
- **O3, momentum identity (M2).** EXACT means ≤ 1e-10 relative on the same packets and on random 3D fields (N = 8,
  h = 1 and 0.5). Anything else is FAIL.
- **Null control.** 34N's stress reproduces its registered r_P (to 2 digits) in the same harness.

**Overall:** PASS if O0–O3 and the null control hold.

**Reported, not gated:**
- the integrated stress tensor against the continuum (as in 34M M3);
- whether the product-form gradient energy changes sign near the cutoff (the zone-boundary mode gives P_a = −G²).

## My expectations, stated before running

Everything passes, with O3 at rounding (about 1e-15), because the balance above is an algebraic identity in the
lattice fields. The risk is a slip in signs or link offsets in the transverse terms; the random-field tests catch that.

## Not in this study

- The coupled matter–gravity run (review step 2). It becomes a separate study if 34O passes.
- The knot; any change to the white paper.

## Files

`src/exgc34o.py run` → `results/EX-GC34O/results.json`. It reuses `exgc34n` (layout, 34N sources, identities) and
`exgc34m.packet`.
