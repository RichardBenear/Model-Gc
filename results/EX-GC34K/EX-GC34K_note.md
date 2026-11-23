# EX-GC34K — the first-order gravitomagnetism error: where it comes from, and whether another re-slicing gradient removes it (EXPLORATORY)

**PI CLOSED 2026-11-23** (Richard: "Close 34K and 34L", after ChatGPT's follow-up review `review/EX-GC34J_followup_chatgpt_review.md`).

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "You can close them and use any inputs you get from the review
for the next test decision. I don't have any bias toward the 3 options. Maybe you should do all 3."
- The plan `spec/EX-GC34K_plan.md` and the `explore/LOG.md` entry were pushed before any run (66f3ab8).
- Code: `src/exgc34k.py`.
- Registered outputs:
  - `diag` and `diag_g1` → `diag.json`;
  - `family` → `family.json`, `derive.json`, `kernels.npz`;
  - `e4` → `e4.json`;
  - `spec` → `candidates.json`;
  - `sync` → `sync.json`;
  - `check` → `check.json`.
- Post hoc output: `posthoc_zb` → `posthoc_zb.json`.

This is our construction, not LaViolette's Model G.

## Conclusion

**No candidate passes as registered.** J1 is the only scored candidate, because no other gradient in the family has
a local solution with better than first-order error. J1 passes everything except the gravitomagnetism benchmark (E4)
and its order test. That now includes a clean, registered synchronous-gauge time domain (G4″).

What each of the three options gives:

- **(a) Finer lattice.** This works, at first order. J1's E4 error is kh/4 at every tested spacing:
  - 3.1e-3 at kh = 0.0125;
  - 7.5e-3 at kh = 0.03;
  - 1.0005e-2 at kh = 0.04;
  - 1.25e-2 at kh = 0.05.

  The largest tested kh with E4 ≤ 1e-2 is 0.03. At kh = 0.04 it misses by 0.05%, so the threshold sits at about
  kh = 0.04.
- **(b) Staggered shift and (c) centred gradient.** Neither works. They are the same thing in the identity (see the
  plan), and the centred gradient has no solution at range 1, 2 or 3. The consistency residual barely falls with
  range (0.32, 0.17, 0.12).

  Post hoc I found why, and it is not a question of range. At every zone-boundary point (each q component 0 or π, not
  all 0), the centred gradient vanishes, so the re-slicing strain is zero. But the lapse–edge coupling on the right of
  the identity is not zero there (|C_en| ≈ 20–40). No bounded kinetic term of any range, local or not, can satisfy
  the identity at those points. **The averaged staggered shift and the centred gradient are excluded outright.**
- **In between.** The other one-sided mixtures (γ = −0.25 and −0.75) and the Kuhn-edge least-squares gradient have no
  local solution at range ≤ 2 (−0.25 also at range 3). Their residuals do fall steadily with range:
  - γ = −0.25: 0.11, 0.015, 3.6e-3;
  - γ = −0.75: 0.028, 4.9e-3;
  - Kuhn: 9.9e-3, 9.8e-4.

  So non-local, decaying solutions probably exist. But they are not local kernels, and their E4 order is unknown.
- **The two that exist.** Only the pure forward (γ = 0) and pure backward (γ = −1) gradients have local solutions.
  Both are first order: E4 ratio 2.00, 1.25e-2 at kh = 0.05.

**Part A (the review's priority ii) narrows 34J's claim.**
- The worst case is a source moving along a lattice axis (q along (1,0,0), velocity parallel), in one Riemann
  component, (01|12). There c1 = 0.250. Every other case is at least 2.6× smaller.
- **The leading coefficient is not fixed** within the identity-obeying space. Over 33 exact solutions (30 random at
  ranges 1 and 2, plus J1, J2 and J3), the worst component's c1 ranges from 0.10 to 0.33, a spread of 91%. So the
  verdict is **A-variable**, and 34J's "the coefficient is forced by the identity" was too strong.
- **But the first order itself did not change in any sample.** All 33 have an order ratio of 2.00.
- **No sample beat J1 on the maximum over cases.** The smallest maximum was 0.252, against J1's 0.250. Lowering the
  worst component pushes the error into another case or component.
- The smallest-c1 kernel keeps G1 on the reduced direction set: no growing root at 56 cases (14 directions × kh 0.05,
  0.2, 1, 3).
- All 9 source and readout placements are first order (ratio 2.0). The best ones (0.250) include the vertex–vertex
  one already used.

## Registered results

**Part A**

| item | result |
|---|---|
| A1 worst case | q ∥ (1,0,0), v parallel; component (01\|12); c1 = 0.2500 |
| A1 next cases (c1, max over components) | (0.8,0.5,0.33) v ∥: 0.038; (1,0,0) oblique: 0.096; rest ≤ 0.016 |
| A2 spread of worst-component c1 (33 kernels) | 0.914 → **A-variable** |
| A2 spread of max-over-cases c1 | 0.58 (range 0.252–0.417; J1 0.250) |
| A2 order ratio, all 33 | 2.000–2.003 |
| A2 smallest c1 (0.101, range 1, scale 0.3): G1 reduced set | PASS (0 growing, 0 unconfirmed) |
| A3 nine placements | all ratio 2.0; c1_max 0.25–0.63 |

**Part B** (existence: consistency residual ≤ 1e-10)

| gradient | R=1 residual | R=2 | R=3 | exists | E4 at kh 0.05 / ratio |
|---|---|---|---|---|---|
| γ = 0 forward | 3.0e-15 (dim 123) | – | – | yes | 1.251e-2 / 2.00 |
| γ = −0.25 | 0.107 | 0.0153 | 0.0036 | no | – |
| γ = −0.5 central (≡ staggered) | 0.317 | 0.168 | 0.123 | no | – |
| γ = −0.75 | 0.028 | 0.0049 | – | no | – |
| γ = −1 backward | 3.6e-15 (dim 177) | – | – | yes | 1.252e-2 / 2.00 |
| Kuhn-edge LS | 0.0099 | 0.00098 | – | no | – |

No member qualifies for Part C, which needs a ratio of at least 3.

**Part C (J1)**

| criterion | J1 |
|---|---|
| I0 identity exact | PASS (1.2e-14, family `g0_R1` re-derivation; 34J record 1e-14) |
| G1 no growing root, full 34F/34G set | PASS (0, 0 unconfirmed) |
| G2 two waves at c | PASS (max 7.5e-4) |
| G3 chains = exact symmetries | PASS |
| G5 static, gauge | PASS (0, 1.4e-16) |
| G5 E4 ≤ 1e-2 at kh = 0.05 | **FAIL** (1.251e-2) |
| E4 order ratio ≥ 3 | **FAIL** (2.00) |
| G4″ synchronous time domain, clean seeds | **PASS** — both directions, both steps: constraints ≤ 7.8e-12; wave amplitude change ≤ 1.4e-10; non-zero seeds 4.1–4.6e-8 of the wave; zero-sector seeds ≤ 6.8e-13; 5 of 5 non-zero modes and 4 of 8 zero-cluster directions satisfy the constraints |
| G7 real space | PASS (5.6e-16) |
| option (a): E4 at kh 0.03 / 0.04 | 7.50e-3 / 1.0005e-2; largest tested kh ≤ 1e-2 is 0.03 |

**My stated expectations:**
- A-fixed: wrong. The coefficient varies, but the order and the best maximum did not improve.
- Forward and backward exist, central does not: right. Central fails at range 3 as well, and the reason is a zero of
  the gradient, not too short a range.
- The intermediate members have no local solution.
- J1 passes G4″ and fails E4: right.

## Post hoc (labelled)

`posthoc_zb.json` checks all 7 non-zero zone-boundary points (q components 0 or π) for γ = −0.25, −0.5 and −0.75:
- For γ = −0.5, |g| ≤ 6e-16 and |σ_ξ| ≤ 7e-16, while |C_en/l0| = 20–40. The identity K σ_ξ = C_en/(μ l0) cannot hold
  there for any bounded K.
- γ = −0.25 and −0.75 have |g| = 5/h there. They have no such zero, which is consistent with their residuals
  shrinking with range.

## What this means

- **The re-slicing-exact route gives a first-order lattice graviton, and no gradient choice tested here changes
  that.** The only gradients with local exact kernels are one-sided. Every exact kernel tested (33 in A, 2 in B) is
  first order in the moving-source benchmark, with a maximum coefficient no better than J1's 1/4.
- **The staggered and centred ideas are dead in this form.** Their gradient vanishes at the zone boundary, where the
  lapse coupling does not. A staggered formulation could only work if the lapse coupling also changed, i.e. a
  different potential term, not just a different shift placement.
- **Option (a) is the practical one.** J1 meets the 1% benchmark at kh ≤ 0.03. It is clean in every other respect,
  now including a registered time-domain pass.
- **Open, not tested, none approved:**
  1. The coefficient varies linearly with the kernel. A targeted minimisation of the maximum over cases within the
     range-2 exact space (not random sampling) could lower 1/4. It could not change the order unless the O(kh) term
     can be zeroed in every case at once.
  2. Non-local (rational) kernels for the in-between gradients, whose residuals shrink with range.
  3. Change the potential side (the lapse coupling) so that a centred gradient becomes admissible.

## Recomputation and shared dependencies

- J1's E4 table in `e4.json` matches 34J's independently recorded `e4.json` to all digits (1.2509733e-2 at
  kh = 0.05). Both use `exgc34e.e4_case`, so this is a reproduction, not an independent check.
- The range-1 forward solution re-derived here (`g0_R1`) has the same space dimension (123) and distance from T
  (0.657) as 34J's J1.
- `diag_g1` regenerated the smallest-c1 kernel by replaying the random stream. Its max-over-cases c1 matched the
  recorded value exactly (0.386665).
- The zone-boundary obstruction was checked by direct evaluation, a separate path from the SVD consistency solve.

## Compute

About 0.53 core-hours, against a 2 core-hour cap and an expected under 1:
- diag: 29 s;
- diag_g1: 7 s;
- family: 1199 s, mostly the two range-3 solves at about 550 s each;
- e4, spec, sync and check: about 30 s.

That total includes a first family attempt (about 600 s) that the cloud background time limit killed before it wrote
anything.

Platform: Linux cloud, Python 3.13, numpy 2.5, scipy 1.18.

## Deviations

- **The plan's A2 follow-up (G1 on the reduced set for the smallest-c1 kernel) was missing from the first `diag`
  code.** I added it as `diag_g1` after `diag` had run, with the registered rule unchanged. It regenerates the kernel
  from the same seeded stream and checks the replay against the recorded value.
- **I0 for J1 in Part C** is taken from the re-derived identical space (`g0_R1`, 1.2e-14) and the 34J record. There is
  no separate `I0` call on the reloaded J1 coefficients.
- **One killed run**, described under Compute. It produced no output, and the rerun is the recorded one.
