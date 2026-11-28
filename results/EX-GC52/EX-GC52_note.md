# EX-GC52 — a moving source in the 51 model: does a single coincidence-carried scalar pick out a preferred frame?

2026-11-27 | EXPLORATORY | Claude | PI CLOSED 2026-11-27 (Richard: "go ahead and close 51, 52, 53 then run 54") | [Plan](../../spec/EX-GC52_plan.md) (amendment A1 before the scored runs) | About 1.0 core-h on the cloud container (Linux, Python 3.13.18, numpy 2.5.3). Scoring and diagnostics: `src/exgc52_summary.py` → `summary.json`. This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**Yes. The 51 model has a preferred frame for moving matter, of the size the continuum analysis predicted.** The field of a moving mass has the right shape: it is Lorentz-contracted, with the clock part (g₀₀) the same as a boosted mass in GR. But two bodies moving side by side at V = 0.2 attract about 11–13% more strongly than relativity allows, after time dilation is accounted for. At V = 0.3 the excess is 28–33%.

The full dynamic run agrees with the semiclassical estimate (R_dyn = 1.115). The missing ingredient is gravity's "magnetic" part (g₀ᵢ). In GR that part cancels this excess, as the magnetic attraction between parallel currents offsets part of the electric force between moving charges. A single scalar field cannot carry it. In PPN terms the model has α₁ ≈ −8 and α₂ ≈ −1, while observations require |α₁| ≲ 10⁻⁴ and |α₂| ≲ 10⁻⁷.

So the 51 model is a mechanism for static gravity only. Moving matter needs a second ingredient, one that carries a direction (a vector or shift). This is in line with the PI's suggestion of torsion in the lattice or clock, or a non-symmetric rewriting rule.

| | What was scored | Result | Verdict |
|---|---|---|---|
| P1 | Q(V) = lattice ⟨∂E/∂N⟩ to 1e-6 | 5.9 × 10⁻¹² | PASS |
| P2 | moving-field amplitude at X = 0 within 2% | 0.996, 0.992, 0.989 (V = 0.1, 0.2, 0.3) | PASS |
| P2 | contraction ratio / γ within 2% | 0.996, 0.991, 0.986 | PASS |
| P3 | R(V) equals the lattice prediction within 1% | R/pred = 0.992, 0.983, 0.971 | FAIL at V = 0.2 and 0.3 |
| P3 | R − 1 exceeds 50% of 4V² (preferred frame) | 0.67 and 0.79 of 4V² at V = 0.2 and 0.3 | PASS |
| P4 | dynamic R_dyn(0.2) equals semiclassical R(0.2) within 3% | 1.115 vs 1.108 (ratio 1.007) | PASS |
| P5 | probe at rest: transverse impulse / GR within 3% (V = 0.2) | 0.962 | FAIL |

## What each part shows

**52A, the charge of a moving packet.**
- Q = ∂H/∂N equals the dispersion's ⟨∂E/∂N⟩ to 6 × 10⁻¹² for r = 0 and r = 1.
- Q(V)/Q(0) is 1.041 at V = 0.2. The continuum value γ(1 + V²) is 1.061.
- The gap comes from the packets' own internal momentum (σm = 3). That spread already gives Q(0) = 1.22 m, which is why A1 moved the continuum clause to reported-only.

**52B, the field of a moving mass.**
- The full nonlinear rule was solved in its steady comoving form at n = 129.
- The field is the Lorentz-contracted point field with charge Q(V): ψ·ρ = −κQγ/4π transversely, and ψ·|X| = −κQ/4π along the motion.
- Along z the amplitude is 1.0002 at every V. Along y it is 0.996–0.989, which passes P2 but is not 1.
- **Post hoc:** the reason is that a moving spin-up packet's charge centre sits off its axis, by −0.178 sites in y at V = 0.2 and −0.259 at V = 0.3. The shift reverses when the spin is flipped (checked directly). It is consistent in size with the relativistic centre-of-mass (Møller) shift of a spinning body, s·V/E ≈ 0.15, but I have not verified that identification.
- Averaging the +y and −y sides gives amplitude 1.0002–1.0004 at every V.
- So g₀₀ of a moving mass is exactly that of a boosted static mass. There is no preferred frame in the clock part.

**52C, forces from those fields.**
- **Co-moving probe.** A probe moving with the source at transverse distance 16 has R(V) = γṗ_y(V)/ṗ_y(0) = 1.022, 1.108 and 1.283. Lorentz invariance requires 1.
  - The lattice prediction (γ²·Q(V)/Q(0)·⟨E_N⟩(V)/⟨E_N⟩(0)) is 1.030, 1.128 and 1.322. The continuum point-particle value is γ⁴(1 + V²)²: 1.041, 1.174 and 1.435.
  - The 1% match fails because the probe at +y sits slightly farther from the shifted charge centre.
  - Post hoc, averaging the field gradient over ±y gives R_sym = 1.034, 1.133 and 1.326. These match the prediction to 0.4%, 0.4% and 0.3%.
- **Probe at rest.** The transverse impulse as the source passes is 0.983, 0.962 and 0.936 of GR's linear value, so P5 fails at V = 0.2. The deficit is the same two lattice effects: Q(V) below γ(1 + V²) because of the packet's spread, and the off-axis charge centre. It is not a model-versus-GR difference. At linear order the model's g₀₀ equals GR's, so the two give the same transverse force on a probe at rest.

**52D, full dynamics.**
- Two packets 14 sites apart move together at V = 0.2. The field is fully dynamic, starting from its steady comoving state.
- Packet 1's transverse force, with the single-packet self-force subtracted, gives R_dyn = γF(0.2)/F(0) = 1.115.
- That is within 0.7% of the semiclassical R(0.2) (P4 passes), and within 1.5% of R_sym (post hoc). Both packets carry the same off-axis shift, so their separation is unchanged.
- **The dynamic lattice model has the preferred-frame force predicted from the continuum Lagrangian.**

## Checks

| | | |
|---|---|---|
| C1 | moving solver at V = 0 vs the 51 static solver | 8.9 × 10⁻¹⁶ — PASS |
| C2 | linear part vs an independent conjugate-gradient solve, to 1e-8 | 1.7 × 10⁻⁷ — FAIL |
| C3 | 52D force F₁(0.2) at dt/2 within 1% | −3.7% — FAIL as run |

**C2 post hoc.** The difference grows as 1/κ (1.8 × 10⁻⁶ at κ = 10⁻⁷ and 1.8 × 10⁻⁵ at 10⁻⁸). That is a floating-point floor from forming f − f(1) with f ≈ −0.5, not a solver error. The check was badly designed.

**C3 post hoc.** My check script set κ from the V = 0.2 packet's charge, while the main run took it from the V = 0 packet's, so κ was 4.0% smaller. Per unit κ, the dt/2 force differs by 0.25%.

## Deviations

1. The P3 match and P5 fail because of the moving packet's off-axis charge centre, which I did not foresee. The registered scores stand. The symmetrised values are post hoc.
2. C2's registered tolerance sat below the floating-point floor of its own construction.
3. C3 was run with the wrong κ (a bug in the check script). The per-κ comparison is post hoc.
4. A1 (before the scored runs) moved the continuum charge clause to reported-only, after a small smoke test showed the packet-spread effect.
5. The field slices used by the summary (`partial_fields_B.npy`, 0.5 MB) are kept out of git, with a SHA-256 checksum committed. 52B/C regenerates them in about 2 minutes.

## What it means

The test did what it was meant to do, and the result is the expected negative one.
- **Already right:** the clock-rate field gets the g₀₀ of a moving mass exactly right. Event counting plus slowness transport already makes the clock part relativistic.
- **Missing:** a moving body must also drag something: a direction-carrying part of the substrate that responds to the mass current, of size g₀ᵢ = −4Vψ in the PPN gauge.
- **The force it must cancel:** the co-moving force excess R − 1 ≈ 4V² (point limit).
- **What the next ingredient must do:** it has to be sourced by the momentum density, ∝ V·(source). The coincidence rule as declared is symmetric in its two sites and cannot tell which way a firing pattern moves, so a handed or non-symmetric coincidence rule is a natural place to look.

**Suggested next (not approved):** declare a minimal vector/shift addition to the 51 rule. One option is coincidences that carry a direction, a lattice torsion, sourced by the matter current. Then rerun 52D unchanged to check whether R returns to 1. 52D is a ready-made test: 0.4 core-h at n = 65.
