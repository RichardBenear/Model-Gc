# EX-GC53 — a link shift for moving matter: clocks read it, matter sources it, and it removes 52's preferred frame

2026-11-27 | EXPLORATORY | Claude | Closure pending | [Plan](../../spec/EX-GC53_plan.md) (amendment A1 after 53A/B/C0, before 53C–E and the checks) | About 1.6 core-h on the cloud container (Linux, Python 3.13.18, numpy 2.5.3). Scoring: `src/exgc53_summary.py` → `summary.json`. This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**Adding a signed timing offset on each link (GR's shift, g₀ᵢ), with GR's strength, removes the preferred frame that 52 found.**
- **Co-moving pair.** In 52 the excess attraction at V = 0.2 was 11.5%. With the shift the pair's sideways pull matches relativity:

  | V | 0.1 | 0.2 | 0.3 |
  |---|---|---|---|
  | R_dyn | 1.000 | 0.997 | 0.989 |

  At V = 0.2 the spin-up and spin-down runs agree to 0.14%.
- **The strength relativity needs.** Comparing 52 (no shift) with 53 (shift at GR's strength) gives λ* = 3.91. GR's value is 4. So the 51 clock-rate field is GR-compatible for moving sources to this accuracy, and the shift was the only missing piece.
- **Rotating ring.** A one-clock light loop around a rotating ring of Dirac matter shows a timing difference that:
  - flips exactly when the rotation flips;
  - is zero (about 10⁻²⁰) for the same ring not rotating;
  - does not depend on how the ring is tilted against the lattice (0.03%).
- **Size of the ring signal (registered fail).** It approaches GR's frame-dragging dipole as the loop grows (11%, 6% and 4% above at half-widths 20, 28 and 36), but misses the registered 5% except for the largest loop. The dipole-pattern tests also fail as registered. Post hoc, the cause is that my ring (radius 10) is too large for d ≤ 40 to be the far field. Against the exact field of the actual current distribution, the lattice loops agree to 1.3–4.8%, except the nearest in-plane loop.

**What this does not show.** The shift's strength was put in from GR (λ = 4). I tried to derive it from the order of firings and could not, as the plan says. So 53 is a calibration that shows the model's other pieces are consistent with GR, not a mechanism for the shift. It answers the PI's question of where something must be added: on the links, as a directed offset that is sourced by the matter current. It does not yet say why the substrate does that.

| | What was scored | Result | Verdict |
|---|---|---|---|
| P1a | loop delay = 2Σβ/N⁴ within 2% (N ≡ 1) | 0.983 | PASS |
| P1a′ | same, weak-field N | 0.826 | FAIL |
| P1b | pure-gauge pattern gives < 10⁻³ of (a) | 7.7 × 10⁻³ | FAIL |
| P1c | one side gives ¼ of (a) within 2% | 0.920 | FAIL |
| P2 | ΣP = ⟨sin k⟩ to 1e-10; ∂H/∂β = −P to 1e-8 | 2.6 × 10⁻¹³; 1 × 10⁻¹⁶ | PASS |
| P3a | ℓ = 0 gives 0; ±ℓ flip; N identical for ±ℓ | 5 × 10⁻²⁰; 8 × 10⁻¹⁶; 9 × 10⁻¹⁶ | PASS |
| P3b | small-loop log-slope −3 ± 0.15 | −3.44 (plane), −2.54 (axis) | FAIL |
| P3c | axis/plane ratio −2 ± 0.1 | −0.72, −1.14, −1.46, −1.63 | FAIL |
| P3d | lattice / analytic dipole within 5% (d ≥ 24 and enclosing) | 0.83–1.31 | FAIL |
| P3e | tilted ratio = aligned ratio within 3% | 0.9997 | PASS |
| P4 | R_dyn = 1 within 2% at V = 0.1, 0.2, 0.3; spin up vs down within 1% | 0.9998, 0.9975, 0.9894; 0.14% | PASS |
| P5 | λ* = 4 ± 0.3 | 3.91 | PASS |

## Details

**53A, the one-clock readout.** These were the first parts run, with the final code, and the registered fails stand. The post-hoc diagnosis was written into A1 before the other parts ran.
- **Weak-field N miss.** This comes from the clock's centroid density. u² + w² is not the conserved energy when N varies. With the conserved energy density, the weak-N case matches the uniform one (0.938 against 0.936).
- **Remaining 6%.** This is lattice dispersion at k₀ = 0.5. The lattice eikonal factor is cos k₀ / cos²(k₀/2) = 0.935, and the wave matches it to 0.1%.
- **Gauge and one-side patterns.** Their misses are of the same few-percent lattice size. The random gauge pattern also scatters the pulse.
- **What the clock reads.** In short, it reads ∮β through a known lattice factor.

**53B, the coupling.** Matter couples to the shift through the link momentum, which is exact to rounding. A spin-polarized packet at rest carries a circulating spin current (largest link value 2.4 × 10⁻⁶ against 1.6 × 10⁻⁴ for the V = 0.2 packet) and zero net momentum. That current is a small spin–spin source of the shift, and a natural lead into the later torsion study.

**53C/D, the rotating ring.**
- **Source.** J = ±1.976 (lattice). q is identical for ±ℓ.
- **Enclosing loops.** Readout over dipole is 1.115, 1.060 and 1.038. The wave delay with the energy centroid is 0.944–0.955 of the readout, against the lattice factor 0.935.
- **Small loops.** The ratios converge on the dipole roughly as 1 + c(a/d)², with c ≈ 1.7 in plane and c ≈ −1.1 on the axis.
- **Post hoc, against the exact field.** I summed the exact free-space shift of the actual lattice current (direct Biot–Savart sum, `posthoc_freespace.json`). Lattice over exact is:
  - 1.026, 1.018 and 1.013 for the enclosing loops;
  - 1.016–1.030 on the axis;
  - 1.03–1.05 in the plane at d ≥ 24.

  The leftover 1–5% is probably the finite Dirichlet box, whose faces carry only the dipole term. I have not checked that.
- **Conclusion.** The scored P3b–d failures measure my choice of ring size, not the physics.

**53E, full dynamics.** The pair runs with N, the three shift components and both packets all dynamic. Results:
- F₁(0) = 4.699 × 10⁻⁴ (52: 4.685 × 10⁻⁴; the small difference is the spin-current shift).
- F₁(0.1) = 4.674 × 10⁻⁴, F₁(0.2) = 4.592 × 10⁻⁴, F₁(0.3) = 4.435 × 10⁻⁴.

So the lab-frame force now falls as 1/γ, as relativity requires. In 52 it rose.

## Checks

| | | |
|---|---|---|
| C0 | steady shift of a moving source vs the boosted GR value at ρ = 8–24, within 3% | mean 0.977, range 0.956–0.990 — FAIL (sign and size right; misses by up to 4.4% at some ρ) |
| C1 | ring shift solve vs conjugate gradient | 1.4 × 10⁻¹⁴ — PASS |
| C2 | 53E V = 0.2 at dt/2, P₁(t) for t ≥ 2 | 7 × 10⁻⁷ — PASS |
| C3 | shift off reproduces 52's P₁(t) | exactly 0 — PASS |

## Deviations

1. 53A ran first, as the scored run. A1 changed how 53D is scored (readout sums instead of waves for the small loops), using what 53A showed. 53D had not run at the time.
2. The ring size (radius 10) was too large for the registered far-field window, so P3b–d fail.
3. C0 is pointwise outside 3%. The likely cause is the compact source (σ = 2) and the box at n = 65, but this is not verified.
4. Compute was about 1.6 core-h, against the 2 requested.

## What it means

- **Structure.**
  - In the static limit (51), a clock-rate field carried by coincidences, with light's slowness as the transported quantity, gives Newton, γ = 1 and β = 1.
  - A directed link offset sourced by the matter current gives the moving-matter corrections (52 → 53).
  - Together they behave like linearized GR for slow matter, on a lattice, with one clock reading frame dragging.
- **What is still put in by hand:** the shift's law and its strength 4.
- **Suggested next (not approved):**
  1. A mechanism for the shift. For example, a non-symmetric coincidence rule where a firing at x followed by one at y within the coincidence window is counted differently from y then x. Test whether its natural strength is 4. 53E is the ready test, at 0.4 core-h per speed.
  2. A cleaner far-field ring test (a smaller ring, or a larger box with multipole boundaries).
  3. Later: torsion from the spin current found in 53B, and directional strain for tensor waves.
