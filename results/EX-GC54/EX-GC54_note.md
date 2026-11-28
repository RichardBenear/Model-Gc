# EX-GC54 — where the link shift comes from: Doppler-squared counting of a pattern's ticks gives frame dragging of about the right strength

2026-11-27 | EXPLORATORY | Claude | PI CLOSED 2026-11-28 (Richard: "close 54 and write up 55 and run it") | [Plan](../../spec/EX-GC54_plan.md) (amendment A1 after the operator check, before every other run) | About 2.7 core-h on the cloud container (Linux, Python 3.13.18, numpy 2.5.3). The PI raised the cap from 2 to about 3 core-h to finish the λ = 2 point and the dt/2 check. Scoring: `src/exgc54_summary.py` → `summary.json`; separate recomputation: `src/exgc54_recompute.py` → `recompute.json`. This is our construction, not LaViolette's Model G and not Chain Theory.

## Conclusion

**One counting rule supplies the link shift with no strength put in, and keeps the co-moving pair Lorentz-invariant within 2%. It misses its registered number, though, and the reason is now visible.**

- **The rule.** Each direction of a link counts its own two-site coincidences. A moving pattern's ticks are Doppler-shifted twice, by (1 ∓ v)², and the asymmetric part sources the shift. The source is 4κ × the energy current ½{Ĥ, v̂}.
  - The 4 is 2 × 2. The first 2 is light's two-site factor (the same one that doubles light bending), and the second comes from the Doppler square.
  - In uniform N this gives strength 3.85 against 53's momentum law. Lorentz invariance in this lattice needs 3.91 (53).
- **The co-moving pair** at V = 0.2:

  | Rule | R_dyn (1 = no preferred frame) | Registered |
  |---|---|---|
  | Energy current (D2-E) | 1.019 | 1.002 ± 0.008: FAIL; within 2% of 1: PASS |
  | Charge current (D2-q) | 0.916 | 0.977 ± 0.008: FAIL; outside 2% of 1, as predicted |
  | 53's law at λ = 2 | 1.056 | 1.056 ± 0.006: PASS |

- **Why D2-E misses (post hoc).** The energy current is smaller inside the packets' own clock-rate wells. Measured in place at t = 0, the pair's source corresponds to strength 3.37, not 3.85. With that strength, the straight line through 52 and 53 predicts R = 1.016, against 1.019 measured. So most of the gap is the well, not the rule's leading-order strength. Whether a counting rule should count per local tick (which would remove this factor) is the next open question.
- **The charge-current reading is out.** Before the runs I found that its operator is exactly twice the energy current (A1). The Dirac mass term carries no current, and Q̂ = 2Ĥ − mβ̂. So it acts as strength 7.7 and overshoots. With its in-place strength (6.96), the straight line gives 0.910, against 0.916 measured.
- **Linearity holds.** 52 (λ = 0), 54C (λ = 2) and 53 (λ = 4) give R = 1.1152, 1.0562 and 0.9975, with equal steps (0.0590 and 0.0587). So 53's λ* = 3.91 stands as a real interpolation (Codex review point).
- **The no-go holds.** A directed offset built from the N history alone (time-ordered coincidences) has a loop sum (curl) at most 1.2 × 10⁻⁴ of the needed shift's. It is a gradient, so the source has to read matter directly.

**What this does not show.** The rule is still a source operator added beside the 51 coincidence rule. It is not yet the coincidence rule itself with directed weights. "Count ticks by energy, not by N charge" is a choice that the outcome favoured, not one derived in advance. And the well factor means the rule as written is not exactly GR-compatible beyond leading order in these packets.

| | What was scored | Result | Verdict |
|---|---|---|---|
| P1 | N-history offset: curl < 2% of the shift's curl, 8 ≤ ρ ≤ 24 | max 1.2 × 10⁻⁴ | PASS |
| P2 | Σs = ⟨Ô⟩ to 1e-10; Σs = k-space value within 1% | 6 × 10⁻¹⁵; D2-E 1.5 × 10⁻¹⁰, D2-q 0.64 | PASS; D2-E PASS, D2-q FAIL |
| P3 | D2-q: R = 0.977 ± 0.008; Lorentz clause fails | 0.916; fails | FAIL (number); as predicted (clause) |
| P4 | D2-E: R = 1.002 ± 0.008; \|R − 1\| < 2% | 1.019 | FAIL; PASS |
| P5 | λ = 2: R = 1.056 ± 0.006 | 1.0562 | PASS |
| P6 | smooth gauge < 10⁻³ of (a) for Λ ≥ 32; linear in A within 10% | 0.011, 0.063, 0.097; ×10.00 and ×10.00 | FAIL; PASS |

## Details

**54B dynamics.** As 53E, with n = 65, separation 14, T = 16 and dt = 0.1, each pair with its single-packet subtraction.
- F₁(0) = 4.697 × 10⁻⁴ (D2-E) and 4.710 × 10⁻⁴ (D2-q). 53 had 4.699 × 10⁻⁴.
- F₁(0.2) = 4.690 × 10⁻⁴ (D2-E) and 4.228 × 10⁻⁴ (D2-q).
- The recompute (independent lstsq) agrees to 10⁻¹⁵. Fitting only t ≥ 2 gives 1.0184 (D2-E) and 0.9164 (D2-q).

**In-place source totals (post hoc, from the run records).** Four times the source over 53's P total, at t = 0:

| | pair | single packet | uniform N |
|---|---|---|---|
| D2-E | 3.37 | 3.51 | 3.85 |
| D2-q | 6.96 | 7.19 | 7.70 |

The energy current depends on the local N (the hops are weighted N_xN_y and the velocity carries the same weights). P does not.

**54D smooth gauge (post hoc diagnosis).** The registered fail is a clock-site effect, not a propagation effect.
- With the same pattern shifted so that β ≈ 0 at the launch and arrival site (cos phase), the delay is 2 × 10⁻¹³ to 9 × 10⁻¹³ at every Λ (`posthoc_gauge_clock.json`).
- With β largest there (sin phase, the registered case), the delay is −20 to −31 times β at the clock.
- So the pulse travels gauge-invariantly, and the one-clock reading picks up the shift at its own site. This is likely also the source of 53A's random-gauge residual. A clock readout should subtract the local β at the clock, or be placed where β = 0.

**Ring (reported, C3).** For the 53C spin-balanced ring the energy-current source has J = 1.907 against P's 1.976 (0.965), and the charge current has exactly twice that. Both vanish for ℓ = 0 (about 10⁻²⁰) and flip exactly with ℓ.

**Operators (P2, A1).** D2-E matches the k-space E·v to 1.5 × 10⁻¹⁰. D2-q is 2.000 × D2-E, because ½{Q̂, v̂} = 2 · ½{Ĥ, v̂} − (m/2){β̂, v̂} and {β̂, α̂} = 0. The plan's table had used the band-diagonal q·v, which leaves out interband terms.

## Checks

| | | |
|---|---|---|
| C1 | 53's law through the 54 code vs 53E's P₁(t), t ≤ 8, to 1e-10 | identical (0) — PASS |
| C2 | D2-E at dt/2 vs dt, P₁(t) for 2 ≤ t ≤ 8, within 1e-5 | 3.3 × 10⁻⁶ — PASS (holding the source over each step is not the cause of the gap) |
| C3 | rule sources zero for ℓ = 0 (< 1e-15) and flip with ℓ | 2.6 × 10⁻²⁰; 2 × 10⁻¹⁶ — PASS |

## Deviations

1. A1 (after the operator check, before every other run). The D2-q operator turned out to be twice the energy current, so P3's registered number was computed from the wrong strength. It is scored as registered.
2. Compute was about 2.7 core-h, against the 2 requested; the PI raised the cap to 3. The two runs paused at the cap (λ = 2 and dt/2) carry about 1960 s of paused wall time in their records, which is subtracted in `summary.json`.
3. The in-place strengths and the 54D phase diagnosis are post hoc.
4. A sign bug in my summary script (dividing by γ instead of multiplying) was fixed before scoring. It would have given R values about 4% low.

## What it means

- **Where frame dragging could come from.** Matter's ticks, counted by two-site directed coincidences, are Doppler-shifted twice. That gives the shift its GR strength without being told it, as long as a tick counts by its energy (rate), not by its N charge. Light's two-site nature fixes both factors of 2, so the static (51) and moving (53/54) parts of gravity come from the same counting.
- **Still to settle:**
  - **Local ticks or lab-time ticks.** The well factor says counting in lab time over-weights the redshift. Counting per local tick may be the fix, and the 54B setup is the ready test.
  - **Building it into the rule.** The directed weights should be part of the coincidence rule itself, not a separate source.
  - **The clock readout** should be corrected for the local β at the clock.
- **Suggested next (not approved):** rerun 54B with the per-local-tick energy current (one extra N factor, registered before running), about 0.9 core-h. If R → 1 within 0.5%, write the directed-weight coincidence rule that produces it.
