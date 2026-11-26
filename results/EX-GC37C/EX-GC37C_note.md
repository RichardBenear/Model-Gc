# EX-GC37C — does the event-counting mechanism survive a moving mass?

2026-11-26 | EXPLORATORY | Claude | [Plan](../../spec/EX-GC37C_plan.md), with amendment A1 (before any run). PI closure pending. Cloud container, Python 3.13. Under 3 CPU-seconds. This is our construction, not LaViolette's Model G.

**Only if the firing-rate field is a wave field.**
- **A diffusing, consumed carrier (Model Gc's picture) fails at once.** A moving sink leaves a wake. Its field is 2.6× weaker ahead than behind at v = 0.1, and the sink feels a backward pull from its own field, i.e. drag.
- **The same consumption feeding a wave field with speed 1 passes.**
  - The field is identical ahead and behind.
  - It is Lorentz-contracted by γ to within 1.3% up to v = 0.4.
  - The self-gradient is exactly zero, so there is no drag.

So the 37B mechanism survives motion only if N obeys a wave equation. "Matter consumes a substance" has to become "matter is a source for a relativistic field."

## Results

Comoving steady state on a 160³ Y1 lattice. The sink has core 1.5 and GM = 0.16 at rest. The readout distance is d = 10 cells.

| v | Diffusing: ahead/behind | Diffusing: side/ahead | Diffusing: self-gradient | Wave: ahead/behind | Wave: side/ahead | γ | Wave: self-gradient |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 0 | 1 | 1 | 1 | 0 |
| 0.05 | 0.607 | 1.283 | 9.2e-4 | 1.000 | 1.0014 | 1.0013 | 0 |
| 0.1 | **0.387** | 1.596 | 1.7e-3 | **1.000** | 1.0058 | 1.0050 | 0 |
| 0.2 | 0.188 | 2.159 | 2.8e-3 | 1.000 | 1.0237 | 1.0206 | 0 |
| 0.4 | 0.108 | 1.995 | 4.1e-3 | 1.000 | 1.1054 | 1.0911 | 0 |

The self-gradient is positive, so N rises ahead of the sink. A particle pulled toward lower N is therefore pulled backward.

**Far field, 40 cells out, relative to rest.**
- Diffusing at v = 0.1: 0.38 ahead and 2.7 behind. The 1/r field ahead is lost beyond about D/v and piles up in a wake behind.
- Wave at v = 0.1: 0.99 both ways.

## Predictions

- **P1 (diffusing has a preferred frame):** fore/aft is off by 61% at v = 0.1 (registered: more than 10%), and the self-gradient is nonzero and rising. **The "grows linearly" clause fails.** It rises sublinearly: 1.83× from v = 0.05 to 0.1, and 1.45× from 0.2 to 0.4, which is Oseen-type saturation. So P1 passes on the asymmetry and drag, and fails on linearity.
- **P2 (wave has none): PASS.**
  - Fore/aft is 1 to machine precision.
  - Side/axial matches γ within 0.07%, 0.3% and 1.3% at v = 0.1, 0.2 and 0.4 (gate 2%).
  - The self-gradient is 0 to roundoff.

**Separate check.** The continuum boosted solution of (∇² − v²∂_x²)c = s/D is c ∝ γ/√(γ²x² + ρ²). That predicts an unchanged deficit ahead and a deficit to the side larger by γ. The lattice gives both: ahead10 changes by 0.8% at v = 0.4, and side/ahead tracks γ. The v = 0 rows reproduce the 37B B2 field.

## What it means (inference)

- **EX-GC25's objection is confirmed quantitatively.** A consumed, diffusing carrier defines a rest frame. Model Gc's literal picture, matter eating a substance, cannot be the source of N for moving masses.
- **The mechanism survives with one change.** Matter must source N through a wave equation with the same speed as light. Then N has no rest frame, and the 37B counting gives the weak-field metric for moving masses too (at first order in this test).
- **That sharpens the open question to: why does the firing-rate deficit propagate as a wave at the speed of light?** One self-consistent answer, untested: a change in a site's firing rate is passed to its neighbours by the same two-site coincidences that carry light. It would then travel at the local light speed, and gravity and light would share one cone, as 34A required.

## Files

[Code](../../src/exgc37c.py) · [Results](results.json)
