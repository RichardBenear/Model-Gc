# EX-GC34O — on Y1, a scalar's own sources now conserve energy and momentum exactly (EXPLORATORY)

**PI CLOSED 2026-11-23** (Richard: "I sign off on 34O and 34P").

Claude, 2026-11-24. PI go-ahead (Richard, 2026-11-23): "Yes, go ahead with 34O".
- Source: `review/EX-GC34N_chatgpt_review.md`.
- The plan `spec/EX-GC34O_plan.md` and the `explore/LOG.md` entry were pushed before any run (ffdbdb6).
- Code: `src/exgc34o.py` (`run`, then `recompute`). Output: `results.json`.

This is our construction, not LaViolette's Model G.

## Conclusion

**PASS on every registered gate.** With one change to how the scalar couples to the spatial metric, the packet's
sources satisfy both of Y1's lattice conservation identities exactly, to rounding, at every spacing, time and
direction. 34N's leftover second-order momentum error is gone.

**The change.** On Y1 the scalar's normal stress along a direction now uses the product of the two gradients next to a
site, G_a(v) G_a(v − e_a), instead of the average of their squares. In the action this is one term,
¼ h_aa(v) G_a(v) G_a(v − e_a) at each site, in place of the parallel part of 34N's link coupling. Everything else is
34N's: the vacuum, the flat scalar evolution, the lapse and shift couplings, the transverse terms and the shear term.
The continuum limit is unchanged.

**Why it works.** The time derivative of the link momentum splits into exact lattice differences, using the exact
lattice identity that the two mixed differences of φ around a plaquette are equal. Only the parallel normal term
needed changing. That confirms ChatGPT's 1D finding and shows the transverse and shear terms were already right.

| | 34N stress | 34O stress |
|---|---|---|
| energy identity | exact | exact (≤ 8.8e-15) |
| momentum identity, packets, h = 0.25 | 0.18–0.37% | ≤ 2.7e-15 |
| momentum identity, random 3D fields | 31–41% | ≤ 2.8e-16 |

## Registered results

| gate | result |
|---|---|
| O0, 1D control: product form | 1.9e-16 / 3.3e-16 (h = 1 / 0.5) |
| O0, 1D control: 34N form (null) | 1.2 / 1.4 |
| **O0** | **PASS** |
| **O1**, derivation check (Richardson) | 1.1e-11 / 8.0e-12 → **PASS** |
| **O2**, energy identity, 18 packet cases + random fields | ≤ 3.9e-14 → **PASS** |
| **O3**, momentum identity, 18 packet cases + random fields | ≤ 2.7e-15 → **EXACT** |
| **null control**: 34N stress in this harness vs its registered r_P | identical to 4 digits in all 18 cases → **PASS** |

**Overall: PASS.**

**My stated expectations:** everything passes, with O3 at rounding. That was right.

**Reported, not gated:**
- The integrated normal stress of the two forms differs by 18% / 4.3% / 1.1% (axis) and 15% / 3.5% / 0.9% (generic)
  at h = 1 / 0.5 / 0.25. That is second order: both have the same continuum limit.
- **Cutoff sign.** For the zone-boundary mode the product G_a(v) G_a(v − e_a) is −G² (the average of squares is +G²).
  So the parallel gradient stress changes sign at the cutoff. That is harmless for a linear source on smooth fields,
  but it would matter for a nonlinear or positivity-based use.

## Recomputation and shared dependencies

`recompute` is a separate path:
- d/dt of the link momentum is taken directly from φ̈ (no polarisation trick, no `midentities`);
- the stress divergence is written out from the plan's formulas.

Result on the h = 0.5 packets (t = 6) and on random fields: ≤ 1.9e-15.

- **Shared:** `exgc34n.lap` and `exgc34m.packet`.
- **Slip:** my first recompute had the shear term's sign wrong, giving O(1). I fixed the sign in the recompute path
  only. The registered run was not touched.

## What this means

- **Y1 now has a matter field whose own sources are fully consistent with both gauge symmetries.** The energy and
  momentum identities hold exactly. That is the condition ChatGPT set before any coupled run.
- **Next (not approved):** review step 2 on Y1. Couple this packet's sources to the Y1 vacuum and test whether the
  packet's energy and momentum change gravity as linearised GR predicts. That is the "missing link" test 34M set out
  to reach.

## Compute

About 23 CPU-seconds, plus about 1 s for the recompute.

Platform: Linux cloud, Python 3.13, numpy 2.5.

## Deviations

None in the registered run. The recompute-path sign slip is described above.
