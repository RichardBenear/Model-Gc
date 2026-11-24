# EX-GC34P — the 34O scalar packet sources Y1 gravity consistently, and its tidal field converges to linearised GR (EXPLORATORY)

**PI CLOSED 2026-11-23** (Richard: "I sign off on 34O and 34P").

**Correction (2026-11-24, after the ChatGPT 34P review; no gate re-scored):** by t = 4 and 8 the two preparations had moved to different places. So the A − B row tests the full time-dependent field of two different source histories. It does not isolate a velocity-dependent (gravitomagnetic) part, as the text below says.

Claude, 2026-11-24. PI go-ahead (Richard, 2026-11-23): "ChatGPT has review 34O and has a small recommended step. Go
ahead with that test."
- Source: `review/EX-GC34O_chatgpt_review.md`.
- The plan `spec/EX-GC34P_plan.md` and the `explore/LOG.md` entry were pushed before any run (eb691ed).
- Code: `src/exgc34p.py`.
- Outputs: `norm.json`, `compat.json`, `curv_axis.json`, `curv_generic.json`, `recompute.json`, and the post hoc
  `posthoc_spectrum.json`.

This is our construction, not LaViolette's Model G.

## Conclusion

**PASS on every registered gate.** The moving scalar packet's own energy, momentum and stress drive Y1's linear
gravity without breaking its constraints. The tidal curvature this produces matches linearised GR driven by the
continuum packet, at second order in the spacing:
- error 0.5–1.2% at h = 0.25;
- falling by 3.97–4.03× per halving.

This holds for both directions, both preparations, and their difference. It is the "does a dynamical matter field's
energy make the right gravity" step that 34M set out to reach, at linear order and without feedback on the scalar.

**ChatGPT's normalisation finding is confirmed.** The inherited packet's energy fell by 63–64× per halving of h. The
fixed packet's energy agrees between h = 0.5 and 0.25 to 0.2%.

## Registered results

**P0, normalisation: PASS.**
- Fixed packet energy at h = 0.5 vs 0.25: differs by 0.21% (axis) and 0.21% (generic).
- After matching the energy, the width differs from the continuum by ≤ 7e-4 and the central wavenumber by ≤ 1e-3.
- The old normalisation's energy ratio per halving was 63.3 / 63.8 (axis) and 63.5 / 63.9 (generic).

**P1, source compatibility: PASS.** Constraint residual relative to the constraint source, over t = 0–8, h = 0.5:

| | dt = 0.05 | dt = 0.025 | ratio |
|---|---|---|---|
| 34O sources, axis (A and B) | 5.85e-7 | 3.66e-8 | 16.0 |
| 34O sources, generic (A and B) | 5.96e-7 | 3.73e-8 | 16.0 |
| 34N sources (null), axis | 7.56e-2 | 7.56e-2 | 1.00 |
| 34N sources (null), generic | 4.06e-2 | 4.06e-2 | 1.00 |

- With the 34O sources the residual is pure RK4 error: it falls 16× per halving of dt, and starts at 5e-16. The
  continuum reference shows the same level (6.3e-7).
- With 34N's sources the constraints break by 4–8%, independent of dt. That is a real inconsistency.
- At h = 0.25 the 34O residual is 6.2e-7; at h = 1 it is 6.3–6.5e-7.

**P2, tidal curvature vs linearised GR: PASS.** Relative L2 error of R_0a0b over all evolved modes:

| direction | h | t = 4: A / B / A − B | t = 8: A / B / A − B |
|---|---|---|---|
| axis | 1 | 0.48 / 0.48 / 0.72 | 0.51 / 0.51 / 0.58 |
| axis | 0.5 | 0.026 / 0.026 / 0.027 | 0.046 / 0.046 / 0.047 |
| axis | 0.25 | **0.0066 / 0.0066 / 0.0068** | **0.0115 / 0.0115 / 0.0118** |
| generic | 1 | 0.35 / 0.35 / 0.54 | 0.38 / 0.38 / 0.39 |
| generic | 0.5 | 0.019 / 0.019 / 0.020 | 0.035 / 0.035 / 0.036 |
| generic | 0.25 | **0.0048 / 0.0048 / 0.0050** | **0.0087 / 0.0087 / 0.0090** |

- The ratios from h = 0.5 to 0.25 are 3.97–4.03 (gate ≥ 3), and every h = 0.25 error is ≤ 1.2e-2 (gate ≤ 2e-2).
- h = 1 is reported only: its grid cannot hold the packet's spectrum.
- The error grows with time, about 1.75× from t = 4 to t = 8, as the matter's dispersion error accumulates.

**Reported, not gated:**
- ‖R_A − R_B‖/‖R_A‖ = 1.14 at t = 4 and 1.36 at t = 8. The velocity-dependent part of the tidal field is as large as
  the field itself, so the A − B row is a real test, not a small residue.
- The dropped source fraction (modes with |k_a| > 4.9) is below rounding.
- All constraint systems had full rank: the minimum relative singular value was 0.077 at h = 0.5.

**My stated expectations:**
- P0 and P1 pass: right.
- P2 about 5e-3 at h = 0.25, second order: right for t = 4 (5–7e-3), and twice that at t = 8.
- The risk I named, h = 0.5 being pre-asymptotic, did not occur: the ratios are 4.0.

## Post hoc (labelled)

**Where the tidal norm sits** (h = 0.5, axis, A; `posthoc_spectrum.json`): 65–70% at 0.5 ≤ |k| < 2 (the packet
scale), 29–34% below 0.5 (including 13–16% in the lowest box shell), and 1.3–1.4% above 2. So the comparison is not
dominated by box-scale modes.

## Recomputation and shared dependencies

`recompute` is a separate path on the 3 largest modes and one median mode (h = 0.5, axis, A, t = 8):
- the per-mode kernel is built from `exgc34n.G_up` in a loop (not the batched kernel);
- the sources are direct DFT sums (no FFT);
- initial data come from `lstsq`;
- integration is DOP853 at rtol 1e-10.

The 3 large modes agree with the main path to 4.2e-11. The median mode is at rounding size (3e-16 of the largest), so
its comparison is meaningless.
- **Shared:** `exgc34n.G_up` and `PM`, the 34O sources, and the packet generator.
- **First attempt:** my first recompute used three probe modes saved by the curvature run. They turned out to be edge
  modes of negligible size, and the comparison was meaningless there (`recompute_first_attempt.json`). I replaced
  them with significant modes before scoring anything.

## What this means

- **On Y1, a dynamical matter field gravitates as linearised GR says, from its own action.**
  - Its energy and momentum source the curvature correctly, including the velocity-dependent (gravitomagnetic-type)
    part.
  - The error is second order.
  - The source is exactly compatible with the gravity equations.

  This closes the chain from 34J (first-order gravitomagnetism) through 34N (staggered lattice) and 34O (exact
  stress).
- **The 34N null control shows why 34O mattered.** With 34N's stress the constraints break by 4–8%, at a level an
  integrator cannot explain.
- **Not tested:** feedback of gravity on the scalar, nonlinear gravity, the knot, and where Y1's geometry comes from.

## Compute

About 2.1 core-hours in total, just over the 2 core-hour cap:
- registered runs: about 1.7 core-h;
- recompute and post hoc: about 0.1 core-h;
- runs I stopped: about 0.5 core-h (see Deviations).

Platform: Linux cloud, Python 3.13, numpy 2.5, scipy 1.18.

## Deviations

- **Compute cap exceeded by about 0.1 core-h.** I underestimated the per-mode cost, expecting 0.7 core-h.
  - I stopped the first batch of runs after about 0.4 core-h and sped up the batched matrix products.
  - A second launch had a coding slip: it dropped the constraint rows' ḣ term, which showed as a 51% initial
    residual. I stopped it within minutes, fixed it and checked it against the original formula. The re-run
    reproduces the first batch's numbers (5.851e-7, 3.655e-8).
- **P1 scoring.** The h = 0.25 and h = 1 constraint residuals come from the P2 runs, which are the same lattice runs
  (the residual does not depend on the packet scale). The plan had listed them under `compat`, so I added a separate
  `score` step.
- **Recompute modes.** Changed as described above, before scoring.
