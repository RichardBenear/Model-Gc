# EX-GC34J — the re-slicing-exact kinetic term is stable and gives two clean waves; its first-order gravitomagnetism error is forced by the identity (EXPLORATORY)

**PI CLOSED 2026-11-23** (Richard: "You can close them", re 34I/34J, after ChatGPT's review `review/EX-GC34I-34J_chatgpt_review.md`). The review narrows one claim: that the first-order E4 error is fixed for all exact solutions is strongly suggested, not established (only three kernels and a max over nine cases were compared). 34K tests it.

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "you can close 34I. Go ahead with the suggested next step."
- The plan `spec/EX-GC34J_plan.md` and the `explore/LOG.md` entry were pushed before any run (b77c1d1).
- Code: `src/exgc34j.py`.
- Registered outputs:
  - `derive` → `derive.json`, `kernels.npz`;
  - `spec` → `candidates.json`;
  - `e4` → `e4.json`;
  - `sync` → `sync.json`;
  - `check` → `check.json`.
- Post hoc output: `sync_posthoc` → `sync_posthoc.json`.

This is our construction, not LaViolette's Model G.

## Conclusion

**No candidate passes as registered.** All three fail the gravitomagnetism benchmark E4 and its order test. All three
also failed the time-domain test G4′, but only because of one seed set the plan should not have used (see post hoc).
Beyond that, the forward-gradient, re-slicing-exact kinetic term (J1 and J2) does everything else asked of it:
- no growing root over the whole 34F/34G case set;
- exactly two isotropic waves within 7.6e-4 of the light speed at kh ≤ 0.05;
- the 7th-edge mode at the cutoff (ωh ≈ 5–6);
- static chains equal to the exact symmetries;
- real-space agreement to 5e-16;
- post hoc, a clean synchronous-gauge time evolution.

**Main new finding: the first-order E4 error cannot be removed by choosing among exact solutions.**
- J1, J2 and J3 give E4 errors equal to 5 digits: 3.1e-3, 6.3e-3, 1.25e-2 and 2.5e-2 at kh = 0.0125, 0.025, 0.05 and
  0.1. That is a ratio of 2.00 per doubling, i.e. exactly first order, about kh/4.
- J1, J2 and J3 are spread across solution spaces of dimension 123 (range 1) and 1263 (range 2), and they differ
  from each other by 44–67% of T.
- So the error is fixed by the identity itself, i.e. by the forward re-slicing gradient, not by the leftover freedom.
- The pre-registration probe showed it is not a source or readout placement either.
- My reading (not proven): displacement gauge puts the shift at vertices, while a forward re-slicing gradient puts it
  half a cell along each axis. No single shift location makes both gauges read as pure gauge beyond first order. A
  centred gradient would fix that, but in 34I it had no local solution at range ≤ 2.
- Matching T's first-order Taylor term (J2, J3) is impossible: the residual is 74% at both ranges.

## Registered results

| criterion | J1 (R=1, closest to T) | J2 (R=1, first-order match) | J3 (R=2, first-order match) |
|---|---|---|---|
| I0 identity exact | PASS (1e-14) | PASS | PASS (4e-14) |
| first-order match residual (reported) | 0.80 | 0.74 | 0.74 |
| G1 no growing root | PASS | PASS | **FAIL** (2 roots, 110 at kh = 3 and the zone boundary) |
| G2 two waves at c | PASS (max 7.5e-4) | PASS (7.6e-4) | PASS (7.1e-4) |
| G3 chains = exact symmetries | PASS | PASS | PASS |
| G5 static, gauge | PASS | PASS | PASS |
| G5 E4 ≤ 1e-2 at kh = 0.05 | **FAIL** (1.25e-2) | **FAIL** (1.25e-2) | **FAIL** (1.25e-2) |
| E4 order ratio ≥ 3 | **FAIL** (2.00) | **FAIL** (2.00) | **FAIL** (2.00) |
| G4′ synchronous time domain | **FAIL** | **FAIL** | **FAIL** |
| G7 real space | PASS (6e-16) | PASS | PASS |

**G4′ in detail** (all three candidates, both directions, both steps):
- The wave's tidal amplitude changed by ≤ 2e-10.
- Non-zero seeds changed the tidal field by 3–5e-8, against an allowed 1e-6, and the result did not depend on the
  step.
- It failed on two counts:
  - the zero-cluster seeds changed the tidal field by 1.5e-6 to 1.2e-5 (limit 1e-6);
  - the constraint residuals of the seeded runs reached 2e-6 to 1.2e-5 (limit 1e-9).
- Both are step-independent.

**My stated expectations:**
- J1 fails only E4: wrong, it also failed G4′, for the seeding reason below.
- I was unsure whether matching was possible: it is not, at range 1 or 2.
- That matching would give second order: moot.
- G4′ passes: failed as registered, clean post hoc.

## Post hoc (labelled)

**Constraint-filtered seeds** (`sync_posthoc.json`). The zero cluster of the synchronous system is 8-dimensional. Only
4 of those directions satisfy the lapse and momentum constraints. The other 4 are static non-solutions, not states of
the theory. My plan seeded all 8. With seeds restricted to the constraint-satisfying subspace:

| measure | result |
|---|---|
| constraint residual, every run | ≤ 8e-12 |
| non-zero seeds' effect on the tidal field | ≤ 4.6e-8, unchanged with the step |
| zero-sector seeds' effect on the tidal field | ≤ 7e-13 |

So the residual gauge carries no tidal field. All non-zero modes already satisfied the constraints. **The time domain
is clean.**

## What this means

- **What is achieved.** Exact re-slicing on the lattice, built in through one linear identity, gives a linear lattice
  graviton:
  - stable;
  - exactly two modes;
  - isotropic to 7.6e-4 at kh ≤ 0.05;
  - a non-metric edge mode pushed to the cutoff;
  - constraints preserved in time.

  That is what 34F–34H were missing.
- **What is not achieved.** Coupling to a moving source (gravitomagnetism) converges to GR only at first order, kh/4.
  At kh = 0.05 that is just over the 1e-2 benchmark, and it would pass at kh ≤ 0.04. The cause is the one-sided
  re-slicing gradient, which the identity forces. It is not the choice among solutions.
- **Options for a later study (none approved):**
  - (a) Accept first order and use finer lattices.
  - (b) A staggered formulation in which the shift lives on half-cells for both gauges. That changes the strain map,
    so the displacement and re-slicing gauges agree on where β lives.
  - (c) Search for a centred gradient at range 3 or more, or with a non-polynomial (rational) symbol, to see whether a
    local centred solution exists.
- J3's growth at very short wavelength (kh = 3) shows that "closest to T" at range 2 is not automatically safe. J1 is
  the cleanest candidate.

## Compute

About 2.8 CPU-minutes:
- derive: 74 s;
- spec: 54 s;
- E4: 2 s;
- sync: 18 s;
- post hoc: 11 s;
- check: under 1 s.

Platform: Linux cloud, Python 3.13, numpy 2.5, scipy 1.18.

## Deviations

- None in the registered runs.
- G4′ seeded constraint-violating zero-cluster directions, as the plan said. I treat that as a plan flaw and report
  it as a FAIL as registered. The filtered rerun is post hoc.
