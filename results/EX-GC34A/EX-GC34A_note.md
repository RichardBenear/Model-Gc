# EX-GC34A — one shared link structure: do all fields share one light cone? (EXPLORATORY)

Claude, 2026-11-22/23. PI go-ahead (Richard, 2026-11-22): "Go ahead with the 34A plan." Plan `spec/EX-GC34A_plan.md`
(afa1f5f, committed and pushed before production). Code `src/exgc34a.py` (solver), `src/exgc34a_rays.py` (predictions),
`src/exgc34a_eval.py` (registered scoring → `summary.json`), `src/exgc34a_check.py` (registered separate recomputation
→ `check.json`), `src/exgc34a_diag_posthoc.py` (post hoc → `diagnostics_posthoc.json`). Our construction, not
LaViolette's Model G. Status: awaiting PI decision.
**PI CLOSED 2026-11-23** (Richard: "go ahead and close 34A to 34E").

## Conclusion

**By its registered readout, 34A did NOT pass: M1 23/54 and M2 0/15 cases met the criteria. The cause is a flaw in my
plan, not in the physics or the solver.** The registered readout was the position of the packet's energy centroid in a
periodic box, and I made the boxes too small for how much the packets spread. Their tails wrapped around the box and
dragged the centroid. I estimated the tails but forgot diffraction: by the end of M2 a packet is about three times
wider than at the start.

**The registered separate recomputation does not have this problem, and it shows what the study was built to show.**
It reads speed and direction from the final field's wave-number spectrum, which doesn't care where the box wraps.

- **One light cone for every field.** On flat and anisotropic backgrounds, in three directions, the massless, massive,
  two-component and nonlinear sigma-model fields all move at the predicted velocity. The match is to within 2e-5
  (linear fields) and 1.4e-4 (sigma model), and the direction matches to within 1e-5 rad.
- **The control fails, as it should.** The field with its own link weights is 4.9% fast against the shared prediction
  and matches its own prediction to 1e-5. This is the GC31A two-cone situation, and the readout catches it.
- **Light bending comes out in the textbook ratios on prescribed backgrounds.** Final propagation angle of massless
  waves:

| background | spectral (field) | rays | ratio to lapse-only (field) |
|---|---:|---:|---:|
| B2 lapse only (PPN γ = 0) | −0.01262 | −0.01238 | 1 |
| B3 GR-like (γ = 1) | −0.02535 | −0.02483 | 2.009 |
| B4 Nordström (conformal) | +0.00003 | 0.00000 | −0.002 |

  Massive fields bend more, as rays of their own mass predict. For example, B2 μ = 1 gives −0.01390, and B4 μ = 1
  gives −0.00133 where the light-like waves give 0.
- **All fields agree with each other under the registered readout too.** The broken centroid values are identical for
  F1 and the F3 symmetric mode, and agree with F4 to 3e-4. The wrap distorts every field in the same way.

The spectral numbers are the registered recomputation path, but the **pass/fail criteria were registered on the
centroid readout**, so I do not re-score 34A on them. A corrected rerun, with the readout fixed before any data, is
proposed below.

## Registered results (`summary.json`)

**M1 (uniform backgrounds, 80³ box, h = 0.2, T = 5).**
- Along x (box face direction), every field reads 1.8–3.7% slow, so all fail.
- Body diagonal: flat B0 passes (0.09–0.30%). On B1/B1′ the massless modes read 0.32–0.33% slow and just fail; the
  massive modes and F4 pass.
- Generic direction: B1/B1′ pass (0.13–0.41%, F4 within its 1%). Flat B0 fails (0.3–0.8%, angle up to 0.0044 rad).
- The pattern follows how close each direction brings the packet to a box face, not the field type.
- F5 control: fails the "≥ 3% off shared" test along x only, because the wrap cancels its true +4.9% excess there.
- Energy is conserved to ≤ 1.5e-5 (linear) and 1.5e-4 (sigma).

**M2 (bumps).** Every field reads α ≈ −0.05 in every background, including Nordström where the true value is 0.
Registered tolerance ≈ 6e-4: 0/15.

**Post hoc check of the cause (`diagnostics_posthoc.json`).** In the final fields, 23% (M1 x-run) and 45% (M2) of the
squared field lies within 2 length units of a periodic face. The ray bundles, which include the packets' spread in
wave number, show the expected spreading.

## Deviations (recorded, mine)

1. **Box too small for diffraction** (M2 diffraction length kσ² = 12 against T = 36). This is the cause of the failures.
2. **Figure not produced.** matplotlib is not installed in this container. Deferred to the rerun.
3. **Compute.** Runs 1.02 core-h, scoring and rays 0.28, recomputation 0.03, post hoc check 0.1: about 1.45 core-h of
   the 2 core-hour cap.

## Proposed: EX-GC34A-r (needs PI go-ahead; would exceed the 2 core-h cap by about 0.6)

Same fields, backgrounds and tolerances. The readout would be fixed before any run: propagation velocity from the
spectral group velocity, sampled at the registered early and late window times (wrap-immune). Packet width and M1 box
would be enlarged so the centroid becomes a valid secondary readout. About 1.1 core-h in the cloud.

The alternative is to accept 34A as a method validation (spectral readout) and move to 34B, which uses the same
solver. The PI decides.

## Platform

Linux cloud container (4 cores), Python 3.13.18, numpy 2.5.3, scipy 1.18.1. Final fields are in `local_raw/EX-GC34A/`
(not committed), with SHA-256 values in `manifest_sha256.json`.
