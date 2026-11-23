# EX-GC34B — what does each part of the shared structure respond to? (EXPLORATORY)

Claude, 2026-11-23. PI go-ahead (Richard, 2026-11-23): "go ahead with 34B". Plan `spec/EX-GC34B_plan.md` (7482f34,
pushed before the runs). Code `src/exgc34b.py`; record `sources.json`; regenerated knot records in `states_h*/`.
Our construction, not LaViolette's Model G.
**PI CLOSED 2026-11-23** (Richard: "go ahead and close 34A-r, 34B, 34C and proceed with 34D").

## Conclusion

**The clock (lapse) part of the shared structure responds to total energy, for every kind of matter tested. The
spatial part responds to stress, which totals zero for a body at rest.** So the shared structure gives a universal
coupling (universal fall, binding included). How much light bends is then decided by the gravity field's own
equations, not by the matter.

- **Lapse source = energy, every case** (relative difference ≤ 1e-6): B = 1 and B = 2 knots at rest, B = 1 at a finer
  lattice, and massless and massive packets before and after propagating. For B = 2 that includes its binding energy
  (−0.094, 7.8% of a single knot's energy).
- **Spatial source of a knot at rest = lattice offset only.** ∫T^{ij} is isotropic for B = 1 (each diagonal
  = (E4 − E2)/3, off-diagonal ≤ 3e-8). For B = 2 (torus) it is nearly isotropic, with the axis component 3.5% higher.
  The total is the lattice Derrick offset E4 − E2, 6.4% of E at h = 0.2 and 3.6% at h = 0.15. That ratio, 0.565,
  matches the h² scaling (0.5625), so the spatial monopole goes to zero in the continuum, as von Laue requires.
- **Radiation and moving matter do source the spatial part.** A massless packet gives ∫T^{xx}/E = 0.973 against the
  lattice expectation v_g/v_p = 0.970. A massive packet gives 0.877 against 0.870 (continuum v² = 0.9). Transverse
  components are ≤ 1.3% of E.

## Registered criteria

| item | result | verdict |
|---|---|---|
| lapse source = H | ≤ 1e-6 in all 7 cases | PASS |
| knot trace of ∫T = E4 − E2 | agrees to ≤ 1e-10 (closed form vs numerical derivative ≤ 7e-11) | PASS |
| knot components ≤ 1.2 × \|E4 − E2\| | max 0.34 × | PASS |
| offset ratio h = 0.2 → 0.15 in [0.40, 0.75] | 0.565 | PASS |
| massless packet ∫T^{xx}/E within 5% of v_g/v_p; transverse ≤ 5% | 0.36%; 1.3% | PASS |
| massive packet ∫T^{xx}/E within 5% of "v_g v_p" | 0.877 vs 0.941: −6.8% | **FAIL as registered** |

**The massive-packet row failed because of a slip in my plan, not in the physics.** The plan derives the expectation
as momentum flux = P v_g with P = E/v_p, which is E v_g/v_p. The table then wrote v_g v_p, and the code used that
too. Against the correctly derived v_g/v_p (0.870) the result is within 0.75%. I record the registered verdict as
FAIL and the slip as a plan error.

**Regenerated knots.** The B = 1, h = 0.2 state is bit-identical to the committed 26E1 state (same SHA-256). B = 2
and B = 1 at h = 0.15 differ in the last digits: energies agree to 3e-8 and 1.4e-6. The committed 26E2 state came from
the PI's PC with different library versions. Both are within the optimiser's stopping tolerance.

## What this means for 34C

On this structure, matter couples universally: everything falls alike, and binding energy weighs. That was 34B's
question, and the answer is yes, by construction. Light bending is a separate matter:

- If the lapse and the spatial weights each obey their own simple field equation, a body at rest produces only a
  lapse field (its stress totals zero). PPN γ = 0, and light bends half as much as GR (the 34A-r lapse-only case).
- GR's γ = 1 needs the spatial metric to be driven by the energy too, which is what the trace structure of the
  linearised Einstein (Fierz–Pauli) action does. Putting that action in by hand would impose GR, not derive it.

So 34C should ask what in Model Gc could select the gravity field's own action. For example: requiring the
gravity variables themselves to move through, and gravitate via, the same structure (a self-consistency or
"bootstrap" requirement), or a consistency condition on a conserved source. Each option is to be stated and tested,
not assumed.

## Compute

Knot optimisations about 12 CPU-min, sources 21 s: about 0.25 core-h. Linux cloud, Python 3.13.18, numpy 2.5.3,
scipy 1.18.1. Rest states are in `local_raw/EX-GC34B/` (not committed), with hashes in `states_h*/static_B*.json`.
