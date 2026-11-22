# EX-GC28H and EX-GC29A — Claude review (EXPLORATORY)

Reviewer: Claude, 2026-11-22. Reviewed ChatGPT's EX-GC28H (ea316c6–8009d35) and EX-GC29A (e170902–ff2bb3f). These
are Skyrme and scalar-loop reference models, not Model Gc. Closure is the PI's.

## Verdict

Both are correct, reproduce exactly, and are reported fairly.

- **28H** confirms the outcomes written into its plan beforehand. The even (dilaton-like) coupling gives same-sign
  sources, but they differ about 3× between species. It needs mass, tadpole and cubic tuning. Screening suppresses the
  fifth force but does not remove it.
- **29A's** split verdict is right. It is positive for an elastic roll-off at the inverse particle size. It is negative
  for a natural frequency cutoff: the response to on-shell waves does not saturate. I think that second result is
  expected for a simple kinematic reason, not a numerical one (point 2 below). That turns it from "not demonstrated"
  into "not available from particle size".

## What I checked

1. **Reruns.** I ran `src/exgc28h.py`, `src/exgc29a.py` and `src/exgc29a_dispersion.py` (numpy 2.5.3, scipy 1.18.1).
   Every number in `results/EX-GC28H/records.json` and the three 29A records matches bit-for-bit (about 1 CPU-minute).
2. **28H loop formulas.** J/E = a·B_a(0)/I_a reproduces 28C's 0.0425 for a = 1 and gives 0.1257 for a = 4 at Λ = 10.
   The cubic, a³∫d⁴q/(2π)⁴(q²+a)⁻³ = a²/(32π²) before the regulator, has the right sign and size. The fifth-force
   ratio 3π(Q/M)²/(ZG) and the bound √(εZG/3π) match my own derivation from the previous round.
3. **29A elastic form factors.** k50·RMS is 1.47 for massive P2 and 1.47 for baryon density, against 1.44
   (= √(3 ln 2)) for a Gaussian of the same RMS. So the roll-off is a generic finite-size effect, as the note says, not
   a topological one.
4. **Protocol.** For both studies the plan and LOG row came before the results. Each has its own go-ahead from you. Use
   was about 2 CPU-s (28H) and 52 CPU-s (29A), and both are labelled EXPLORATORY. The 29A dispersion diagnostic was
   registered after the first results and before it ran, and the JSON-serialisation rerun is disclosed. Both are fine.

## Comments

**1. 28H: screening behaves like a conducting sphere, so suppression is only about 1/α (adds to the note).** From the
committed screening rows (q = .25, μ = 0, R24, h = .05), the local slope d ln Q/d ln α between successive α values
settles near −0.85 (massless core) and −0.89 (massive core). It is not −1 or steeper. The reason is simple. Where αρ ≫ 1,
the field locks to φ ≈ q/α. Outside the radius r* where αρ ≈ 1, it is a free 1/r tail. So the charge seen outside is
about (q/α)·r*.

- **Massless core:** f ~ r⁻² gives ρ ~ r⁻⁸ and r* ~ α^{1/8}, so Q ~ α^{−7/8}. The slopes trend towards −0.875.
- **Massive core:** r* grows like ln α, so Q ~ (ln α)/α, approaching −1 slowly.

Two consequences:

- Each factor of 10 in force suppression costs about a factor of 3 in α, and real-G-level suppression needs absurd
  polarisabilities.
- The residual charge depends on each particle's outer profile through r*. So even at huge α, Q/M differs between
  species. Screening cannot make the force universal, only smaller. This is consistent with the note's conclusion and
  explains it.

**2. 29A: the on-shell growth is forward scattering (main point).** An incident wave of frequency ω exciting
continuum modes at the same ω is elastic scattering. That always includes zero momentum transfer, where the form
factor is 1. The radial eigenmodes and j₀(ωr) share the same asymptotic phase (up to the phase shift), so their product
has a non-oscillating part. Each mode's weight then falls only as a power of ω. The form factor never enters, because
it suppresses momentum transfer, not frequency.

This is the standard result that a finite-size object's cross-section tends to its geometric size at short
wavelength, not to zero. It also explains why the matched vacuum control grows too. The finest matched-dispersion runs
show M4 still roughly doubling from 32 to 64, after the numerical phase error was removed. So I'd read 29A as a clear
negative: particle size gives a form factor for momentum transfer, not a frequency cutoff for absorbing or scattering
vacuum modes. The note says this in words ("small coherent form factor does not bound the full excitation spectrum").
I'd state it as the expected kinematic result.

**3. 29A probe channel (minor).** The probe is an s-wave j₀, but the hedgehog breathing fluctuation lives in the
grand-spin-0 channel, whose vacuum operator has the l = 1 centrifugal term 2/r². The overlap is still well defined
because J is localised. But this is not the physical s-wave scattering of a scalar off a Skyrmion. A note sentence
would avoid it being read as one.

**4. Minor.**
- HANDOFF.md has no entries yet for 28H or 29A (or the rest of the 28 series). It should be brought up to date so a
  fresh session sees where things stand.
- 28H's "dressing vanishes at fixed q" is right (field energy ∝ qQ). Worth adding: the polarisation term αρφ²/2 is an
  explicit internal mass for φ, so screening itself breaks the shift symmetry. That is the same point as F in 28F.

## Bottom line for the PI

- 28H and 29A are correct and reproduce exactly, and both close the routes they tested.
- **28H:** neither a dilaton-like loop coupling nor internal screening gives universal sourcing. Screening only
  shrinks the extra force, roughly as 1/α, and leaves species differences.
- **29A:** a particle's size gives an elastic form factor at the inverse size, as any smooth object would. It does not
  give a frequency cutoff, because forward scattering is never suppressed.
- Both are clean candidates for closure.
- What remains is the GC30 question, where common metric coupling comes from, which neither study addresses.
