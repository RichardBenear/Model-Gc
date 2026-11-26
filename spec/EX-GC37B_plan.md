# EX-GC37B — gravity as an event-counting effect: clocks are one-site events, propagation needs two

2026-11-25 | EXPLORATORY | Claude | **DRAFT: awaiting PI go-ahead. Nothing has been run.** Default cap 2 core-hours; expected under 0.5. This is our construction, not LaViolette's Model G.

## The idea in plain words

Suppose the substrate is a lattice of sites that "fire" (update) at a local rate N(x). Near mass, N is lower: something matter does, such as consuming the Model Gc carrier, leaves fewer resources for vacuum updates.

There are two kinds of event:
- **Clock events** need one site. A clock at x ticks at rate N(x). That is gravitational time dilation, and 35I already showed it.
- **Propagation events** move a signal across a link, so they need **both ends awake at once**. If sites fire independently at rates N_x and N_y, the coincidence rate is proportional to N_x N_y ≈ N². So a signal's coordinate speed is c(x) ≈ N(x)².

Why that matters:
- **Weak-field GR in isotropic coordinates** has coordinate light speed (1 + 2Φ) and clock rate (1 + Φ). With N = 1 + Φ (Φ < 0 near mass), that is exactly **speed = N², clocks = N**.
- **Rulers come out right without extra assumptions.** A ruler built from local physics, for example a light round-trip timed by a local clock, has coordinate length ∝ c × period = N² × (1/N) = N. So rulers shrink by N, which is GR's spatial factor A = 1/N ≈ 1 − Φ.
- **Light bending is 4GM/b, γ = 1.** The refractive index is 1/N² ≈ 1 − 2Φ, twice the lapse-only value. 35I and 35K got only the lapse half from clocks alone, and needed the spatial metric put in by hand for the other half.
- **The equivalence principle is built in.** All matter and light move by the same handshake events, so all fall alike.

So one scalar field N, the local firing rate, would give the whole weak-field metric with γ = 1. The time part comes from one-site events and the space part from two-site events. No Einstein equation is put in. What remains is to set N from matter, for example by the Model Gc carrier deficit, which is 1/r around a consuming core (EX008).

This is a first-order (weak-field) coincidence until tested. GR's second-order terms (PPN β) and moving sources are not covered by the argument.

## Questions

- **B1 (rule):** with N(x) prescribed as 1 − GM/r, do token- and wave-level versions of the two-tier counting rule give:
  - clock rate N;
  - coordinate signal speed N²;
  - bending 4GM/b, Shapiro delay;
  - local light speed 1 measured with local clocks and light-defined rulers?
- **B2 (source):** if N is the local density of a conserved, diffusing carrier that matter consumes (Model Gc), relative to its far value, is N = 1 − GM/r with GM ∝ consumption rate? And is the same constant the same for every particle type?
- **B3 (frame), design only, not run:** a moving sink. Does the carrier field give drag or a preferred frame (EX-GC25's objection), and does a wave-like carrier fix it?

## Method

**B1a, tokens.** This is the mechanism shown directly.
- 2D square lattice. Each site fires as an independent Poisson process at rate N(x).
- A signal token sitting on x with direction d hops to x + d only when x and x + d both fire within a window τ (τ ≪ 1/N).
- A static clock at x counts its own firings.
- Fixed seeds, many tokens.
- Readouts:
  - mean coordinate speed vs N, fitting the exponent p in v ∝ N^p (expect p = 2);
  - beam deflection at impact parameters b = 10–40 cells;
  - round-trip (Shapiro) delay;
  - redshift between static clocks.
- **Control:** the one-site rule (hop when x fires) must give p = 1 and half the bending.

**B1b, waves.** Uses the 35I event framework.
- Site events at proper rate N.
- Each face's flux is applied only through a handshake factor equal to the neighbour's rate.
- Choose the face coefficient so that the coordinate wave operator is the one for speed N_x N_y. The exact discrete form is written in the code before any run.
- Readout: the 35K group-angle measurement of bending vs the lapse-only control and vs 4GM/b.
- The ratio full/lapse-only should be 2, as in 35K, but here **without** a separate spatial metric.

**B2, source.**
- 3D lattice, conserved diffusing carrier, steady state around a sink, as in the Gc1 radial solver.
- Set N = c(x)/c∞.
- Fit N vs r to 1 − GM/r and test universality across two sink types with different core sizes and the same consumption.
- Feed that N into B1a for one bending run.

## Internal checks

- **G1:** token-rule speeds at uniform N must match the Poisson coincidence formula (rate ∝ N_x N_y τ) to statistical error.
- **G2:** B1b must reproduce 35I's lapse-only bending when the handshake factor is switched off.
- **G3, separate path:** an eikonal ray tracer with the predicted speed law N² gives the expected deflection. Tokens and waves must agree with it.

## Predictions (scored)

- **P1:** token speed exponent p = 2.00 ± 0.05; control p = 1.00 ± 0.05.
- **P2:** token and wave bending = 2 × lapse-only within 3%, and within 5% of 4GM/b at the largest b.
- **P3:** the local light speed, measured with local clocks and light-round-trip rulers, is 1 within statistical error at all r.

  This is close to built in, so it is reported as a consistency check rather than evidence.
- **P4 (B2):** carrier N fits 1 − GM/r for r ≥ 3 core radii within 2%. GM is the same for both sink types within 2% at equal consumption.

## Known risks, stated before running

- **Second order:** N² for speed matches GR at first order only. PPN β (perihelion) would need N(Φ) to second order; not tested here.
- **Rest frame:** a diffusing carrier has a rest frame. A moving mass would drag and the field would not Lorentz-contract (EX-GC25). This is the deepest risk, and B3 only designs a check.
- **Lorentz invariance of the counting rule:** Poisson firing on a fixed lattice is not boost invariant. Like 37A, this tests the static weak-field mechanism, not relativity.

## What each outcome means

- **B1 and B2 pass:** a concrete, intuitive mechanism for the full weak-field metric.
  - Matter starves the vacuum's update rate.
  - Clocks feel that once (one-site events) and signals twice (two-site events).
  - That is why light bends twice the Newtonian amount.
  - It ties Model Gc's consumed carrier to Y1's lapse.
- **B1 passes, B2 fails universality:** the counting mechanism holds, but the carrier picture needs work.
- **B1 fails:** the handshake rate is not N² in a way that survives a real rule. The idea is ruled out cheaply.

## Outputs

`src/exgc37b.py`, `src/exgc37b_eikonal.py` (G3), `results/EX-GC37B/`.

Before running: an `explore/LOG.md` entry and PI go-ahead.
