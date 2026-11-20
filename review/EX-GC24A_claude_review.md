# EX-GC24A — Claude review (EXPLORATORY)

Reviewer: Claude. Date: 2026-11-20. Reviewed commit: 07afc7b (ChatGPT, `spec/EX-GC24A_plan.md`,
`src/exgc24_impedance.py`, `results/EX-GC24/EX-GC24_note.md`, `results/EX-GC24/diagnostics.json`).
Model Gc is our construction.

## Verdict

The result holds. The wave medium does give the moving source a reactive (inertia-like) response, but it is far too
small to matter against the head-on absorption drag inherited from the recoil closure. Coasting time comes out at about
0.04 (strong source) and 0.002 (weak source), against the hundreds of time units needed. The added mass grows as S²
while the gravitational source grows as S, so it could not be a universal inertial mass even if it were large. I agree
with the conditional negative, with stopping Stage B, and with deferring the 14 core-hour h = 0.125 run.

## What I checked

1. **Reran the whole script on the cloud container** (Python 3.11, newer numpy than the PI's PC).
   - Spectral impedance rows: identical to the committed file.
   - Time-domain (RK4) fit: matches to 5e-11.
   - Linearised 1/G-weighted harmonic: matches to 3e-11.
   - Steady moving-source controls: match to 3e-13 for every v ≠ 0. A first comparison showed a relative difference
     of 1.4, but that came only from the v = 0 rows, whose values are round-off (about 1e-17). It is not a
     discrepancy.
   - The A5 centroid check was not rerun: it needs a K = 0.01 settled state that exists only on ChatGPT's side
     (`local_raw/` is not in git). I accept the committed value (ratio 0.999999) unverified.
2. **Separate path.** Before the report I had computed m and b at τJ = 20, L = 32 with my own FFT code, written
   separately from the shared script. I got m = 0.2308, b = 0.01158; the note has 0.23074 and 0.01158. The only shared
   dependency is the parameter set and the Gaussian source shape.
3. **Analytics.** The finite-damping mass formula m0 = Σ W c²/a² [1 − 1/(a τJ²)] agrees with a small-ω expansion of
   the impedance. It gives m_add → S²/(6π^{3/2} c² σ) as τJ → ∞, and b_wake ≈ m_add/τJ, as in my EX-GC24 proposal
   review.
4. **Sign at strong damping.** At τJ = 0.25 the reactive mass is negative (m/b = −5.3). This matches EX-GC17, where
   the purely diffusive token gas gave an anti-inertial response. Two separate constructions agree on this, which is a
   useful cross-check: the diffusive limit pushes the body forward when it accelerates, and only a wave-like medium
   gives positive mass.

## Points of agreement and caution

- **A4 is decided by the closure, not by the dynamics.** The head-on term b_bare = S·v is the inherited closure: a
  body moving through an isotropic token gas meets more tokens from ahead. The gate just compares two numbers that are
  both fixed analytically, so A4 confirms the arithmetic rather than testing anything new. The note does say this.
  In my view (reasoning, not tested) any absorber at rest in an isotropic medium has this Doppler-type imbalance,
  so the drag is not specific to our closure. The drag goes away only if the medium has no rest frame.
- **Retaining the 1/G weighting roughly doubles the wake damping** (b rises from 0.0116 to 0.025, m barely moves) at
  S = 5.32. It does not change the conclusion.
- **The h = 0.125 full run.** Deferring it is right: it would only refine a slip that is already small. If the PI
  wants it, it fits in the PC allowance (about 14 of the remaining ~35 core-hours).

## Protocol notes

- There is **no `explore/LOG.md` entry for EX-GC24A**. The plan (`spec/EX-GC24A_plan.md`) was committed in the same
  commit as the results, so the git history can't show that the plan came before the run. This is a protocol
  deviation, and I think it is minor because the gate (b_total against τ·m) was already written in the EX-GC24
  proposal and my review of it before the run. It should be recorded in the log as a deviation, not fixed afterwards.
- Compute used is small (well under the 2 core-hour cap). The note should state the core-hours explicitly.

## What this means for the programme (the PI's decision)

EX-GC16 to EX-GC24 now give a consistent picture. Gravity from net consumption needs a medium with a rest frame, and
that medium then gives drag. No added mass the medium produces comes close to cancelling the drag. Inertia that works
in every frame needs the opposite: no preferred frame. I see two directions for EX-GC25:

- **(a) Internal inertia.** Keep the medium as it is and give the particle its own inertia, for example from the
  knot's internal (topological) dynamics. The drag stays, so the particle would still slow down over a set time. The
  test is whether that time can be long and the same for every particle type.
- **(b) A frame-free (Lorentz-compatible) reformulation.** Make the medium relativistic, so that a moving soliton is
  just a boosted copy of one at rest. Inertia then comes for free (E/c²). The price is that gravity can no longer come
  from net absorption: it would have to couple to the soliton's field energy. This is a larger change to Model Gc.

(a) is cheaper and keeps the existing results. (b) goes more directly at the programme's goal of a link to quantum
mechanics and gravity. My recommendation is a short analytic look at (b) before investing in (a), because if gravity
can't be kept in (b), that settles the choice.
