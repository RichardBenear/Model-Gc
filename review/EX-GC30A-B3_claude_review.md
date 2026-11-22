# EX-GC30A, 30B, 30B2, 30B3 — Claude review (EXPLORATORY)

Reviewer: Claude, 2026-11-22. Reviewed ChatGPT's EX-GC30A (f9aedde–5f58068), 30B (e13555f–fe4f42c), 30B2
(56da099–c470084) and 30B3 (de961dc–c1846a3). These are reference analogue-medium, lattice and field-theory models,
not Model Gc. Closure is the PI's.

## Verdict

All four are correct, reproducible and fairly reported. I found nothing that should change the notes before closure.
The series is a clean sequence of negative results, each with an explicit mechanism:

- **30B (static tuning):** a common metric in a two-component condensate needs λ = n·g12, which background density
  breaks.
- **30B2 (topology):** Weyl windings survive, but the metric tensors split.
- **30B3 (RG locking):** speeds converge only logarithmically in a weakly coupled theory, and not at all for a
  decoupled sector.

30A's framing is sound and its no-go discussion is accurate. As before, most results reproduce known literature
(Liberati–Visser–Weinfurtner, Volovik/Zubkov, Anber–Donoghue). The notes say so.

## What I checked

1. **Reruns.** I ran `src/exgc30b.py`, `exgc30b2.py` and `exgc30b3.py` (numpy 2.5.3, scipy 1.18.1). 30B3 Yukawa
   matches bit-for-bit. Everything else agrees to rounding: the largest difference is ~3e-12 absolute, on quantities
   that are themselves fit errors of order 1e-12.
2. **30A two-cone counterexample.** c±² = (n/m)(g ± g12), and the ratio at g12/g = 0.1 is √(1.1/0.9) = 1.1055.
   Correct. The Weinberg–Witten statement and the list of possible escapes are accurate, including that GR itself
   escapes through the absence of a gauge-invariant local gravitational stress tensor.
3. **30B formulas.** Expanding (ε+2λ)(ε+2λ+2n(g−g12)) against ε(ε+2n(g+g12)) gives k² coefficients 2λ + n(g−g12)
   and n(g+g12). They are equal iff λ = n·g12, and the mismatch is D_B = 2|n·g12 − λ|/(ng+λ). That gives 8.0% at
   n = 0.8, as tabulated.
4. **30B2 formulas.** B_z = 4t_z²a²[1−(m±a)²] at the two root families, with difference 16t_z²a³|m|. At m = .025
   that is 0.05 on a tensor norm of ≈1.60, i.e. 3.1%, matching the reported 3.124%. Flavor B_z = 1−m² gives .96 and
   .64 untuned, as stated.
5. **30B3 RG.**
   - α = g²/16π²c³ with β_g = 5g³/(16π²c³) gives α′ = −10α² in the IR.
   - The relative-speed eigenvalue 7g²/(24π²c³) is 14α/3.
   - So η = η₀(1+10α₀t)^(−7/15).
   - At α₀ = .03 and 40 decades (t = 92.1) the factor is 0.209, matching 0.00209 for η₀ = .01.
   - QED is infrared-free, so its weakening coupling and the stall of the weakly charged sector are correct.
6. **Protocol.** Each study has a plan and LOG row before its results, a quoted go-ahead from you, and negligible
   compute (about 23 CPU-s for the whole series). All are labelled EXPLORATORY. Two deviations were registered before
   their reruns: the 30B integer-dtype bug, a genuine error caught by the pre-registered analytic control, which is the
   protocol working as intended, and the 30B2 import fix.

## Comments (no change required)

**1. 30B3 is missing one prior reference.** The "Lorentz invariance as an infrared fixed point" idea goes back to
Chadha & Nielsen (Nucl. Phys. B 217, 125 (1983)), who also found the logarithmically slow approach. The sharper
objection is Collins, Perez, Sudarsky, Urrutia & Vucetich (PRL 93, 191301 (2004)). With a physical cutoff, loops feed
Lorentz violation into low-energy propagation at order one rather than suppressing it. So RG locking has to beat both
slow running and radiative re-generation. That strengthens 30B3's negative conclusion and is worth citing in the next
paper-only screen.

**2. The gates were strict but fair.** "All symmetry-preserving variations within 1%" is demanding. It is the right
test, though, because the real equivalence-principle bounds are far tighter. The notes are careful to separate
conditional kinematic successes (shared shear or flow at the tuned point, symmetry-related Weyl nodes) from the failed
robustness gate.

**3. The common thread for GC31A.** In each model, a common metric holds only where a symmetry maps one excitation onto
another. Inequivalent sectors never share one robustly. That suggests GC31A's screen should ask one question first:
does the candidate pregeometric structure make all excitations propagate on the same underlying relation (one causal
order or one graph), so that a shared cone is automatic? Or does it still allow sector-specific propagation rules on
that structure? 30A already gave a counterexample to the second case: two walkers on the same graph with different
amplitudes. A causally invariant rewriting rule, the parent ModelG programme's original goal, falls in the first
class. It may be worth stating that link in the GC31A plan.

## Bottom line for the PI

- All four studies are correct and reproduce, and I see no reason to change the notes. They can be closed as you
  proposed.
- Together they show that a robust common metric for genuinely different excitations does not come from any of:
  - tuning a medium,
  - momentum-space topology,
  - weakly coupled RG flow.
- The candidates left are a strongly interacting infrared fixed point and pregeometric structures. ChatGPT's proposed
  next step (GC31A) is the latter.
