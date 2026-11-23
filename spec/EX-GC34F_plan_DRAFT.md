# EX-GC34F (PROPOSED, NOT APPROVED) — stability of the 34E lattice geometry (EXPLORATORY)

Claude, 2026-11-23. Draft for the PI's go-ahead; nothing here has been run. Prompted by ChatGPT's review of white
paper v0.5 and the post hoc confirmation in `results/EX-GC34E/review_unstable.json`. Our construction, not
LaViolette's Model G. Proposed cap: 2 core-hours on the cloud (expected < 0.5: Fourier-space linear algebra plus one
small time-domain run).

## Question

The 34E vacuum equations have a growing mode in the re-slicing direction (lapse plus lengthwise shift), with
Im ω ≈ 1.4–2.1 k off the lattice axes at every spacing tested. Is this real in the time domain, which term admits it,
and can a stated change remove it without losing 34D's static result and 34E's gravitomagnetism?

## Parts (criteria to be fixed before any run)

**F1 — complete spectrum.** A constraint-reduced complex-frequency calculation: deflate the exact null directions
(the displacement gauge, plus the exact re-slicing direction where present) explicitly rather than through a pencil,
and keep every complex root. Directions: the axes, the face and body diagonals, and a grid over the irreducible
1/48 of directions. Spacings kh = 0.025, 0.05, 0.1, 0.2. Report each root's multiplicity and polarisation, and check
every eigenvector against all 11 equations (residual ≤ 1e-10). Registered expectation: the growing branch reproduces
`review_unstable.json` to 1e-6, and the axis directions are decided.

**F2 — time-domain check.** Integrate the linear equations for a single Fourier mode (generic direction, kh = 0.05)
in continuous time with lapse and shift solved from their own equations at each step, starting on (i) the growing
eigenvector and (ii) a gravitational-wave eigenvector. Expectation: (i) grows at Im ω within 1%; (ii) oscillates at
the F1 frequency with no growth above round-off, or picks up the growing mode at a rate set by its round-off seed.

**F3 — which term.** Turn off or vary, one at a time: the trace coefficient in the DeWitt form (K:K − λ K²), the
per-tetrahedron strain fit (averaging over the 6 tetrahedra of a cube instead of each separately), and the
lapse–Regge coupling. Report how Im ω and the re-slicing breaking (34E E2) change.

**F4 — a candidate fix (only if F3 points to one).** A stated alternative kinetic action or a constraint-preserving
treatment of the re-slicing direction. It must keep: no growing mode at any direction or spacing tested, the 34D
static kernel unchanged (≤ 1e-12), 34E E4 gravitomagnetism ≤ 1e-2 at kh = 0.05, and two wave modes within 1e-3 of the
light speed.

## Scope

Linear, flat background, continuous time. No knot, no non-linear dynamics. A negative F4 is an acceptable outcome; it
would mean the Regge-plus-DeWitt route needs a different discretisation of time slicing before any non-linear work.
