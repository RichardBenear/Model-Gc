# Model Gc definition (version 1, 2026-11-18)

Status: reference definition, written by Claude at the PI's request (Richard Benear, 2026-11-18: "set up what needs to
be moved and added"). The name was proposed by Claude, and the repository name adopts it. Model Gc is our own
construction, not LaViolette's.

## Design principle
A field with 1/r range in 3D needs three things:
1. The carrier (G) is conserved in the bulk: there is no source or loss in the ether, and the supply comes from far
   away (G → G∞).
2. Net consumption happens only in the particle core, through a species that is zero in the ether (W).
3. The particle drifts toward lower G, so the far-field well attracts.

The steady G equation outside the core is then ∇²G = 0. The only decaying solution is G∞ − M/(4πr), where
M = ∫k·G·W dV is the consumption rate. A pair approaches at v(d) = μM/(4πd²) per particle, where μ (< 0) is the l = 1
mobility to a G gradient.

## Gc1 (gravity only)
```
X_t = D_X ∇²X + A − (1+B)X + X²Y + c_x·k·G·W
Y_t = D_Y ∇²Y + B·X − X²Y + c_y·k·G·W
W_t = D_W ∇²W − W + c·Y·W² − W³
G_t = D_G ∇²G − k·G·W
```
| Symbol | Value used so far | Meaning |
|---|---|---|
| A | 140/9 | X supply (Brusselator); ether X₀ = A |
| B | 29 | X → Y conversion; ether Y₀ = B/A |
| D_X, D_Y, D_W, D_G | 1, 12, 0.1, 1 | diffusion constants |
| c | 1.0 | W autocatalysis on Y. W = 0 is stable in the ether because c·Y₀ = 1.86 < 2 |
| k | 0.05–0.5 (≥ 0.001 tested) | G consumption rate in the core |
| (c_x, c_y) | (1, 0) or (0, 1) | consumed G becomes X ("G feeds X") or Y ("G feeds Y") |
| G∞ | 1 | far background |

Boundary conditions: zero flux for X, Y and W; G = G∞ at the outer boundary.

Reaction reading: A → X, X → Ω, B + X → Y + Z, 2X + Y → 3X (as in Model G, without the G step). W is sustained by Y
(gain Y·W²), with losses W and W³. G + W → W + X (or Y) at rate k·G·W.

## Gc2 (gravity plus charge, EX010)
Adds a conserved signed carrier C, fixed at C∞ far away, and two mutually exclusive core species W₊ and W₋.
- "+" cores make C: Y + W₊ ⇌ C + W₊, with tuned kf = 0.053588 and ε = 1.37951.
- "−" cores consume C: C + X + W₋ → W₋.
- G is consumed by both types, as in Gc1.

At the tuned point Q₊ = −Q₋ and μ₊ = −μ₋. Reciprocity is tuned at one operating point, not built in.
Code: `src/ex010_network.py`.

## Known properties (with sources)
- A stable 3D particle exists at the values above. Leading eigenvalues are negative apart from the translation mode
  (EX008; MG195 radial checks).
- There is no 1D localized particle for this chemistry at B = 29–30.1 (EX008).
- The far field is G∞ − M/(4πr) to 1e-6, independent of domain size (EX008).
- The pair law holds to 1–2% for d = 22–40. Below d ≈ 12 the short-range Turing force dominates, and it is still
  visible at d = 18–22 at weak k (MG195).
- The force scales as k² at fixed 1/r range (EX009).
- Screening control: adding λ(G∞ − G) gives the Yukawa form, matched to 1% (MG195 C2).
