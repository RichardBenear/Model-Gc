# EX-GC61 — the shift as a travelling field: local events, arrival, settling and the detector without refitting

2026-11-28 | EXPLORATORY | Claude | Plan only. PI asked for it to be written ("Yes, do 1 and write the 2 plan and also 3 addon", 13:52 UTC). PI go-ahead 2026-11-28 13:58 UTC ("go ahead with 61"). Requested cap: 2 core-h (expected about 1). Our construction; not LaViolette's Model G and not Chain Theory.

## Why

In 57–60 the shift β was always solved statically (a Poisson solve in a zero-boundary box) and then handed to the probe. A static solve is instantaneous action at a distance, so it says nothing about how a spinning source's influence gets to the probe. The Codex 60 review's main recommendation was:
- switch on a localised source current;
- let its shift propagate by a local law;
- measure when it arrives at several distances;
- compare the settled field and the axis and diagonal detector rates with the static results, without refitting;
- vary the box to separate propagation from reflection.

For the readout, total J is the main benchmark and raw Σ is secondary. The review also asked for three detector checks that 60 left open (61D).

What this does and does not test:
- **It tests** whether the shift sector can be carried entirely by local, order-free events at the light speed, and still give 59's and 60's static results where they should hold, and free-space GR where the static box could not.
- **It does not** derive λ = 4 or the source. Those stay as put in (53, 58). The wave law below is 53's `evolve_b` law, which 53 wrote by hand as a PDE. Here it is recast as events, and its consequences are scored for the first time.

## The rule (worked before any run)

Each link l = (x, a) carries an offset β_l and a rate π_l. Weights are N̄_l = ½(N_x + N_{x+a}), as in 53, and the coincidence weight is N̄_l N̄_l′, as in 49's N rule. Each tick of length dt has three classes of event:

1. **Drift (one-link events).** β_l += (dt/2) N̄_l π_l.
2. **Coincidence (two-link events).** For each pair of parallel neighbouring links l, l′ = l + e_b (any b), with Δ = β_l′ − β_l:
   - π_l += dt N̄_l² N̄_l′ Δ;
   - π_l′ −= dt N̄_l N̄_l′² Δ.

   **Source (one-link events).** π_l −= dt N̄_l λκ j_l(t).
3. **Drift again**, as in class 1.

Within each class every event reads only variables that the class does not change, so the events commute: any order within a class gives the same result exactly. The tick is a Strang (leapfrog) step. Its continuum limit is 53's law:
- ∂_t β = N̄π;
- ∂_t π = N̄(Σ N̄N̄′ (β′ − β) − λκ j).

In the static limit this is 59's ∇²β = λκ j, and the signal speed is N̄² = light's (37B, 49).

Faces: β = π = 0 on the box faces, the same Dirichlet condition as 59's static solve.

**Source.** 59's polarised packets (n = 97, m = 0.6, σ = 3, spin ±z), with j = P + ¼Cᵀs (rule PS), are held fixed. Their current is ramped as j(t) = j·r(t), with r = sin²(πt/2τ) for t < τ and 1 afterwards, and τ = 8. N is 59's static N for the source. The probe is 59/60's (σ = 3, m = 0.6, spin along x, the PS coupling), and the ±S odd part removes everything that is even in the source spin.

**Expected behaviour in a 3D continuum.** Behind the front (t > r + τ) the field of a ramped static source is exactly static. Before any reflection gets back, it is the **free-space** field, not the zero-boundary box field that 59 used. Reflections from the nearest face reach the axis site d = 16 (source at the centre, n = 97) by the image path, a distance of 80, at about t ≈ 80 + τ/2.

So the dynamic run makes a sharp prediction: the 2.3% zero-box deficit in 59's on-axis rate is absent before t ≈ 80 and appears only afterwards.

## Parts

**61A. Event rule against the PDE, and order freedom** (n = 65, small source).
- Compare the event leapfrog with 53's RK4 `evolve_b` (same N̄, λ, κ, j, initial zero field) at t = 20, at dt = 0.1 and 0.05.
- Run the same ticks with events shuffled at random within each class, with two seeds.

**61B. Arrival and settling** (n = 97 and n = 129, field only; dt = 0.1 to t = 100).
- Record β (the x and y components) at sites on the axis (0, 0, d) and in the plane (d, 0, 0), for d = 8, 12, 16, 24.
- Arrival time t½(d): the first time the odd part of β at the site reaches half its value at t = 50.
- Settled field: compare the dynamic β at t = 50 (n = 97) with the static solve in an n = 193 zero box (near free space; the image deficit there is about 0.3%), over shells 8 ≤ |r| ≤ 20.
- Boxes: compare n = 97 and n = 129 at the axis site d = 16.

**61C. The detector in the travelling field, with nothing refitted** (n = 97 and n = 129; probe window t ∈ [40, 50]).
- Probe runs at the axis site (c, c, c + 16) and the diagonal site (c + 12, c, c + 12), with source up and down. β is taken from the evolving field, interpolated linearly between ticks for the RK4 stages.
- Readouts: the J angle (primary) and the Σ angle (secondary) about z, for the ±S odd rate.
- Predictions use only frozen numbers: 60A's σ = 3 coefficients (Σ: k_rot −0.99758, k_sh +0.09274; J: −0.99754, +0.00198) with the probe-density-weighted lattice gradients of the dynamic β at t = 45, and free-space GR from 59's `omega_gr`.
- A late window, t ∈ [90, 100] at n = 97, on the axis only, after the reflection has arrived.

**61D. Detector add-on** (60's flat 48³ box and probe; dt = 0.1, T = 10).
- **D1. Matched physical readout.**
  - The resync-description detector is Σ′ = U Σ U†, read as ⟨U†ψ|Σ|U†ψ⟩, with U = e^{−iG}, G = ½{χ, H}, χ = −s xy (60B).
  - It is applied to four runs at σ = 3: the H_g run, the P-coupled xy shear, the sum field (rotation + shear), and the H_g and P shear at 2s.
  - Every sample is saved: each Σ component and the angle at all 101 readout times, under both the raw and the transformed detector.
  - The same matched readout is run for H_g and P shear at σ = 4 and 5.
- **D2. Larger box.** σ = 5 shear, P and PS, at NB = 64.
- **D3. Tilted spin, general gradient, off-symmetry site.**
  - (i) Flat box. The probe spin is along n̂ = (1, 1, 1)/√3. The field is a predeclared uniform gradient β_a = s G_ab x_b with all nine entries non-zero. G is fixed now as
    ```
    G = [[ 0.3, -0.7,  0.4],
         [ 0.5,  0.2, -0.6],
         [-0.1,  0.8, -0.5]]
    ```
    Run it with PS and with a β = 0 reference, and read the Σ and J vectors.
  - (ii) 59's static source field (n = 97, zero boundaries) at the off-symmetry site (c + 9, c + 5, c + 13), with probe spin along n̂ and ±S. The 3 × 3 density-weighted gradient matrix of β at the site is computed before the run.
  - Vector prediction from 60's two coefficients, assuming isotropy:
    - dS/dt = k_rot Ω̃ × S + k_sh [Sym·S − (Ŝ·Sym·Ŝ) S], where Ω̃ = −½∇×β is the rotation that 60 normalised to −1 per unit s, and Sym = ½(∂_b β_a + ∂_a β_b);
    - this is the only form linear in Sym and S whose perpendicular part is covariant;
    - the same with the J coefficients.

    A cubic lattice allows the diagonal strain (Sym_xx-type) to have its own coefficient. If D3 fails, that coefficient is the first suspect, and the failure measures it.

## Predictions (scored)

- **P1 (event rule = the law).**
  - (a) The event leapfrog matches 53's RK4 PDE solution at t = 20 to 10⁻³ of max|β|, and the difference falls by 4 ± 1 at dt/2.
  - (b) Shuffled orders give identical fields, to 10⁻¹² relative.
- **P2 (arrival at the light speed).**
  - t½(d) − d is the same at d = 8, 12, 16 and 24 to within ±0.75, on the axis and in the plane.
  - A straight-line fit of t½ against d has slope 1.00 ± 0.03.
- **P3 (settles to the free-space field).** At t = 50 (n = 97), the rms relative difference between the dynamic β and the n = 193 static field over 8 ≤ |r| ≤ 20 is below 1%.
- **P4 (boxes: propagation vs reflection).**
  - At the axis site d = 16, the n = 97 and n = 129 fields agree to 10⁻⁴ relative for t ≤ 70.
  - Their difference first exceeds 10⁻³ relative at a time between 76 and 92.
- **P5 (on-axis rate without refitting).**
  - In the t ∈ [40, 50] window (n = 97), the J and Σ odd rates at the axis site are 1.00 ± 0.015 of free-space Ω_GR.
  - The same rates are 1.010 to 1.035 times 59's static zero-box rate (59 main axis_16: 0.977 of GR).
  - n = 129 gives the same window rate to 0.2%.
- **P6 (diagonal rate without refitting).**
  - In the same window, the J and Σ rates at the diagonal site equal the frozen two-coefficient prediction from the dynamic field's gradients at t = 45, within 1%.
  - J also equals the rotation-only prediction within 1%.
- **P7 (after reflection).** In the late window t ∈ [90, 100] (n = 97), the axis J rate differs from the [40, 50] rate by more than 0.5%. This shows the detector sees the reflected field, i.e. that P5's agreement is not insensitivity.
- **P8 (matched readout, D1).**
  - (a) In the H_g run, the transformed detector's Σ angle stays at the free (β = 0) history: the maximum over all 101 samples is below 10⁻³ of the raw shear signal at t = T.
  - (b) That maximum residual grows by 4 ± 1 when s doubles (second order).
  - (c) Under the transformed detector, the P-shear rate is k_P − k_Hg (from the raw rates of the same runs) to within 5% of that difference.
  - (d) Under the transformed detector, the sum-field rate is k_rot + (k_P − k_Hg) to 10⁻³ of k_rot. Sensitivity to rotation is kept.
  - (e) The residual fraction (k_P − k_Hg)/k_P falls with width as σ⁻²: ratios at σ = 4 and 5 to σ = 3 of 0.56 and 0.36, within 25% each. This is the first independent prediction of the 7% residual, from sin k − ½ sin 2k ≈ k³/2.
- **P9 (larger box, D2).**
  - At NB = 64, |k_sh(P) − k_sh(PS)|/k_sh < 10⁻⁷ (it was 7.3 × 10⁻⁶ at NB = 48).
  - k_sh(σ = 5) at NB = 64 equals NB = 48 within 0.1%.
- **P10 (vector prediction, D3).** For (i) and (ii), the measured perpendicular rate vector of J matches the predicted one to within 3% of the predicted norm. The same holds for Σ.

## Checks

- **C1.** With j = 0, the field stays at zero to 10⁻¹⁴.
- **C2.** The field energy ½Σ(π² + Σ_pairs N̄N̄′Δ²) is conserved after the ramp (t > τ), with the source term counted, to 10⁻⁴ relative over t ∈ [10, 100].
- **C3.** With the dynamic field replaced by 59's static zero-box β, the 61C code reproduces 59's axis d = 16 Σ rate to 10⁻¹⁰.
- **C4.** The ratio ∇·j / |j| of the source current, reported. The wave law does not need ∇·j = 0, but the static Poisson comparison assumes the same j.

## What the outcomes mean

- **P1–P7 pass.** The shift needs no instantaneous solve.
  - A local, order-free event rule carries it at the light speed, and it settles to the free-space field.
  - Before any reflection, a probe reads the free-space GR precession without refitting. 59's 2.3% box deficit is shown to be a static-solve artefact, and it returns once the reflection arrives.
  - The whole gravitomagnetic chain (source → events → probe) is then causal and local, with λ and the source law still put in by hand.
- **P2 fails.** If the slope is below 1, lattice dispersion is slowing the front at the ramp's bandwidth. If the offset is not constant, the source's size is setting the arrival.
- **P3, P5 fail together.** The dynamic field has not settled by t = 50: a lattice transient or a longitudinal (∇·j) part. C4 and the shells show which.
- **P8 passes.** The transformed detector removes the pure-coordinate drift and keeps rotation. The 60 claim then becomes a demonstrated statement for this detector. If (e) also passes, the residual is the lattice current difference.
- **P8 (a) or (d) fails.** The "readout convention" reading is wrong in part, and Σ′ is not the right gyroscope.
- **P10 fails.** The detector needs a third, cubic, coefficient, and D3 measures it.

## Cost

Measured from 60: a 48³ probe run takes about 60 s, and a matched readout about twice that if read at all 101 samples.
- 61A and 61B: field only, about 5 min.
- 61C: 10 probe runs, about 12 min.
- 61D: about 25 min.
- Total: about 0.8–1 core-h. Requested cap: 2 core-h. If it would go past the cap, I stop and ask.

## Outputs

`src/exgc61.py`, `src/exgc61_summary.py`, `src/exgc61_recompute.py` (P5 and P6 from an independent spin readout and a polyfit; arrival from a cubic-spline crossing), `results/EX-GC61/`, note, LOG row, HANDOFF section.
