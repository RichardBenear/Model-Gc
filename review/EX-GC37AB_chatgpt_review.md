# Independent review of EX-GC37A and EX-GC37B

ChatGPT | [GC37A results](../results/EX-GC37A/EX-GC37A_note.md) · [GC37B results](../results/EX-GC37B/EX-GC37B_note.md)

**GC37A finds a positive tensor kinetic response from the scalar vacuum, with a speed that misses the GR test. GC37B gives a credible static mechanism for nearly doubling light bending relative to lapse alone.** Wave propagation under realized coincidence events and propagation coupled to the computed carrier field remain open.

## GC37A · Induced gravity

The plus-polarization vertex differs from the GC34O action specified in the plan. GC37A uses `K_x^2-K_y^2`. For a wavevector along z, GC34O's forward/backward product gives `K_x^2 cos^2(k_x/2) - K_y^2 cos^2(k_y/2)`. The static diagonalization checks the implemented action and can agree even with this difference. [GC37A plan](../spec/EX-GC37A_plan.md) · [GC37A code](../src/exgc37a.py) · [GC34O action](../src/exgc34o.py)

I independently evaluated the small-momentum bubble on a 64³ grid at m=0.1. Using the GC34O vertex and retaining the momentum-independent contact assumption moved the plus result from c_g²≈0.512 to **0.716**. The cross result stayed near **0.607**. The ordering of the two kinetic coefficients also reversed. This was a derivative check. The corrected contact term and static diagonalization were outside that check. The speed failure survives, and the reported 2.7× split pertains to the implemented vertex.

The large constant is a zero-momentum TT response. Interpreting it as a physical graviton mass requires a consistent background subtraction and checks of the full constraint and gauge structure. The continuum control supports frequency–momentum symmetry in the tested TT sector. The lapse, shift, mode count and attractive static response remain untested. [GC37A results](../results/EX-GC37A/EX-GC37A_note.md)

## GC37B · Event counting

The stochastic token result is the strongest part. The simulator samples Poisson firings and finds v∝N^1.993, with the straight-path delay following the coincidence-rate prediction. The small but statistically significant waiting-time miss matters for precision. The wave test uses a deterministic face weight N_x N_y and advances site clocks at N. The resulting effective wave operator bends the packet **2.04×** as much as lapse alone. Wave propagation under realized random coincidences remains untested. [Token and carrier code](../src/exgc37b.py) · [Wave code](../src/exgc37b_wave.py)

The absolute bend missed its registered gate. Applying the stated leading dispersion correction to the **2.873°** finite-path prediction gives about **3.05°**, with the measured **3.183°** about **4.5% higher** relative to the corrected prediction. Dispersion explains part of the gap. The static carrier fit compares two prescribed core widths at equal consumption. The computed field was not fed into the propagation run, and the test does not establish source universality across matter types. [GC37B results](../results/EX-GC37B/EX-GC37B_note.md)

## Priorities

1. Rerun GC37A's TT comparison with the GC34O vertex and a matched contact term before changing the cutoff or opening the full metric sector.
2. Feed the actual GC37B lattice carrier solution into both propagation controls and compare bend and delay without refitting GM. If that passes, test a moving source to probe the preferred-frame question.
