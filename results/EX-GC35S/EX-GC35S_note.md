# EX-GC35S — exclusion as an emergent clock

2026-11-24 | EXPLORATORY | ChatGPT | [Amended scope](../../spec/EX-GC35S_plan.md). PI closure pending.

**The ring rule gives a reproducible density–rate relation, but the density step does not support a lapse interpretation under the two proposed checks.** Exclusion is an interaction through blocking, even though legal (10→01) matches do not share a consumed site. This run amends the earlier Step 2 stop to test that rule.

The experiments use periodic rings with 32 and 64 sites, balanced and shuffled positions for each uniform density, and eight randomized schedules per fixture. Each event belongs to both sites it rewrites. The reference is a common causal cut, retaining events through the same forward-depth limit. All schedules reached identical dependency records and states on all 22 main fixtures. For the half-filled fixed-count diagnostic, eight schedules produced eight different serial prefixes at 97 events. Equal event counts do not guarantee the same cut. The GC35L negative control gives event depths (1,1) versus (1,2) on the two disputed histories.

| Uniform occupancy | .125 | .25 | .50 | .75 | .875 |
|---|---:|---:|---:|---:|---:|
| Local ticks per forward or reverse layer | .25 | .50 | 1.00 | .50 | .25 |

On the sampled densities, the rate equals 2 min(ρ,1−ρ) under this layer normalization. Every uniform site rate agrees between the two depth directions on the selected windows. This is a graph-depth rate. The result does not recover the proposed ρ(1−ρ) functional form, though it has the predicted nonmonotonic sign. A completely empty ring has no events for normalization.

The periodic step starts with .25 and .75 plateaus, with no reservoir or absorber. The measured mean occupancy moves to about .304/.696 on the 32 ring and .308/.692 on the 64 ring in the measurement windows. Interpolating the sampled uniform rates at **observed** occupancy predicts .6071 and .6154. Forward rates are .5982 and .6120, both within 1.5%. Reverse-layer rates in those same regions are .7232 and .7383, about 21% higher than forward. The full-ring forward/reverse differences were .1004 and .1107, with a half-window sensitivity of .1063 and .0915 respectively. The longer ring's difference exceeds that finite-window diagnostic, so foliation independence is not established at a 5% scale. The halves mix over time and these are finite transient windows.

The probe moved four sites while spending exactly four of its own residence-site ticks in every qualifying crossing. This is a degeneracy of the readout: the occupied site's next event is its departure. Converting crossing depths using the independent uniform calibration gave 2.00 versus 1.16 sites per calibrated tick in the two 64-ring regions (30 and 15 overlapping segments). The smaller ring did not yield a qualifying high-density crossing. Those calibrated speeds differ by 42%, above 5%. The equal residence-tick speeds alone therefore do not distinguish exclusion crowding from a lapse. Absorber-dependent clock signs were not measured.

[Code](../../src/exgc35s.py) · [Raw cuts, plateau rates and crossings](results.json)
