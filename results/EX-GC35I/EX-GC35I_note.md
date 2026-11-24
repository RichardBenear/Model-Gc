# EX-GC35I — local event rates set by a prescribed lapse

2026-11-24 | EXPLORATORY | Design, implementation and runs ChatGPT. [Plan](../../spec/EX-GC35I_plan.md) · [Claude review](../../review/EX-GC35EF_claude_review.md). PI closure pending.

**Schedule, finite-field, frequency and angle checks passed.** For the stepped lapse one side generates 1.5× as many events over the fixed coordinate-time interval. The two static observers measure the expected redshift, and the transmitted mode turns by the expected angle. The lapse and the permitted versions are prescribed inputs.

## Event rule

Use a static time metric ds²=−N(x)²dt²+dx²+dy² and one transverse Fourier component of a massless scalar packet. The spatial discretization uses harmonic face lapse for the x flux. Writing the update in a common proper step Δτ gives a centered second-difference event for the local field. Its spatial operator is (1/N)∂ₓ(N∂ₓφ)−k̂ᵧ²φ, with the original centered transverse symbol k̂ᵧ. The update coefficients within either constant-lapse region are identical. The local versions occur every three coordinate quanta on the N=1 side and every two on the N=1.5 side. Across the interface the face flux uses both lapse values.

At each event, read the latest neighbor version at or before the current local time. If those times differ, extrapolate linearly from that version and the one immediately before it. Both inputs precede the new output time. The rule is acyclic, uses one writer per version, and the ready scheduler needs no timestep barrier. In one fine run a site had completed tick 12 while another interior site was still at tick zero. The coordinate-time reference simply orders the same events by their output timestamps.

## Measured comparison

For the static comparison the coordinate carrier frequency is conserved and a static clock measures ω/N. With conserved transverse wavenumber the continuum wave-normal angle obeys sin θ₂=(N₂/N₁)sin θ₁. Time-phase fits use local observer windows around passage of the packet. Angle estimates use the spatial phase of the incident and transmitted packets in the fixed transverse sector.

| Lapse | h | Events per run | Six randomized schedules match reference | Site event counts left/right | Measured proper-frequency ratio | Continuum ratio | Outgoing wave-normal angle | Continuum angle |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 to 1.5 | .2 | 300,900 | 6/6 | 600/900 | 0.667310 | 0.666667 | 32.911° | 32.785° |
| 1 to 1.5 | .1 | 1,201,800 | 6/6 | 1,200/1,800 | 0.667113 | 0.666667 | 32.851° | 32.785° |
| Uniform 1 | .1 | 961,200 | 6/6 | 1,200/1,200 | 1.008587 | 1.000000 | 21.159° | 21.161° |

The incident angle on the fine stepped grid was 21.162°. The measured coordinate frequencies on either side were 1.66089 and 1.66200. The measured outgoing angle differed from the continuum prediction by 0.066°. All histories were finite and all 18 randomized schedules matched their reference bit for bit at every version. The uniform control showed essentially no angle change and a 0.86% frequency-fit offset. This offset illustrates finite-packet and sampling sensitivity. The experiment estimates a wave-normal direction in one transverse sector and does not track the center of a localized beam in two dimensions.

The result connects a prescribed lapse to event counts and propagation within this fixed-version rule. Dynamic lapse, reciprocal coupling, and competing rewrite matches remain open. The run took 7.5 CPU seconds. [Source](../../src/exgc35i.py) · [Raw results](results.json)
