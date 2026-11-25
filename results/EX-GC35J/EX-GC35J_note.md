# EX-GC35J — soft-constraint Skyrme compatibility gate

2026-11-24 | EXPLORATORY | Design and runs ChatGPT | [Plan](../../spec/EX-GC35J_plan.md) · [Claude review](../../review/EX-GC35GHI_claude_review.md). PI closure pending.

**The direct soft-constraint momentum identity failed.** Relaxation and moving knots did not proceed. The potential can be differentiated locally in time, but the spatial Ward identity also requires an exact discrete chain rule across each link.

For a link x→y, the 34O endpoint average in the momentum-density time derivative supplies (∇Vₓ+∇Vᵧ)·(uᵧ−uₓ)/(2h). The volume-factor stress contributes −(Vᵧ−Vₓ)/h. These cancel for quadratic V. For V=μ(|u|²−1)²/4, take uₓ=e₁ and uᵧ=2e₁. The first term is 3μ/h and the second is −2.25μ/h. The remainder is 0.75μ/h even though both endpoints are nonzero. An exact discrete gradient could reproduce the potential difference but would change the construction and would still need a compatible action and momentum law.

Across random four-component fields and two spacings, the quadratic control had relative momentum defects at most 5.51e−15. The soft potential's defect was 0.713 relative to the larger norm of the two potential-dependent terms, at each of μ=1,4,16. The terms were tested against the independently computed local-acceleration contribution. The temporal chain rule agreed to 1.11e−10 against centered finite differences. The norm here is a rough-field identity screen and is not a prediction for smooth relaxed knots. In the 35H collinear sign-flip control both the soft force and potential differences vanish. The unrestricted sigma evolution then has no hard-constraint reaction term, though that profile is no longer a static solution.

The 35H quartic candidate also fails a separate static variation test against the 34O link generator on random fields. The normalized divergence mismatch was 1.003 and the rank-one null control was exact. This test uses the flat quartic Euler variation, rather than assuming a dynamical quartic momentum density. The potential failure already blocks the proposed direct coupling before that additional derivation is needed.

[Source](../../src/exgc35j.py) · [Raw results](results.json)
