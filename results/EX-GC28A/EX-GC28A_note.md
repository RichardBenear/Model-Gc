# EX-GC28A: quantum vacuum fluctuation screen

**Exploratory, completed; independent review done (review/EX-GC27AB_28A_claude_review.md); PI CLOSED 2026-11-22 (Richard: "go ahead and close 28 and 29A and we will defer 29B").** Scope: one 1+1D kink and the spherically symmetric fluctuation channel of B=1 Skyrme, not a full 3D quantum soliton. Production CPU 7.98 s (<0.003 core-h), including canonical diagnostic. Source `src/exgc28a.py`; run normally and with `canonical` argument. Plans in `spec/EX-GC28A_plan.md` and `spec/EX-GC28A_regulator_check.md`.

## Findings

The zero-point difference Z_raw = 1/2 sum(omega_defect - omega_vacuum), hbar=1, is regulator dependent. Subtract first-order vacuum perturbation Z_Born = Tr[D_v^(-1/2)(D-D_v)]/4. The remainder Z is a no-tadpole energy in the kink benchmark; in Skyrme it is only a single-channel diagnostic, not a physical fully renormalized mass. Collective translation mode in kink is identified by overlap >0.99 and assigned zero frequency. Hessian finite difference error 5.27e-11; dense/generalized eigenvalue check 4.26e-14.

| h | phi4 kink (R20) Z | Skyrme midpoint/nodal C0 Z (R24) | canonical C0 Z (R24) | canonical C=.178396 Z (R24) |
|---:|---:|---:|---:|---:|
| .2 | -.673439 | -.193113 | -.175984 | -.185057 |
| .1 | -.668918 | -.254072 | -.169261 | -.177393 |
| .05 | -.667131 | -.357393 | -.166827 | -.174593 |
| .025 | -.666526 | -.554010 | -.165998 | -.173631 |

Kink exact no-tadpole result = mu*(1/(2 sqrt(3))-3/pi) = -.666255 at mu=1. Finest grid error .041%; raw sum changes from -3.5091 to -5.5122. The original discrete Skyrme midpoint Hessian with nodal kinetic metric has a residual power-like divergence; do not interpret it as physics. Canonical diagnostic solves continuum static BVP with f(0)=pi, f(R)=0, then u=sqrt(K) delta f, K=A r²+2B sin²f. Operator -u''+[(sqrt K)''/sqrt K+(U_ff-K_ff f'^2-2K_f f'')/(2K)]u, vacuum -u''+(2/r²+C/(2A))u. Identical centered second derivative in both channels. Raw C0 result changes -1.779 to -3.288 while Z stabilizes; R12 vs R24 fine-grid difference 8.55e-5. Massive case difference 3.26e-6. This supports regulator control in one radial channel only.

## Scale response

Eta perturbation ±.02; A=exp(2 eta), B=1; consistent C=Cbase exp(4 eta), independent C=Cbase. Centered derivatives dZ/deta at fixed h,R24:

| Cbase | h | consistent | independent | comoving consistent |
|---:|---:|---:|---:|---:|
| 0 | .1 | -.174865 | -.174865 | -.169273 |
| 0 | .05 | -.168769 | -.168769 | -.166839 |
| .178396 | .1 | -.183871 | -.154304 | -.177405 |
| .178396 | .05 | -.176890 | -.148570 | -.174605 |

The comoving response equals Z times sinh(.02)/.02 as dictated by the imposed scale transformation; this is a consistency check, not evidence of gravitational universality. Independent potential scale produces a modest but resolved change in the derivative (~16% vs consistent at h=.05), and the baseline radial energy shifts ~4.6% between cases. No comparison to classical energy in physical units: quantization normalization and renormalized counterterms must be fixed first. No B=2, full angular channels, translational/rotational collective quantization, stress tensor, active sourcing, or gravitational potential measured. No evidence yet that ZPF is the source of gravity.

## Recommendation

EX-GC28B: test the *connected* vacuum energy of two separated localized scattering defects with matched regulator, subtracting both isolated-defect terms and adding vacuum back. Begin with a controlled 1D Gaussian potential pair in a linear scalar field; compare sign, distance dependence, massless vs massive mediator, and strength dependence, then proceed to 3D compact defects if justified. This probes whether a vacuum-mediated interaction emerges without imposing attraction or a coupling to rest energy. Casimir-like attraction is not itself gravity: test for long range, additivity, and source/inertial-energy proportionality rather than assuming them. Do not tune a regulator-dependent self-energy to match gravity.
