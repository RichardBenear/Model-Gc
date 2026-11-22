# EX-GC28E — self-consistent secondary-field dressing (EXPLORATORY)

**Summary:** A topological B1 radial Skyrmion dressed with its own massless scalar field has a conserved static energy including positive field energy and negative interaction energy. Local polarization changes the far monopole by up to ~76% suppression or ~28% enhancement at q=.25, with <.2% total energy changes. Translation kinetic inertia agrees with rest energy on the stationary continuum branch within finite-box errors. Feedback does not enforce gravitational source proportionality, and exact sign symmetry permits equal-energy opposite sources or zero source. This is a phenomenological classical model of vacuum response, not computed ZPF or real gravitational shielding.

Design/run ChatGPT, PI closure/review pending. Plan `spec/EX-GC28DEF_plan.md`; source `python src/exgc28def.py E`. 40 converged BVPs, 12.45 CPU-s one thread. C=0 or .178396, rho=(1-cos f)^2, Z=1 rescaled, static energy as plan. No extra GC28A quantum self-energy added. f endpoints pi,0, regular scalar derivative and Coulomb Robin outer boundary; include analytic exterior gradient energy. Boundary small-r cutoff 1e-4, R12/24. BVP tolerance 2e-7 and tighter 5e-8 subset. Residual ≤2e-7. Radial scalar source Q=integral r² rho(q-alpha phi)/(3pi), amplitude A=3pi Q. Absolute flux discrepancy ≤3.51e-10, integrated scalar energy ledger 2(Efield+Epol)+Eint within 8.97e-12. Negative alpha sufficient scalar positivity bound ≤.341<1 (fixed profile; no full coupled stability). No nonradial or dynamical stability assumed.

| C=.178396, q=.25, R24 | alpha=0 | alpha=1 | alpha=5 | alpha=-.25 |
|---|---:|---:|---:|---:|
| Etotal |1.26089|1.26159|1.26216|1.26050|
| Q |.0183838|.0108263|.00444351|.0225878|
| Qbare (profile only)|.0183838|.0182189|.0181194|.0184976|
| Efield |.00163948|.000530475|.0000759445|.00252516|
| Einteraction |-.00327896|-.00184815|-.000683795|-.00409026|
| Epol |0|.000393598|.000265953|-.000480030|
| Mtranslation |1.26089|1.26159|1.26216|1.26050|

Undressed classical core E=1.26251, Q=0. Negative q partner has bit-identical E and opposite Q (symmetry f unchanged, phi→-phi); q=0 stays zero even for alpha=5. Q changes by factor 5.08 between strongest screening and anti-screening at fixed q. Etotal changes ~.132%. Core profile relaxes, so Qbare varies slightly; strong suppression is mostly polarization, not core energy loss. Scalar field self-energy is positive; total dressing is negative due to interaction and polarization. Do not compare just gradient energy to mass. q=.5 cases confirm nonlinearity of profile and feedback. C=0 gives same qualitative outcome. Mtranslation=(2E2+4E4+2Efield)/3 is calculated from the common Lorentz-invariant kinetic action for a translated stationary profile, not velocity/impulse measurement. Difference from E given by Derrick; max relative R24 error .0121%, R12 .0729% (massless boundary tail). R12 vs R24 Q changes up to ~1.04%, massive <.04%; finite boundary effects disclosed. Full dynamics including acceleration, radiation and nonradial modes remains open.

Interpretation: the correct accounting is the dressed particle plus field and interaction stress. In this nonmetric scalar response, source remains a signed internal response rather than universal energy. Dressing cannot repair an exact odd-source symmetry on this branch. Under genuine metric coupling the total conserved stress would source, but that is not derived by this model. Gc27/28C's universality problem survives self-energy.
