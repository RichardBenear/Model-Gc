# EX-GC35N — absolute lens deflection at h=.08

2026-11-24 | EXPLORATORY | ChatGPT | [Plan](../../spec/EX-GC35N_plan.md). PI closure pending.

**The full-to-lapse bend ratio remains 2.007. The excess after the leading dispersion correction did not decrease materially.** The same static mass, beam, clock reconstruction and finite-path comparator as GC35K were used at packet spacing .08 and proper step .013333. The full-metric bend was 3.039°, against 2.903° in the leading asymptotic weak-field expression. The dedicated .08 runs were a lapse/full pair. Raw group angles are used below, with the prior flat direction within 10⁻⁶ degrees of zero.

| Spacing | Metric | Measured group bend | Finite-path continuum | Finite path × sec²(kh/2) | Residual excess |
|---:|---|---:|---:|---:|---:|
| 0.12 | lapse-only | 1.560° | 1.437° | 1.523° | 2.44% |
| 0.12 | full | 3.135° | 2.834° | 3.004° | 4.37% |
| 0.08 | lapse-only | 1.514° | 1.439° | 1.477° | 2.53% |
| 0.08 | full | 3.039° | 2.838° | 2.911° | 4.38% |

The finite-path comparator uses the measured longitudinal span in the same periodic metric. The full residual after the leading carrier correction stays around 4.4%. The data are consistent with additional finite-width, path or lattice contributions. Two samples do not establish an asymptotic floor. The wave norm ratios were .99846 and .99830, with 1,938,354,455 events per run. The pair used 781 CPU seconds.

[Data](results.json) · [Driver](../../src/exgc35k.py). Reproduce with `--h .08 --step .0133333333333333 --mode lapse|full`.
