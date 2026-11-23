"""EX-GC34D POST HOC diagnostic (not registered): straight-ray bending ratio for a ray along the face diagonal (1,1,0),
because rays along a lattice axis turned out to obey h_xx = -2n exactly for every k perpendicular to the ray.
Uses the stored full Gaussian solution. -> results/EX-GC34D/posthoc_diag.json"""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / 'src'))
import exgc34d as D
z = np.load(ROOT / 'local_raw/EX-GC34D/gauss_fields.npz'); nh, hh = z['nh'], z['hh']
N = D.NBOX; H = D.H
kz = 2 * np.pi * np.fft.rfftfreq(N[2], d=H)
out = {}
for lab, rows, phi_h in (('axis_x', (np.zeros(N[1], int), np.arange(N[1])), lambda hx: 0.5 * hx[0]),
                         ('diag_110', ((-2 * np.arange(N[1])) % N[0], np.arange(N[1])),
                          lambda hx: 0.25 * (hx[0] + hx[1] + 2 * hx[3]))):
    i, j = rows
    if lab == 'axis_x':
        sel = lambda F: F[0]                       # kx = 0 plane: all ky, kz
    else:
        sel = lambda F: F[i, j]                    # kx = -ky line, all kz
    n = sel(nh); hx = [sel(hh[c]) for c in range(6)]
    res = {}
    for b in (6.0, 8.0, 12.0):
        ph = 1j * kz * np.exp(1j * kz * b)
        a_l = 2 * np.real((-n * ph).sum()); a_f = 2 * np.real(((-n + phi_h(hx)) * ph).sum())
        res[str(b)] = dict(lapse=float(a_l), full=float(a_f), ratio=float(a_f / a_l))
    out[lab] = res
(ROOT / 'results/EX-GC34D/posthoc_diag.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))

# POST HOC 2: is the lattice solution exactly "each edge shrinks by the mean lapse of its ends" plus a lattice gauge
# mode, at every k? Least-squares fit of the 7 edge responses to l0_e psi (1 + e^{i q.d})/2 + G(q) c.
import mpmath as mp
mp.mp.dps = 50
terms, _ = D.build_terms()
rng = np.random.default_rng(3401)
qs = np.vstack([rng.uniform(-np.pi, np.pi, (20, 3)), [[0.8, 0.8, 0], [0.5, 0.5, 0.5], [0.05, 0.02, 0.01],
                                                       [np.pi, np.pi, 0.3]]])
l0 = np.array([D.H * np.linalg.norm(d) for d in D.ET]); rows = []
for q in qs:
    x = D.response(terms, q[None])[0]; G = D.gauge_modes(q[None])[0][1:]
    m = np.array([(1 + np.exp(1j * q @ d)) / 2 for d in D.ET])
    Bm = np.column_stack([l0 * m, G]); c = np.linalg.lstsq(Bm, x[1:], rcond=None)[0]
    rows.append(dict(q=q.tolist(), resid=float(np.linalg.norm(Bm @ c - x[1:]) / np.linalg.norm(x[1:])),
                     psi_over_n=[float((c[0] / x[0]).real), float((c[0] / x[0]).imag)]))
out['conformal_fit'] = dict(max_resid=max(r['resid'] for r in rows),
                            max_dev_psi_over_n=max(abs(r['psi_over_n'][0] + 1) + abs(r['psi_over_n'][1]) for r in rows),
                            rows=rows)
(ROOT / 'results/EX-GC34D/posthoc_diag.json').write_text(json.dumps(out, indent=1) + '\n')
print('conformal fit', out['conformal_fit']['max_resid'], out['conformal_fit']['max_dev_psi_over_n'])
