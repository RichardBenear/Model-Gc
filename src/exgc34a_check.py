"""EX-GC34A separate recomputation: propagation direction and speed from the final fields' spectra.

For each saved final state (local_raw, not committed) form the positive-frequency amplitude a(q) = f^(q) + i v^(q)/w(q)
with the lattice dispersion w(q) of the *uniform* background (M1: exact; M2: flat dispersion, since at the end the packet
is at x ~ 15 where |Phi| < 0.004), weight by w^2 |a|^2 and average the group velocity dw/dq (finite difference).
Compared with the centroid-fit velocity (M1) and late-window direction (M2) in summary.json. Writes check.json and
manifest_sha256.json. Shared dependency: the solver's fields and the background weight split.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[k] = '1'
import json, hashlib
from pathlib import Path
import numpy as np
import exgc34a_rays as R

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC34A'; RAW = ROOT / 'local_raw/EX-GC34A'
H = 0.2; DT = 0.25 * H
MU = {'F1': 0.0, 'F2': 1.0, 'F4': 0.0, 'F5': 0.0, 'sym': 0.0, 'anti': np.sqrt(2 * 0.6 / 0.6)}


def gam():
    import exgc34a
    return exgc34a.gam_B1()


def coeffs(bgname, scale):
    if bgname in ('B1', 'B1p'):
        g = gam(); sg = np.sqrt(np.linalg.det(g))
        w = R.weights_from_A(np.linalg.inv(g) * sg * scale, 12.0 if bgname == 'B1' else 16.0)
        return sg, sg, w
    return 1.0, 1.0, R.weights_from_A(np.eye(3) * scale, 12.0)


def omega(K, Ms, Vs, w, mu):
    q = (w * 4 * np.sin(H * (K @ R.DIRS.T) / 2) ** 2).sum(-1) / H ** 2
    om = np.sqrt((q + Vs * mu ** 2) / Ms)
    return 2 / DT * np.arcsin(om * DT / 2)


def mean_vg(f, v, Ms, Vs, w, mu):
    n = f.shape
    K = np.stack(np.meshgrid(*[2 * np.pi * np.fft.fftfreq(n[a], d=H) for a in range(3)], indexing='ij'), -1)
    om = omega(K, Ms, Vs, w, mu)
    a = np.fft.fftn(f) + 1j * np.fft.fftn(v) / np.maximum(om, 1e-12)
    wt = om ** 2 * abs(a) ** 2
    # keep the forward branch only (energy is concentrated near +k0); backward residue reported
    vg = np.zeros(K.shape); e = 1e-5
    for i in range(3):
        d = np.zeros(3); d[i] = e
        vg[..., i] = (omega(K + d, Ms, Vs, w, mu) - omega(K - d, Ms, Vs, w, mu)) / (2 * e)
    return (wt[..., None] * vg).sum((0, 1, 2)) / wt.sum()


def main():
    summ = json.loads((OUT / 'summary.json').read_text())
    rows = []; man = {}
    for p in sorted(RAW.glob('M[12]_*.npz')):
        tag = p.stem; man[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        rec = json.loads((OUT / (tag + '.json')).read_text()); fid = rec['field']; bg = rec['bg']
        z = np.load(p); f, v = z['f'], z['v']
        if fid == 'F3':
            comps = {'sym': ((f[0] + f[1]) / 2, (v[0] + v[1]) / 2), 'anti': ((f[0] - f[1]) / 2, (v[0] - v[1]) / 2)}
        elif fid == 'F4':
            comps = {'F4': (f[1:], v[1:])}
        else:
            comps = {fid: (f, v)}
        for m, (ff, vv) in comps.items():
            Ms, Vs, w = coeffs(bg if tag.startswith('M1') else 'B0', 1.10 if fid == 'F5' else 1.0)
            if ff.ndim == 4:
                vg = sum(mean_vg(ff[c], vv[c], Ms, Vs, w, MU[m]) * np.sum(ff[c] ** 2) for c in range(ff.shape[0])) \
                    / sum(np.sum(ff[c] ** 2) for c in range(ff.shape[0]))
            else:
                vg = mean_vg(ff, vv, Ms, Vs, w, MU[m])
            s = [r for r in summ['M1' if tag.startswith('M1') else 'M2'] if r['tag'] == tag and r['mode'] == m][0]
            row = dict(tag=tag, mode=m, vg_spectral=vg.tolist())
            if tag.startswith('M1'):
                vc = np.array(s['v']); row['rel_vs_centroid'] = float(np.linalg.norm(vg) / np.linalg.norm(vc) - 1)
                row['ang_vs_centroid'] = float(np.arccos(np.clip(vg @ vc / np.linalg.norm(vg) / np.linalg.norm(vc), -1, 1)))
            else:
                row['theta_late_spectral'] = float(np.arctan2(vg[1], vg[0]))
            rows.append(row)
    (OUT / 'check.json').write_text(json.dumps(rows, indent=1) + '\n')
    (OUT / 'manifest_sha256.json').write_text(json.dumps(man, indent=1) + '\n')
    m1 = [r for r in rows if r['tag'].startswith('M1')]
    print('M1 max |rel| %.2e  max ang %.2e' % (max(abs(r['rel_vs_centroid']) for r in m1),
                                               max(r['ang_vs_centroid'] for r in m1)))


if __name__ == '__main__':
    main()
