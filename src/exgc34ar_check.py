"""EX-GC34A-r separate recomputation (registered in spec/EX-GC34Ar_plan.md).

M1: packet position = peak of FFT-Gaussian-smoothed energy proxy of the final field (parabolic sub-cell fit), periodic
displacement from x0 = 0, v = dx/T; compared with the circular-centroid fit velocity in summary.json.
M2: late spectral direction (exgc34a_check.mean_vg, flat dispersion) vs the late-window circular-centroid direction.
Writes results/EX-GC34Ar/check.json and manifest_sha256.json.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[k] = '1'
import json, hashlib
from pathlib import Path
import numpy as np
import exgc34a_check as C

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC34Ar'; RAW = ROOT / 'local_raw/EX-GC34Ar'
H = 0.2; LO1 = np.array([-8.0, -8.0, -8.0]); L1 = 16.0


def proxy(f, v, fid, mode):
    if fid == 'F3':
        g = (f[0] + f[1]) / 2 if mode == 'sym' else (f[0] - f[1]) / 2
        gv = (v[0] + v[1]) / 2 if mode == 'sym' else (v[0] - v[1]) / 2
        return gv ** 2 + sum((np.roll(g, -1, a) - g) ** 2 for a in range(3)) / H ** 2
    if fid == 'F4':
        return (v[1:] ** 2).sum(0) + sum(((np.roll(f[1:], -1, a + 1) - f[1:]) ** 2).sum(0) for a in range(3)) / H ** 2
    return v ** 2 + sum((np.roll(f, -1, a) - f) ** 2 for a in range(3)) / H ** 2


def peak(e):
    n = e.shape
    K2 = sum(np.meshgrid(*[(2 * np.pi * np.fft.fftfreq(n[a], d=H)) ** 2 for a in range(3)], indexing='ij'))
    s = np.fft.ifftn(np.fft.fftn(e) * np.exp(-K2 / 2)).real          # Gaussian smoothing, width 1
    i = np.unravel_index(np.argmax(s), n); pos = []
    for a in range(3):
        im = list(i); ip = list(i); im[a] = (i[a] - 1) % n[a]; ip[a] = (i[a] + 1) % n[a]
        ym, y0, yp = s[tuple(im)], s[i], s[tuple(ip)]
        pos.append(LO1[a] + H * (i[a] + 0.5 * (ym - yp) / (ym - 2 * y0 + yp)))
    return np.array(pos)


def main():
    summ = json.loads((OUT / 'summary.json').read_text()); rows = []; man = {}
    for p in sorted(RAW.glob('M[12]_*.npz')):
        tag = p.stem; man[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        rec = json.loads((OUT / (tag + '.json')).read_text()); fid = rec['field']
        z = np.load(p); f, v = z['f'], z['v']
        modes = ['sym', 'anti'] if fid == 'F3' else [fid]
        for m in modes:
            if tag.startswith('M1'):
                d = peak(proxy(f, v, fid, m)); d = (d + L1 / 2) % L1 - L1 / 2      # periodic displacement from 0
                vp = d / rec['T']
                s = [r for r in summ['M1'] if r['tag'] == tag and r['mode'] == m][0]; vc = np.array(s['v'])
                rows.append(dict(tag=tag, mode=m, v_peak=vp.tolist(),
                                 rel_vs_centroid=float(np.linalg.norm(vp) / np.linalg.norm(vc) - 1),
                                 ang_vs_centroid=float(np.arccos(np.clip(vp @ vc / np.linalg.norm(vp) / np.linalg.norm(vc), -1, 1)))))
            else:
                ff, vv = ({'sym': ((f[0] + f[1]) / 2, (v[0] + v[1]) / 2), 'anti': ((f[0] - f[1]) / 2, (v[0] - v[1]) / 2)}[m]
                          if fid == 'F3' else ((f[1:], v[1:]) if fid == 'F4' else (f, v)))
                Ms, Vs, w = C.coeffs('B0', 1.0)
                if ff.ndim == 4:
                    vg = sum(C.mean_vg(ff[c], vv[c], Ms, Vs, w, C.MU[m]) * np.sum(ff[c] ** 2) for c in range(3)) \
                        / sum(np.sum(ff[c] ** 2) for c in range(3))
                else:
                    vg = C.mean_vg(ff, vv, Ms, Vs, w, C.MU[m])
                s = [r for r in summ['M2'] if r['tag'] == tag and r['mode'] == m][0]
                rows.append(dict(tag=tag, mode=m, theta_late_spectral=float(np.arctan2(vg[1], vg[0])),
                                 theta_late_centroid=s['theta_late']))
    (OUT / 'check.json').write_text(json.dumps(rows, indent=1) + '\n')
    (OUT / 'manifest_sha256.json').write_text(json.dumps(man, indent=1) + '\n')


if __name__ == '__main__':
    main()
