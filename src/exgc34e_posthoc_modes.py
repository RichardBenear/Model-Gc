"""EX-GC34E POST HOC (not registered): vacuum mode search that tolerates an exactly null re-slicing direction.
The registered E3 method (linearised quadratic eigenproblem after dropping axis-edge columns/rows) breaks when the
pencil is singular, which happens exactly where the lattice keeps the re-slicing symmetry (q along an axis).
Here: project the kernel off the displacement gauge at each w, scan w, and find the frequencies where an extra
singular value dips to zero.   -> results/EX-GC34E/posthoc_modes.json"""
import json, sys
from pathlib import Path
import numpy as np
import scipy.linalg as sl
from scipy.optimize import minimize_scalar
import mpmath as mp
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / 'src'))
import exgc34e as E
mp.mp.dps = 50
st = E.tet_static()


def svals(q, w):
    Mm = E.full_kernel(st, q[None], w)[0]; G = E.gauge_modes(q[None], w)[0]
    Qc = sl.null_space(G.conj().T)
    s = np.linalg.svd(Qc.conj().T @ Mm @ Qc, compute_uv=False)
    return s[::-1] / s[0]          # ascending, relative to largest


out = []
for kh in (0.05, 0.1, 0.2):
    for dn, dv in (('100', (1, 0, 0)), ('110', (1, 1, 0)), ('111', (1, 1, 1))):
        q = kh * np.array(dv, float) / np.linalg.norm(dv); k = kh / E.H
        ws = np.linspace(0.02, 3.0, 1500) * k
        S = np.array([svals(q, w) for w in ws])
        base = S[:, 0]                                  # re-slicing direction (exact or near) at every w
        dips = []
        for j in (1, 2, 3):
            c = np.log(S[:, j])
            for i in range(1, len(ws) - 1):
                if c[i] < c[i - 1] and c[i] <= c[i + 1] and c[i] < np.median(c) - 4:
                    r = minimize_scalar(lambda w: np.log(svals(q, w)[j]), bracket=(ws[i - 1], ws[i], ws[i + 1]),
                                        tol=1e-12)
                    dips.append(dict(col=j, w=float(r.x), speed=float(r.x / k), depth=float(np.exp(r.fun))))
        # fine scan near w = k: every deep zero of the second singular value is one mode
        wf = np.linspace(0.985, 1.005, 4001) * k
        c1 = np.log(np.array([svals(q, w)[1] for w in wf]))
        fine = []
        for i in range(1, len(wf) - 1):
            if c1[i] < c1[i - 1] and c1[i] <= c1[i + 1]:
                r = minimize_scalar(lambda w: np.log(svals(q, w)[1]), bracket=(wf[i - 1], wf[i], wf[i + 1]), tol=1e-14)
                sv = svals(q, r.x)
                fine.append(dict(speed=float(r.x / k), sv1=float(sv[1]), sv2=float(sv[2])))
        speeds = [f['speed'] for f in fine if f['sv1'] < 1e-11]
        out.append(dict(kh=kh, dir=dn, reslicing_sv_median=float(np.median(base)),
                        reslicing_sv_max=float(base.max()), dips=dips, fine=fine, speeds=speeds, c_light=float(E.light_speed(q))))
        print(kh, dn, 'reslice sv ~%.1e' % np.median(base), 'dips', [(d['col'], round(d['speed'], 6), '%.0e' % d['depth'])
                                                                     for d in dips], 'c_light %.6f' % E.light_speed(q), 'modes', [(round(f['speed'], 6), '%.0e' % f['sv1'], '%.0e' % f['sv2']) for f in fine])
(ROOT / 'results/EX-GC34E/posthoc_modes.json').write_text(json.dumps(out, indent=1) + '\n')

# POST HOC 2: does the near-null re-slicing direction (q off-axis) hide a slow extra mode? Scan the smallest
# singular value from w = 1e-4 k to 3 k on a log grid and report its minimum relative to its median.
extra = []
for kh in (0.05, 0.1, 0.2):
    for dn, dv in (('110', (1, 1, 0)), ('111', (1, 1, 1)), ('gen', (0.8, 0.5, 0.33))):
        q = kh * np.array(dv, float) / np.linalg.norm(dv); k = kh / E.H
        ws = np.geomspace(1e-4, 3.0, 2000) * k
        s0 = np.array([svals(q, w)[0] for w in ws])
        i = int(np.argmin(s0))
        extra.append(dict(kh=kh, dir=dn, s0_min=float(s0[i]), w_at_min_over_k=float(ws[i] / k),
                          s0_median=float(np.median(s0)), s0_at_w0=float(svals(q, 0.0)[0])))
        print('extra', kh, dn, '%.2e at w/k=%.4f, median %.2e, static %.2e' % (s0[i], ws[i] / k, np.median(s0), svals(q, 0.0)[0]))
d = json.loads((ROOT / 'results/EX-GC34E/posthoc_modes.json').read_text())
(ROOT / 'results/EX-GC34E/posthoc_modes.json').write_text(json.dumps(dict(modes=d, reslicing_scan=extra), indent=1) + '\n')
