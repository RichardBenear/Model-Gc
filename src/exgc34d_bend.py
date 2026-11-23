"""EX-GC34D part D4: 34A-r light-bending machinery on the Regge-solved background vs lapse-only (same n).

Fields (exgc34a solver, F1 massless) on the 34A-r M2 grid; rays (exgc34a_rays Hamiltonian) with the same metric via
cubic-spline interpolation of the site fields. Readout identical to src/exgc34ar_eval.py (circular centroid,
early t 2-8 vs late t >= 28 direction).

  python src/exgc34d.py first (writes local_raw/EX-GC34D/gauss_fields.npz)
  python src/exgc34d_bend.py      -> results/EX-GC34D/bend.json
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = '1'
import json, sys, time
from pathlib import Path
from types import SimpleNamespace
from multiprocessing import Pool
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
os.environ['EXGC34A_STUDY'] = 'EX-GC34D'
import exgc34a as S
import exgc34a_rays as R
import exgc34d as D

OUT = ROOT / 'results/EX-GC34D'; RAW = ROOT / 'local_raw/EX-GC34D'
H = 0.2; DT = 0.25 * H; SIG = 2.0
LO = (-24.0, -6.0, -6.0); NM2 = (240, 100, 60)
BOX = (LO, (48.0, 20.0, 12.0))
X0 = np.array([-18.0, 6.0, 0.0]); K0 = np.array([3.0, 0.0, 0.0]); T = 36.0
COMP = ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))


def fields_at(nh, hh, s):
    """n, h (3,3) on the M2 grid at site offset s (lattice units), from origin-centred rfft coefficients."""
    N = D.NBOX
    q = np.meshgrid(2 * np.pi * np.fft.fftfreq(N[0]), 2 * np.pi * np.fft.fftfreq(N[1]),
                    2 * np.pi * np.fft.rfftfreq(N[2]), indexing='ij')
    ph = np.exp(1j * (q[0] * s[0] + q[1] * s[1] + q[2] * s[2]))
    idx = [(np.arange(NM2[a]) + int(round(LO[a] / H))) % N[a] for a in range(3)]
    cut = lambda F: np.fft.irfftn(F * ph, s=N)[np.ix_(*idx)]
    n = cut(nh); hm = np.zeros(NM2 + (3, 3))
    for c, (a, b) in enumerate(COMP):
        hm[..., a, b] = hm[..., b, a] = cut(hh[c])
    return n, hm


def A_of(n, hm, spatial):
    N = 1 + n
    g = np.eye(3) + (hm if spatial else 0 * hm)
    sg = np.sqrt(np.linalg.det(g))
    return N, sg, (N * sg)[..., None, None] * np.linalg.inv(g)


def build_backgrounds():
    z = np.load(RAW / 'gauss_fields.npz'); nh, hh = z['nh'], z['hh']
    G = S.Grid(NM2, LO, H)
    site = fields_at(nh, hh, (0, 0, 0))
    links = [fields_at(nh, hh, 0.5 * S.DIRS[d]) for d in range(9)]
    bgs = {}
    for name, spatial in (('regge', True), ('lapse', False)):
        N, sg, _ = A_of(*site, spatial)
        w = np.stack([S.split_weights(A_of(*links[d], spatial)[2], 12.0)[d] for d in range(9)])
        bgs[name] = SimpleNamespace(G=G, Ms=sg / N, Vs=N * sg, w=w)
    np.savez(RAW / 'm2_site_fields.npz', n=site[0], h=site[1])
    return bgs


class InterpBg(R.Bg):
    """Ray background from the solved site fields (cubic spline, periodic on the M2 box)."""
    def __init__(self, spatial):
        super().__init__('flat', dt_lf=DT)
        z = np.load(RAW / 'm2_site_fields.npz')
        self.spatial = spatial
        self.c = [ndimage.spline_filter(z['n'], order=3, mode='grid-wrap')]
        self.c += [ndimage.spline_filter(z['h'][..., a, b], order=3, mode='grid-wrap') for a, b in COMP]

    def coeffs(self, x):
        shp = x.shape[:-1]
        idx = [((x[..., a] - LO[a]) / H).ravel() for a in range(3)]
        v = [ndimage.map_coordinates(c, idx, order=3, mode='grid-wrap', prefilter=False).reshape(shp)
             for c in self.c]
        n = v[0]; hm = np.zeros(shp + (3, 3))
        if self.spatial:
            for c, (a, b) in enumerate(COMP):
                hm[..., a, b] = hm[..., b, a] = v[1 + c]
        N, sg, A = A_of(n, hm, True)
        return sg / N, N * sg, R.weights_from_A(A, 12.0)


def ray_job(name):
    bg = InterpBg(name == 'regge')
    rec = json.loads((OUT / f'D4_{name}_F1.json').read_text())
    t = np.array(rec['hist']['t'])
    c = R.rays(bg, X0, K0, SIG, 0.0, H, T, DT, nq=4, times=t, box=BOX)
    return name, c.tolist()


def unwrap(c):
    lo, Lb = np.asarray(BOX[0]), np.asarray(BOX[1])
    return lo + Lb * np.unwrap(2 * np.pi * (np.asarray(c) - lo) / Lb, axis=0) / (2 * np.pi)


def alpha(t, cc):
    ang = lambda v: float(np.arctan2(v[1], v[0]))
    ve = R.fit_v(t[(t >= 2) & (t <= 8)], cc[(t >= 2) & (t <= 8)], 0)
    vl = R.fit_v(t[t >= 28], cc[t >= 28], 0)
    return ang(vl) - ang(ve)


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True); RAW.mkdir(parents=True, exist_ok=True)
    bgs = build_backgrounds()
    for name, bg in bgs.items():
        tag = f'D4_{name}_F1'
        if (OUT / (tag + '.json')).exists():
            continue
        hist = S.run('scalar', bg, X0, K0, SIG, T, DT, 10, save_final=RAW / (tag + '.npz'))
        S.save(tag, dict(tag=tag, bg=name, field='F1', h=H, dt=DT, sigma=SIG, T=T, hist=hist,
                         platform=S.platform_info()))
        print(tag, 'wall %.1f' % hist['wall_s'], flush=True)
    with Pool(2) as P:
        rays = dict(P.map(ray_job, ['regge', 'lapse']))
    out = {}
    for name in ('regge', 'lapse'):
        rec = json.loads((OUT / f'D4_{name}_F1.json').read_text()); hist = rec['hist']
        t = np.array(hist['t'])
        out[name] = dict(alpha=alpha(t, unwrap(hist['cc'])), alpha_ray=alpha(t, unwrap(rays[name])),
                         E_drift=float(hist['c_E'][-1] / hist['c_E'][0] - 1))
    tol = 0.05 * abs(out['lapse']['alpha_ray'])
    for name in out:
        out[name]['pass_P6'] = bool(abs(out[name]['alpha'] - out[name]['alpha_ray']) <= tol)
    out['tol_P6'] = tol
    out['ratio_field'] = out['regge']['alpha'] / out['lapse']['alpha']
    out['ratio_ray'] = out['regge']['alpha_ray'] / out['lapse']['alpha_ray']
    out['pass_P5'] = bool(abs(out['ratio_field'] / 2 - 1) <= 0.03 and abs(out['ratio_ray'] / 2 - 1) <= 0.03)
    out['wall_s'] = time.time() - t0
    (OUT / 'bend.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
