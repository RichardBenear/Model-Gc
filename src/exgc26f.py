"""EX-GC26F (exploratory): two B = 1 knots attracting through their self-generated Nordström field.

Base: the closed EX-GC26C coupled solver (common action) and the EX-GC26D inputs (closed 26C end state, hash-checked).
Two copies of the 26C knot at x = -d/2, +d/2 (same orientation: product ansatz), each with its static field; the field
wall is held at the exact superposed static value -A/|x - x1| - A/|x - x2| (no taper). Control: same matter start,
G_N = 0, Phi = 0. Gravity = coupled minus control separation. Massless Skyrme reference model, not Model Gc.

  python src/exgc26f.py run --d 6 --GN 0.01      # coupled
  python src/exgc26f.py run --d 6 --GN 0         # control
  python src/exgc26f.py analyse
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import argparse, hashlib, json, platform, time
from pathlib import Path
import numpy as np
import scipy
from scipy.ndimage import map_coordinates, spline_filter
import exgc26c as c
import exgc26d as d26
import exgc26e1 as e1

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC26F'
RAW = ROOT / 'local_raw/EX-GC26F'
INP = ROOT / 'local_raw/EX-GC26D'          # the 26D input arrays are reused unchanged
T_RUN = 10.0
MARGIN = 7.0                               # matter half-width = d/2 + 7 (the single-knot box beyond each knot)
FIELD_EXTRA = 4.0                          # field half-width = matter half-width + 4


def save(tag, r):
    OUT.mkdir(parents=True, exist_ok=True); p = OUT / (tag + '.json'); q = p.with_suffix('.tmp')
    q.write_text(json.dumps(r, indent=1, allow_nan=False) + '\n'); q.replace(p)


def inputs():
    man = json.loads((d26.OUT / 'inputs.json').read_text())
    def inp(name):
        path = INP / name; assert hashlib.sha256(path.read_bytes()).hexdigest() == man[name], name; return np.load(path)
    return man, inp('c26_p.npy'), inp('c26_Fcrop.npy')


def qmul(a, b):
    a0, av, b0, bv = a[0], a[1:], b[0], b[1:]
    return np.concatenate([[a0 * b0 - np.sum(av * bv, 0)], a0 * bv + b0 * av + np.cross(av, bv, axis=0)])


def place(G, G1, p1, cx):
    """The single-knot state p1 (grid G1) shifted to x = cx on grid G; vacuum outside its box. cx is grid-aligned."""
    X = np.array(G.X, float); X[0] -= cx; idx = X / G1.h + (G1.n - 1) / 2
    q = np.stack([map_coordinates(p1[k], idx, order=1, mode='constant', cval=(1. if k == 0 else 0.)) for k in range(4)])
    return q / np.sqrt(np.sum(q * q, 0))


def rest_fn(Fd1, Frest, A):
    """26D rest field, blended over r = 9.5-10.5 into -A/r, without the far taper."""
    coeff = spline_filter(Frest, order=3, mode='constant')
    def fn(X):
        r = np.sqrt(np.sum(X ** 2, axis=0)); idx = X / Fd1.h + (Fd1.n - 1) / 2
        inner = map_coordinates(coeff, idx, order=3, mode='nearest', prefilter=False)
        w = np.clip((10.5 - r) / 1.0, 0, 1); w = w * w * (3 - 2 * w)
        return w * inner + (1 - w) * (-A / np.maximum(r, 1e-9))
    return fn


def halves(G, dens):
    """x-centroids of a density on the x < 0 and x > 0 halves (the mid-plane x = 0 is split equally)."""
    X = G.X[0]; out = []
    for sgn in (-1, 1):
        w = np.where(sgn * X > 0, 1., np.where(X == 0, .5, 0.)) * dens
        out.append(float(np.sum(w * X) / np.sum(w)))
    return out


def setup(dsep, GN):
    e1.configure(0.2); e1.load(); h = e1.H
    assert abs(dsep / 2 / h - round(dsep / 2 / h)) < 1e-9, 'd/2 must be grid-aligned'
    L = dsep / 2 + MARGIN; G = e1.ref.Grid(h, L); Fd = c.Field(G, L + FIELD_EXTRA)
    man, p1, Fc = inputs(); G1 = e1.ref.Grid(h, e1.L)
    p = G.fix(qmul(place(G, G1, p1, -dsep / 2), place(G, G1, p1, dsep / 2)))
    G26 = e1.ref.Grid(h, e1.L); Fd1 = c.Field(G26, c.L_PHI); k = man['crop_half_cells']; cc = Fd1.n // 2
    Frest = np.zeros(Fd1.r.shape); Frest[cc - k:cc + k + 1, cc - k:cc + k + 1, cc - k:cc + k + 1] = Fc
    A = d26.fit_A(Fd1, Frest)
    if GN > 0:
        fn = rest_fn(Fd1, Frest, A); X = np.array(np.meshgrid(Fd.ax, Fd.ax, Fd.ax, indexing='ij'))
        F = np.zeros(Fd.r.shape)
        for cx in (-dsep / 2, dsep / 2):
            Y = X.copy(); Y[0] -= cx; F += fn(Y)
    else:
        F = np.zeros(Fd.r.shape)
    Pi = np.zeros_like(F)
    return G, Fd, p, F, Pi, A, man


def run(dsep, GN):
    t0 = time.process_time(); G, Fd, p, F, Pi, A, man = setup(dsep, GN); tag = f'd{dsep:g}_GN{GN:g}'
    v = np.zeros_like(p); s = (p, v, F, Pi); dt = .05 * G.h; steps = round(T_RUN / dt); every = round(.25 / dt); rows = []
    d0 = e1.a2.geometric(p); fac = G.h ** 3 / e1.NORM
    for i in range(steps + 1):
        t = i * dt
        if i % every == 0 or i == steps:
            p, v, F, Pi = s; W = np.exp(2 * F[Fd.sub]); Hd, P, Bd, _ = e1.fields(G, p, v, W)
            rows.append(dict(t=t, H_m=float(Hd.sum() * fac), H_phi=Fd.energy(F, Pi, GN) if GN > 0 else 0., P_m=(P.sum((1, 2, 3)) * fac).tolist(),
                             P_phi=d26.field_momentum(Fd, F, Pi, GN), centroids=halves(G, Hd), bcentroids=halves(G, Bd),
                             phi_at_knots=[float(F[Fd.sub][np.argmin(np.abs(G.ax - x)), G.n // 2, G.n // 2]) for x in (-dsep / 2, dsep / 2)]))
            if (i // every) % 8 == 0: print(tag, round(t, 2), round(time.process_time() - t0), rows[-1]['centroids'], flush=True)
        if i < steps:
            s = c.step(G, Fd, s, dt, GN) if GN > 0 else d26.step0(G, s, dt)
    p, v, F, Pi = s
    RAW.mkdir(parents=True, exist_ok=True); end = RAW / f'{tag}_end.npz'; np.savez(end, p=p, v=v, F=F, Pi=Pi)
    r = dict(d=dsep, GN=GN, L=float(G.ax[-1]), L_phi=float(Fd.ax[-1]), T=steps * dt, dt=dt, A_single=A, source_sha256=man['c26_p.npy'],
             field_sha256=man['c26_Fcrop.npy'], rows=rows, degree_initial=d0, degree_final=e1.a2.geometric(p),
             end_sha256=hashlib.sha256(end.read_bytes()).hexdigest(), cpu_s=time.process_time() - t0,
             platform=f'Python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__} {platform.platform()}')
    save(tag, r); print('DONE', tag, r['cpu_s'], flush=True)


def accel(t, D, t_min=0.):
    sel = t >= t_min; c2 = np.polyfit(t[sel], D[sel], 2); return float(2 * c2[0])


def analyse():
    c26 = json.loads((ROOT / 'results/EX-GC26C/B1_GN0.01.json').read_text()); ratio_p = 2 * c26['E2'] / c26['E']
    res = dict(passive_ratio_2E2_over_E=ratio_p)
    for dsep in (6.0, 9.2):
        a, b = OUT / f'd{dsep:g}_GN0.01.json', OUT / f'd{dsep:g}_GN0.json'
        if not (a.exists() and b.exists()): continue
        ra, rb = json.loads(a.read_text()), json.loads(b.read_text()); t = np.array([x['t'] for x in ra['rows']])
        out = dict()
        for key in ('centroids', 'bcentroids'):
            Dc = np.array([x[key][1] - x[key][0] for x in ra['rows']]); Dn = np.array([x[key][1] - x[key][0] for x in rb['rows']])
            pred = -2 * ratio_p * ra['A_single'] / Dn ** 2           # relative acceleration, Newtonian with lattice passive charge
            aG = accel(t, Dc - Dn); aG1 = accel(t, Dc - Dn, 1.0)
            out[key] = dict(a_rel_gravity=aG, a_rel_gravity_from_t1=aG1, predicted=float(pred.mean()), ratio=aG / float(pred.mean()),
                            ratio_from_t1=aG1 / float(pred.mean()), a_rel_control=accel(t, Dn), D0=float(Dn[0]), D_control_end=float(Dn[-1]),
                            gravity_displacement_end=float((Dc - Dn)[-1] - (Dc - Dn)[0]))
        H = np.array([x['H_m'] + x['H_phi'] for x in ra['rows']]); Pt = np.array([np.array(x['P_m']) + [x['P_phi'], 0, 0] for x in ra['rows']])
        out.update(energy_drift=float(np.max(np.abs(H / H[0] - 1))), total_P_max_over_E=float(np.max(np.abs(Pt)) / H[0]),
                   degree=[ra['degree_initial']['degrees'], ra['degree_final']['degrees']], cpu_s=[ra['cpu_s'], rb['cpu_s']])
        res[f'd{dsep:g}'] = out
    if 'd6' in res and 'd9.2' in res:
        res['inverse_square_check'] = {k: dict(measured=res['d6'][k]['a_rel_gravity'] / res['d9.2'][k]['a_rel_gravity'],
                                               predicted=res['d6'][k]['predicted'] / res['d9.2'][k]['predicted']) for k in ('centroids', 'bcentroids')}
    (OUT / 'summary.json').write_text(json.dumps(res, indent=1) + '\n'); print(json.dumps(res, indent=1))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('phase', choices=['run', 'analyse']); ap.add_argument('--d', type=float, default=6.)
    ap.add_argument('--GN', type=float, default=.01); ap.add_argument('--T', type=float, default=None); a = ap.parse_args()
    global T_RUN
    if a.T: T_RUN = a.T
    if a.phase == 'run': run(a.d, a.GN)
    else: analyse()


if __name__ == '__main__':
    main()
