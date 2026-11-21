"""EX-GC26D (exploratory): boosted Skyrmion together with its own Nordström field; is there drag?

Uses the closed EX-GC26C coupled solver (common action) and the 26B Lorentz-boost initial data. Massless Skyrme
reference model, not Model Gc; the construction is ours.

  inputs: local_raw/EX-GC26D/{c26_p.npy, c26_Fcrop.npy, rest_B1_h0.2_L7.npy}, hashes in results/EX-GC26D/inputs.json
  python src/exgc26d.py run --u 0.3 --GN 0.01     # coupled (needs prep)
  python src/exgc26d.py run --u 0.3 --GN 0        # uncoupled control (26E1 rest state, Phi = 0)
  python src/exgc26d.py analyse
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
import exgc26e1 as e1
import exgc26b_mass as bm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC26D'
RAW = ROOT / 'local_raw/EX-GC26D'
T_RUN = 8.0


def save(tag, r):
    OUT.mkdir(parents=True, exist_ok=True); p = OUT / (tag + '.json'); q = p.with_suffix('.tmp')
    q.write_text(json.dumps(r, indent=1, allow_nan=False) + '\n'); q.replace(p)


def prep():
    c.OUT = OUT; c.RAW = RAW
    c.run(1, 0.01, 12.0)                      # writes results/EX-GC26D/B1_GN0.01.json and local_raw/EX-GC26D/B1_GN0.01_end.npz
    (OUT / 'B1_GN0.01.json').replace(OUT / 'prep_B1.json')


def fit_A(Fd, F):
    prof = c.shells(Fd, F, rmax=11.0); rr = np.array([x['r'] for x in prof]); ph = np.array([x['mean'] for x in prof])
    sel = (rr >= 3) & (rr <= 9.5)
    return -float(np.sum(ph[sel] / rr[sel]) / np.sum(1 / rr[sel] ** 2))


def rest_field_fn(Fd, F, A):
    coeff = spline_filter(F, order=3, mode='constant')
    def fn(X):            # X: (3, ...) physical coordinates in the rest frame
        r = np.sqrt(np.sum(X ** 2, axis=0)); idx = X / Fd.h + (Fd.n - 1) / 2
        inner = map_coordinates(coeff, idx, order=3, mode='nearest', prefilter=False)
        w = np.clip((10.5 - r) / 1.0, 0, 1); w = w * w * (3 - 2 * w)          # 1 for r <= 9.5, 0 for r >= 10.5
        tp = np.clip((16.0 - r) / 4.0, 0, 1); tp = tp * tp * (3 - 2 * tp)   # far taper 12 -> 16 so Phi meets the wall value 0 smoothly
        return w * inner + (1 - w) * tp * (-A / np.maximum(r, 1e-9))
    return fn


def boost_field(Fd, fn, u, x0):
    g = 1 / np.sqrt(1 - u * u); X = np.array(np.meshgrid(Fd.ax, Fd.ax, Fd.ax, indexing='ij'))
    Xr = X.copy(); Xr[0] = g * (X[0] - x0); eps = 1e-4
    F = fn(Xr); dP = np.zeros_like(Xr); dP[0] = eps
    der = (fn(Xr + dP) - fn(Xr - dP)) / (2 * eps)
    Pi = -g * u * der
    for A_ in (F, Pi):
        A_[0, :, :] = A_[-1, :, :] = A_[:, 0, :] = A_[:, -1, :] = A_[:, :, 0] = A_[:, :, -1] = 0.
    return F, Pi


def field_momentum(Fd, F, Pi, GN):
    if GN == 0: return 0.
    dF = np.zeros_like(F); dF[1:-1] = (F[2:] - F[:-2]) / (2 * Fd.h)
    return float(-np.sum(Pi * dF) * Fd.h ** 3 / (4 * np.pi * GN))


def run(u, GN):
    t0 = time.process_time(); e1.configure(0.2); e1.load(); G = e1.ref.Grid(e1.H, e1.L); Fd = c.Field(G, c.L_PHI)
    x0 = -u * T_RUN / 2; n = np.array([1., 0., 0.]); tag = f'u{u:g}_GN{GN:g}'
    man = json.loads((OUT / 'inputs.json').read_text())
    def inp(name):
        path = RAW / name; assert hashlib.sha256(path.read_bytes()).hexdigest() == man[name], name; return np.load(path)
    if GN > 0:
        # closed 26C end state (cloud), shipped as the matter field p and the field Phi cropped to +-12 (see inputs.json)
        p_rest = inp('c26_p.npy'); Fc = inp('c26_Fcrop.npy'); k = man['crop_half_cells']; cc = Fd.n // 2
        Frest = np.zeros(Fd.r.shape); Frest[cc - k:cc + k + 1, cc - k:cc + k + 1, cc - k:cc + k + 1] = Fc
        resid = dict(v_max_26C_end=man['v_max'])
        A = fit_A(Fd, Frest); fn = rest_field_fn(Fd, Frest, A); F, Pi = boost_field(Fd, fn, u, x0); src = man['c26_p.npy']
    else:
        p_rest = inp('rest_B1_h0.2_L7.npy'); F = np.zeros(Fd.r.shape); Pi = np.zeros_like(F); A = 0.; resid = {}; src = man['rest_B1_h0.2_L7.npy']
    p, v = bm.init_boost(G, p_rest, u, n, t=x0 / u)
    s = (p, v, F, Pi); dt = .05 * G.h; steps = round(T_RUN / dt); every = round(.25 / dt); rows = []; d0 = e1.a2.geometric(p)
    fac = G.h ** 3 / e1.NORM
    for i in range(steps + 1):
        t = i * dt
        if i % every == 0 or i == steps:
            p, v, F, Pi = s; W = np.exp(2 * F[Fd.sub]); Hd, P, Bd, _ = e1.fields(G, p, v, W)
            Hphi = Fd.energy(F, Pi, GN) if GN > 0 else 0.
            rows.append(dict(t=t, H_m=float(Hd.sum() * fac), H_phi=Hphi, P_m=(P.sum((1, 2, 3)) * fac).tolist(), P_phi=field_momentum(Fd, F, Pi, GN),
                             centroid=[float(np.sum(Hd * x) / Hd.sum()) for x in G.X], bcentroid=[float(np.sum(Bd * x) / Bd.sum()) for x in G.X]))
            if (i // every) % 8 == 0: print(tag, round(t, 2), round(time.process_time() - t0), rows[-1]['centroid'][0], flush=True)
        if i < steps:
            s = c.step(G, Fd, s, dt, GN) if GN > 0 else step0(G, s, dt)
    p, v, F, Pi = s
    # field co-motion at the end, on lines through the knot's energy centroid
    X = rows[-1]['centroid'][0]; comov = []
    if GN > 0:
        coeff = spline_filter(F, order=3, mode='constant')
        def at(x, y, z): return float(map_coordinates(coeff, np.array([[x], [y], [z]]) / Fd.h + (Fd.n - 1) / 2, order=3, prefilter=False)[0])
        for d in (2., 3., 4., 5.):
            fore, aft, side = at(X + d, 0, 0), at(X - d, 0, 0), at(X, d, 0)
            comov.append(dict(d=d, fore=fore, aft=aft, side=side, fore_over_aft=fore / aft, axial_over_side=(fore + aft) / 2 / side))
    RAW.mkdir(parents=True, exist_ok=True); end = RAW / f'{tag}_end.npz'; np.savez(end, p=p, v=v, F=F, Pi=Pi)
    r = dict(u=u, GN=GN, x0=x0, T=steps * dt, dt=dt, A_rest=A, rest_residual=resid, source_sha256=src, rows=rows, comotion=comov,
             degree_initial=d0, degree_final=e1.a2.geometric(p), end_sha256=hashlib.sha256(end.read_bytes()).hexdigest(), cpu_s=time.process_time() - t0,
             platform=f'Python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__} {platform.platform()}')
    save(tag, r); print('DONE', tag, r['cpu_s'], flush=True)


def step0(G, s, dt):
    """Uncoupled control: same RK4 and kernel with W = 1 (Phi = 0 throughout)."""
    p, v, F, Pi = s; one = np.ones(G.r.shape); zero = np.zeros_like(one)
    a = lambda q, w: e1.accel(G, q, w, one, zero, 0., 2)
    k1v = a(p, v); k2p = v + dt * k1v / 2; k2v = a(G.fix(p + dt * v / 2), k2p); k3p = v + dt * k2v / 2; k3v = a(G.fix(p + dt * k2p / 2), k3p)
    k4p = v + dt * k3v; k4v = a(G.fix(p + dt * k3p), k4p)
    q = G.fix(p + dt * (v + 2 * k2p + 2 * k3p + k4p) / 6)
    return q, e1.ref.tangent(q, v + dt * (k1v + 2 * k2v + 2 * k3v + k4v) / 6) * G.inner, F, Pi


def analyse():
    res = {}
    for f in sorted(OUT.glob('u*_GN*.json')):
        r = json.loads(f.read_text()); t = np.array([x['t'] for x in r['rows']]); out = {}
        for key in ('centroid', 'bcentroid'):
            C = np.array([x[key][0] for x in r['rows']]); e = (t >= .5) & (t <= 2); l = t >= 6
            ve, vl = np.polyfit(t[e], C[e], 1)[0], np.polyfit(t[l], C[l], 1)[0]; vall = np.polyfit(t[t >= .5], C[t >= .5], 1)[0]
            out[key] = dict(v_early=float(ve), v_late=float(vl), rel_change=float(vl / ve - 1), v_all=float(vall))
        H = np.array([x['H_m'] + x['H_phi'] for x in r['rows']]); Pm = np.array([x['P_m'][0] for x in r['rows']]); Pt = Pm + np.array([x['P_phi'] for x in r['rows']])
        out.update(speed_over_P_H=float(out['centroid']['v_all'] / np.mean(Pt / H)), energy_drift=float(np.max(np.abs(H / H[0] - 1))),
                   total_P_drift_over_E=float(np.max(np.abs(Pt - Pt[0])) / H[0]), matter_P_drift_over_E=float(np.max(np.abs(Pm - Pm[0])) / H[0]),
                   field_P_share=float(r['rows'][0]['P_phi'] / Pt[0]), comotion=r['comotion'], degree=[r['degree_initial']['degrees'], r['degree_final']['degrees']], cpu_s=r['cpu_s'])
        res[f.stem] = out
    for u in (0.1, 0.3):
        a, b = res.get(f'u{u:g}_GN0.01'), res.get(f'u{u:g}_GN0')
        if a and b:
            res[f'field_drag_u{u:g}'] = {k: a[k]['rel_change'] - b[k]['rel_change'] for k in ('centroid', 'bcentroid')}
            res[f'speed_ratio_coupled_over_control_u{u:g}'] = a['centroid']['v_all'] / b['centroid']['v_all']
    (OUT / 'summary.json').write_text(json.dumps(res, indent=1) + '\n'); print(json.dumps(res, indent=1))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('phase', choices=['prep', 'run', 'analyse']); ap.add_argument('--u', type=float, default=.3)
    ap.add_argument('--GN', type=float, default=.01); a = ap.parse_args()
    if a.phase == 'prep': prep()
    elif a.phase == 'run': run(a.u, a.GN)
    else: analyse()


if __name__ == '__main__':
    main()
