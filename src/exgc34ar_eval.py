"""EX-GC34A-r evaluator: as exgc34a_eval.py but the readout is the periodic (circular) centroid, unwrapped in time,
applied identically to the field (cc series) and to the rays (box=...). Registered in spec/EX-GC34Ar_plan.md.

  python src/exgc34ar_eval.py      -> results/EX-GC34Ar/summary.json
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = '1'
import json, sys
from pathlib import Path
import numpy as np
from multiprocessing import Pool
import exgc34a_rays as R

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC34Ar'
BOX = {'M1': ((-8.0, -8.0, -8.0), (16.0, 16.0, 16.0)), 'M2': ((-24.0, -6.0, -6.0), (48.0, 20.0, 12.0))}
H = 0.2; DT = 0.25 * H; SIG = 2.0
MU = {'F1': 0.0, 'F2': 1.0, 'F4': 0.0, 'F5': 0.0, 'sym': 0.0, 'anti': np.sqrt(2 * 0.6 / (1 - 0.4))}
GAM = [[1.1875, 0.101725, 0.037025], [0.101725, 1.031793, 0.084366], [0.037025, 0.084366, 0.830707]]


def gam_exact():
    sys.path.insert(0, str(ROOT / 'src')); import exgc34a
    return exgc34a.gam_B1()


def bg_for(name, scale=1.0):
    if name == 'B0':
        return R.Bg('flat', scale=scale, dt_lf=DT)
    if name in ('B1', 'B1p'):
        return R.Bg('aniso', beta_div=12.0 if name == 'B1' else 16.0, gam=gam_exact(), scale=scale, dt_lf=DT)
    if name == 'B2':
        return R.Bg('bump', gppn=0.0, dt_lf=DT)
    if name == 'B3':
        return R.Bg('bump', gppn=1.0, dt_lf=DT)
    if name == 'B4':
        return R.Bg('bump', nordstrom=True, dt_lf=DT)


def angle(v):
    return float(np.arctan2(v[1], v[0]))


def unwrap(c, box):
    lo, Lb = np.asarray(box[0]), np.asarray(box[1])
    ph = np.unwrap(2 * np.pi * (np.asarray(c) - lo) / Lb, axis=0)
    return lo + Lb * ph / (2 * np.pi)


def pred(args):
    bgname, scale, x0, k0, mu, T, times, box = args
    c = R.rays(bg_for(bgname, scale), x0, k0, SIG, mu, H, T, DT, nq=4, times=times, box=box)
    return unwrap(c, box)


def main():
    recs = {p.stem: json.loads(p.read_text()) for p in sorted(OUT.glob('M[12]_*.json'))}
    jobs = {}
    for tag, r in recs.items():
        hist = r['hist']; t = np.array(hist['t'])
        if tag.startswith('M1'):
            x0 = [0, 0, 0]; k0 = r['k0']; T = 5.0
        else:
            x0 = [-18.0, 6.0, 0.0]; k0 = [3.0, 0, 0]; T = 36.0
        modes = ['sym', 'anti'] if r['field'] == 'F3' else [r['field']]
        for m in modes:
            jobs[(tag, m, 1.0)] = (r['bg'], 1.0, x0, k0, MU[m], T, t, BOX[tag[:2]])
            if m == 'F5':
                jobs[(tag, m, 1.10)] = (r['bg'], 1.10, x0, k0, MU[m], T, t, BOX[tag[:2]])
    keys = list(jobs)
    with Pool(int(os.environ.get('MODELG_WORKERS', 4))) as P:
        res = P.map(pred, [jobs[k] for k in keys])
    preds = dict(zip(keys, res))
    out = {'M1': [], 'M2': []}
    for tag, r in recs.items():
        hist = r['hist']; t = np.array(hist['t'])
        modes = ['sym', 'anti'] if r['field'] == 'F3' else [r['field']]
        for m in modes:
            key = 'c_' + m if r['field'] == 'F3' else 'c'
            c = unwrap(hist['c' + key], BOX[tag[:2]])
            if tag.startswith('M1'):
                v = np.array(R.fit_v(t, c, 1.0)); vp = np.array(R.fit_v(t, preds[(tag, m, 1.0)], 1.0))
                rel = float(np.linalg.norm(v) / np.linalg.norm(vp) - 1)
                ang = float(np.arccos(np.clip(v @ vp / np.linalg.norm(v) / np.linalg.norm(vp), -1, 1)))
                row = dict(tag=tag, mode=m, bg=r['bg'], dir=r['dir'], v=v.tolist(), v_ray=vp.tolist(), rel=rel, ang=ang)
                if m == 'F5':
                    vo = np.array(R.fit_v(t, preds[(tag, m, 1.10)], 1.0))
                    row['v_ray_own'] = vo.tolist()
                    row['rel_own'] = float(np.linalg.norm(v) / np.linalg.norm(vo) - 1)
                    row['ang_own'] = float(np.arccos(np.clip(v @ vo / np.linalg.norm(v) / np.linalg.norm(vo), -1, 1)))
                    row['pass'] = bool(abs(rel) >= 0.03 and abs(row['rel_own']) <= 0.003 and row['ang_own'] <= 0.003)
                elif m == 'F4':
                    row['pass'] = bool(abs(rel) <= 0.01 and ang <= 0.01)
                else:
                    row['pass'] = bool(abs(rel) <= 0.003 and ang <= 0.003)
                row['E_drift'] = float(hist[key + '_E'][-1] / hist[key + '_E'][0] - 1)
                out['M1'].append(row)
            else:
                def alpha(cc):
                    ve = R.fit_v(t[(t >= 2) & (t <= 8)], cc[(t >= 2) & (t <= 8)], 0)
                    vl = R.fit_v(t[t >= 28], cc[t >= 28], 0)
                    return angle(vl) - angle(ve)
                a = alpha(c); ap = alpha(preds[(tag, m, 1.0)])
                th_late = angle(R.fit_v(t[t >= 28], c[t >= 28], 0))
                out['M2'].append(dict(tag=tag, mode=m, bg=r['bg'], mu=float(MU[m]), alpha=a, alpha_ray=ap,
                                      theta_late=th_late,
                                      E_drift=float(hist[key + '_E'][-1] / hist[key + '_E'][0] - 1)))
    # M2 tolerance: 5% of |alpha_ray(B2)| for the same mass
    refB2 = {}
    for row in out['M2']:
        if row['bg'] == 'B2':
            refB2[round(row['mu'], 6)] = abs(row['alpha_ray'])
    for row in out['M2']:
        tol = 0.05 * refB2.get(round(row['mu'], 6), np.nan)
        row['tol'] = tol; row['pass'] = bool(abs(row['alpha'] - row['alpha_ray']) <= tol)
    out['M1_pass'] = f"{sum(r['pass'] for r in out['M1'])}/{len(out['M1'])}"
    out['M2_pass'] = f"{sum(r['pass'] for r in out['M2'])}/{len(out['M2'])}"
    (OUT / 'summary.json').write_text(json.dumps(out, indent=1) + '\n')
    print(out['M1_pass'], out['M2_pass'])


if __name__ == '__main__':
    main()
