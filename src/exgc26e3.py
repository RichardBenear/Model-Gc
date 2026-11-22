"""EX-GC26E3 (exploratory): does internal (isorotational) kinetic energy fall and resist force like rest energy?

The closed 26E1 B = 1 rest state is set isospinning at angular frequency omega (rotation in the (phi1, phi2) plane,
which for the hedgehog is a spatial rotation about z) and evolved in the same ramped uniform Nordstrom field as 26E1
(Phi = lambda(t) g x, g = 0.0125, h = 0.2, L = 7, T = 5). The spin-up is sudden, so the knot also breathes; both are
internal motion. Fall rate is compared with the static B = 1 knot of 26E1 (same code, field, fits).

Massless Skyrme reference model, not Model Gc; construction ours. Kernel, field, ramp, integrator: 26E1, unchanged.

  python src/exgc26e3.py run --omega 0.5 --mode null|energy|control
  python src/exgc26e3.py analyse
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import argparse, hashlib, json, platform, time
from pathlib import Path
import numpy as np
import scipy
import exgc26e1 as e1

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC26E3'
RAW = ROOT / 'local_raw/EX-GC26E3'
E1OUT = ROOT / 'results/EX-GC26E1'


def save(tag, r):
    OUT.mkdir(parents=True, exist_ok=True); p = OUT / (tag + '.json'); q = p.with_suffix('.tmp')
    q.write_text(json.dumps(r, indent=1, allow_nan=False) + '\n'); q.replace(p)


def spin_velocity(p, omega):
    v = np.zeros_like(p); v[1] = -omega * p[2]; v[2] = omega * p[1]; return v      # tangent: p.v = 0


def extra(G, p, v, t, g):
    """Quadratic-term kinetic energy K2 and gradient energy E2, unweighted and W-weighted (W = e^{2 Phi})."""
    W, _ = e1.weights(G, t, g); fac = G.h ** 3 / e1.NORM; vv = np.sum(v * v, axis=0)
    e2 = e1.fields(G, p, v, W)[3]
    return dict(K2=float(vv.sum() * fac), E2u=float(e2.sum() * fac), K2W=float((W * vv).sum() * fac), E2W=float((W * e2).sum() * fac),
                F_energy_full=float(-2 * g * e1.ramp(t) * np.sum(W * (e2 - vv)) * fac))


def run(omega, mode):
    t0 = time.process_time(); e1.configure(0.2); e1.load(); G = e1.ref.Grid(e1.H, e1.L)
    st = json.loads((E1OUT / 'static_B1.json').read_text()); rp = e1.restpath(1)
    assert hashlib.sha256(rp.read_bytes()).hexdigest() == st['sha256']
    p = np.load(rp); v = spin_velocity(p, omega) * G.inner
    fctl = 2 * e1.GRAV * st['E2'] / st['Q']                          # the 26E1 control constant, unchanged
    g = e1.GRAV if mode == 'energy' else 0.; f = -fctl if mode == 'control' else 0.
    dt = .05 * e1.H; steps = round(e1.T_END / dt); every = max(1, round(.1 / dt)); rows = []; d0 = e1.a2.geometric(p)
    tag = f'w{omega:g}_{mode}'
    for i in range(steps + 1):
        t = i * dt
        if i % every == 0 or i == steps:
            row = e1.measure(G, p, v, t, g, f); row.update(extra(G, p, v, t, g)); rows.append(row)
            if (i // every) % 10 == 0: print(tag, round(t, 2), round(time.process_time() - t0), row['bcentroid'][0], row['H'], flush=True)
        if i < steps:
            p, v = e1.step(G, p, v, t, dt, g, f)
    RAW.mkdir(parents=True, exist_ok=True); end = RAW / (tag + '_end.npz'); np.savez(end, p=p, v=v)
    r = dict(omega=omega, mode=mode, h=e1.H, L=e1.L, T=steps * dt, dt=dt, g=g, f=f, f_control_constant=fctl, E_static=st['E'], E2_static=st['E2'],
             Q_static=st['Q'], rows=rows, degree_initial=d0, degree_final=e1.a2.geometric(p), rest_sha256=st['sha256'],
             end_sha256=hashlib.sha256(end.read_bytes()).hexdigest(), cpu_s=time.process_time() - t0,
             platform=f'Python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__} {platform.platform()}')
    save(tag, r); print('DONE', tag, r['cpu_s'], flush=True)


def qfit(t, y, t0):
    s = t >= t0 - 1e-9; return float(2 * np.polyfit(t[s], y[s], 2)[0])


def fall(r, null, key='bcentroid'):
    t = np.array([x['t'] for x in r['rows']]); y = np.array([x[key][0] for x in r['rows']]); yn = np.array([x[key][0] for x in null['rows']])
    d = (y - y[0]) - (yn - yn[0]); return {f't{t0}': qfit(t, d, t0) / e1.GRAV for t0 in (1, 2)}


def predictions(null, st):
    """From the no-field spinning run (t >= 1 averages): lattice energy-coupling prediction and the alternatives."""
    R = [x for x in null['rows'] if x['t'] >= 1]
    H = np.mean([x['H'] for x in R]); src = np.mean([2 * (x['E2u'] - x['K2']) for x in R]); Q = np.mean([x['Q'] for x in R])
    static_pred = 2 * st['E2'] / st['E']
    fctl = 2 * e1.GRAV * st['E2'] / st['Q']
    return dict(E_total=H, kinetic_fraction=float(1 - st['E'] / H), source_2_E2_minus_K2=float(src),
                a_over_g_energy_lattice=float(-src / H), ratio_energy_lattice=float(src / H / static_pred),
                ratio_no_weight=float(st['E'] / H),                    # internal motion has inertia but no weight
                ratio_no_inertia=float(src / st['E'] / static_pred),   # weight as the action gives, inertia = rest energy only
                a_over_g_control=float(-fctl * Q / (e1.GRAV * H)), ratio_control=float((fctl * Q / (e1.GRAV * H)) / (fctl * st['Q'] / (e1.GRAV * st['E']))))


def analyse():
    st = json.loads((E1OUT / 'static_B1.json').read_text())
    b1e, b1n = (json.loads((E1OUT / f'B1_{m}.json').read_text()) for m in ('energy', 'null'))
    stat = {k: fall(b1e, b1n, k) for k in ('bcentroid', 'centroid')}
    res = dict(static_B1_26E1=stat, static_pred_a_over_g=-2 * st['E2'] / st['E'])
    for f in sorted(OUT.glob('w*_null.json')):
        null = json.loads(f.read_text()); w = null['omega']; out = dict(prediction=predictions(null, st))
        Hn = np.array([x['H'] for x in null['rows']])
        out['null'] = dict(energy_drift=float(np.max(np.abs(Hn / Hn[0] - 1))), max_bcentroid_drift=float(max(abs(x['bcentroid'][0] - null['rows'][0]['bcentroid'][0]) for x in null['rows'])),
                           degree=[null['degree_initial']['degrees'], null['degree_final']['degrees']],
                           K2_range=[min(x['K2'] for x in null['rows']), max(x['K2'] for x in null['rows'])], cpu_s=null['cpu_s'])
        for mode in ('energy', 'control'):
            q = OUT / f'w{w:g}_{mode}.json'
            if not q.exists(): continue
            r = json.loads(q.read_text()); fl = {k: fall(r, null, k) for k in ('bcentroid', 'centroid')}
            out[mode] = dict(a_over_g=fl, ratio_to_static={k: {t: fl[k][t] / stat[k][t] for t in fl[k]} for k in fl},
                             degree=[r['degree_initial']['degrees'], r['degree_final']['degrees']], cpu_s=r['cpu_s'])
        res[f'w{w:g}'] = out
    (OUT / 'summary.json').write_text(json.dumps(res, indent=1) + '\n'); print(json.dumps(res, indent=1))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('phase', choices=['run', 'analyse']); ap.add_argument('--omega', type=float, default=.5)
    ap.add_argument('--mode', choices=['null', 'energy', 'control'], default='null'); a = ap.parse_args()
    run(a.omega, a.mode) if a.phase == 'run' else analyse()


if __name__ == '__main__':
    main()
