"""EX-GC26E1 analysis and figure, from the committed compact records only (separate path from the run driver).

  python src/exgc26e1_analysis.py   -> results/EX-GC26E1/summary.json, paper/figures/exgc26e1_fall.png
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC26E1'


def load(B, mode):
    p = OUT / f'B{B}_{mode}.json'
    return json.loads(p.read_text()) if p.exists() else None


def series(r, key):
    return np.array([row[key][0] for row in r['rows']])


def qfit(t, y, t0):
    s = t >= t0 - 1e-9; c = np.polyfit(t[s], y[s], 2)
    return 2 * c[0], float(np.max(np.abs(np.polyval(c, t[s]) - y[s])))


def main():
    st = {B: json.loads((OUT / f'static_B{B}.json').read_text()) for B in (1, 2)}
    g = None; res = dict(static={}, runs={})
    for B in (1, 2):
        s = st[B]; res['static'][B] = dict(E=s['E'], E2=s['E2'], E4=s['E4'], derrick_ratio=s['ratio'], Q=s['Q'],
                                           degree=s['geometric']['degrees'], residual=s['projected_gradient_max'],
                                           wall_energy_fraction=s['wall_energy_fraction'])
    Ebind = 2 * st[1]['E'] - st[2]['E']
    res['binding'] = dict(E_bind=Ebind, fraction=Ebind / (2 * st[1]['E']))
    fctl = 2 * 0.0125 * st[1]['E2'] / st[1]['Q']
    pred = {}
    for B in (1, 2):
        pred[('energy', B)] = -2 * st[B]['E2'] / st[B]['E']
        pred[('control', B)] = -fctl * st[B]['Q'] / (0.0125 * st[B]['E'])
    res['predicted_a_over_g'] = {f'{m}_B{B}': v for (m, B), v in pred.items()}
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for B in (1, 2):
        null = load(B, 'null')
        for mode in ('energy', 'control'):
            r = load(B, mode)
            if r is None:
                continue
            g = 0.0125
            t = np.array([row['t'] for row in r['rows']])
            out = {}
            for cen in ('bcentroid', 'centroid'):
                y = series(r, cen) - series(r, cen)[0]
                if null is not None:
                    y = y - (series(null, cen) - series(null, cen)[0])
                for t0 in (1, 2):
                    a, resid = qfit(t, y, t0); out[f'{cen}_a_over_g_t{t0}'] = a / g; out[f'{cen}_fit_resid_t{t0}'] = resid
                if cen == 'bcentroid':
                    axs[0].plot(t, y, label=f'B={B} {mode}', ls='-' if mode == 'energy' else '--')
            P = np.array([row['P'][0] for row in r['rows']]); F = np.array([row['F_energy' if mode == 'energy' else 'F_control'] for row in r['rows']])
            imp = np.r_[0., np.cumsum((F[1:] + F[:-1]) / 2 * np.diff(t))]
            Pn = P - P[0] - ((np.array([row['P'][0] for row in null['rows']]) - null['rows'][0]['P'][0]) if null else 0)
            Htot = np.array([row['H'] + row['V'] for row in r['rows']]); post = t >= 1
            Ptr = max(np.hypot(row['P'][1], row['P'][2]) for row in r['rows']) / st[B]['E']
            out.update(predicted_a_over_g=pred[(mode, B)], momentum_over_impulse=float(Pn[-1] / imp[-1]),
                       impulse_over_E=float(imp[-1] / st[B]['E']),
                       energy_drift_after_ramp=float(np.max(np.abs(Htot[post] / Htot[post][0] - 1))),
                       transverse_P_over_E=float(Ptr), degree_initial=r['degree_initial']['degrees'], degree_final=r['degree_final']['degrees'],
                       baryon_number_density_end=r['rows'][-1]['B'], cpu_s=r['cpu_s'])
            out['ratio_to_prediction'] = out['bcentroid_a_over_g_t1'] / pred[(mode, B)]
            res['runs'][f'{mode}_B{B}'] = out
        if null is not None:
            t = np.array([row['t'] for row in null['rows']]); y = series(null, 'bcentroid')
            res['runs'][f'null_B{B}'] = dict(max_bcentroid_drift=float(np.max(np.abs(y - y[0]))), max_centroid_drift=float(np.max(np.abs(series(null, 'centroid') - series(null, 'centroid')[0]))),
                                             P_drift_over_E=float(max(abs(row['P'][0] - null['rows'][0]['P'][0]) for row in null['rows']) / st[B]['E']),
                                             energy_drift=float(max(abs(row['H'] / null['rows'][0]['H'] - 1) for row in null['rows'])), cpu_s=null['cpu_s'],
                                             degree_final=null['degree_final']['degrees'])
    R = {}
    for mode in ('energy', 'control'):
        k1, k2 = f'{mode}_B1', f'{mode}_B2'
        if k1 in res['runs'] and k2 in res['runs']:
            for cen in ('bcentroid', 'centroid'):
                for t0 in (1, 2):
                    R[f'{mode}_{cen}_t{t0}'] = res['runs'][k2][f'{cen}_a_over_g_t{t0}'] / res['runs'][k1][f'{cen}_a_over_g_t{t0}']
            R[f'{mode}_predicted'] = pred[(mode, 2)] / pred[(mode, 1)]
    R['continuum_energy'] = 1.0
    res['ratio_B2_over_B1'] = R
    res['cpu_s_total'] = float(sum(v.get('cpu_s', 0) for v in res['runs'].values()) + sum(st[B]['cpu_total'] for B in (1, 2)))
    (OUT / 'summary.json').write_text(json.dumps(res, indent=1) + '\n')
    # figure
    axs[0].set(xlabel='time', ylabel='baryon centroid shift (null subtracted)', title='Fall along x (g = 0.0125)'); axs[0].legend(fontsize=8); axs[0].grid(alpha=.2)
    labels, meas, prd = [], [], []
    for mode in ('energy', 'control'):
        for B in (1, 2):
            k = f'{mode}_B{B}'
            if k in res['runs']:
                labels.append(f'{mode}\nB={B}'); meas.append(-res['runs'][k]['bcentroid_a_over_g_t1']); prd.append(-pred[(mode, B)])
    x = np.arange(len(labels)); axs[1].bar(x - .18, meas, .36, label='measured (baryon centroid, t≥1)'); axs[1].bar(x + .18, prd, .36, label='lattice prediction')
    axs[1].axhline(1, color='k', lw=.7, ls=':'); axs[1].set_xticks(x); axs[1].set_xticklabels(labels); axs[1].set(ylabel='fall acceleration / g'); axs[1].legend(fontsize=8)
    axs[1].set_ylim(0, max(meas + prd + [1]) * 1.15); axs[1].grid(alpha=.2, axis='y')
    fig.suptitle('EX-GC26E1 — single Skyrmion vs bound pair in a Nordström field (exploratory)'); fig.tight_layout()
    fig.savefig(ROOT / 'paper/figures/exgc26e1_fall.png', dpi=150); plt.close(fig)
    print(json.dumps(dict(binding=res['binding'], predicted=res['predicted_a_over_g'], ratios=R,
                          runs={k: {kk: v[kk] for kk in v if 'a_over_g_t1' in kk or kk in ('momentum_over_impulse', 'ratio_to_prediction', 'energy_drift_after_ramp', 'max_bcentroid_drift', 'P_drift_over_E', 'degree_final')} for k, v in res['runs'].items()}), indent=1))


if __name__ == '__main__':
    main()
