"""EX-GC26E3 separate recomputation (momentum ledger, inertia from momentum vs centroid velocity) and figure.

  python src/exgc26e3_figure.py -> results/EX-GC26E3/recompute.json, paper/figures/exgc26e3_spin.png
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; E1 = ROOT / 'results/EX-GC26E1'; E3 = ROOT / 'results/EX-GC26E3'; G = 0.0125


def ld(p): return json.loads(p.read_text())


def main():
    cases = [('static (26E1)', E1 / 'B1_energy.json', E1 / 'B1_null.json', 'F_energy'),
             ('ω = 0.5', E3 / 'w0.5_energy.json', E3 / 'w0.5_null.json', 'F_energy_full'),
             ('ω = 0.7', E3 / 'w0.7_energy.json', E3 / 'w0.7_null.json', 'F_energy_full'),
             ('ω = 0.7 control', E3 / 'w0.7_control.json', E3 / 'w0.7_null.json', 'F_control')]
    res = {}; fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for (name, fe, fn, fk), col in zip(cases, ('k', 'C0', 'C3', 'C2')):
        r, n = ld(fe), ld(fn); t = np.array([x['t'] for x in r['rows']]); s = t >= 1
        H = float(np.mean([x['H'] for x in n['rows'] if x['t'] >= 1]))
        P = np.array([x['P'][0] for x in r['rows']]) - np.array([x['P'][0] for x in n['rows']]); F = np.array([x[fk] for x in r['rows']])
        I = np.concatenate([[0], np.cumsum((F[1:] + F[:-1]) / 2 * np.diff(t))])
        out = dict(E_total=H, force_over_gE=float(F[s].mean() / (G * H)), dPdt_over_force=float(np.polyfit(t[s], P[s], 1)[0] / F[s].mean()),
                   momentum_over_impulse_end=float(P[-1] / I[-1]))
        for key, ls in (('bcentroid', '-'), ('centroid', '--')):
            y = np.array([x[key][0] for x in r['rows']]) - np.array([x[key][0] for x in n['rows']]); y -= y[0]
            c = np.polyfit(t[s], y[s], 2); vel = np.polyval(np.polyder(c), t[s])
            out[key] = dict(a_over_g=float(2 * c[0] / G), inertia_dPdv_over_E=float(np.polyfit(vel, P[s], 1)[0] / H))
            if 'control' not in name: axs[0].plot(t, y, ls=ls, color=col, label=f'{name}, {"baryon" if key[0] == "b" else "energy"} centroid')
        if 'control' not in name:
            Bc = np.array([x['bcentroid'][0] - x['centroid'][0] for x in r['rows']]) - np.array([x['bcentroid'][0] - x['centroid'][0] for x in n['rows']])
            axs[1].plot(t, Bc - Bc[0], color=col, label=name); out['baryon_minus_energy_centroid_end'] = float(Bc[-1] - Bc[0])
        res[name] = out
    axs[0].set(xlabel='time', ylabel='fall displacement (null subtracted)', title='Fall of static and spinning knots'); axs[0].legend(fontsize=6.5); axs[0].grid(alpha=.2)
    axs[1].set(xlabel='time', ylabel='baryon − energy centroid (null subtracted)', title='The two centroids separate when the knot spins'); axs[1].legend(fontsize=7); axs[1].grid(alpha=.2)
    fig.suptitle('EX-GC26E3 — isospinning Skyrmion in a uniform Nordström field (exploratory)'); fig.tight_layout()
    fig.savefig(ROOT / 'paper/figures/exgc26e3_spin.png', dpi=150)
    (E3 / 'recompute.json').write_text(json.dumps(res, indent=1) + '\n'); print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
