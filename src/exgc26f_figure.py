"""EX-GC26F figure and separate recomputation from committed records (no quadratic fit).

Gravity displacement g(t) = (D_coupled - D_control)(t) - (same at t = 0) is compared with the double time-integral Q(t) of
the predicted relative acceleration -2 (2E2/E) A / D_control(t)^2. Ratio = least squares of g on Q over t >= 5.
  python src/exgc26f_figure.py -> results/EX-GC26F/recompute.json, paper/figures/exgc26f_attraction.png
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / 'results/EX-GC26F'


def cumint(y, t):
    return np.concatenate([[0], np.cumsum((y[1:] + y[:-1]) / 2 * np.diff(t))])


def main():
    c26 = json.loads((ROOT / 'results/EX-GC26C/B1_GN0.01.json').read_text()); pr = 2 * c26['E2'] / c26['E']
    res = dict(passive_ratio=pr); fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for d, col in (('6', 'C3'), ('9.2', 'C0')):
        a = json.loads((OUT / f'd{d}_GN0.01.json').read_text()); b = json.loads((OUT / f'd{d}_GN0.json').read_text())
        t = np.array([x['t'] for x in a['rows']]); res[f'd{d}'] = {}
        for key, ls in (('bcentroids', '-'), ('centroids', '--')):
            Dc = np.array([x[key][1] - x[key][0] for x in a['rows']]); Dn = np.array([x[key][1] - x[key][0] for x in b['rows']])
            g = Dc - Dn; g -= g[0]; Q = cumint(cumint(-2 * pr * a['A_single'] / Dn ** 2, t), t); sel = t >= 5
            ratio = float(np.sum(g[sel] * Q[sel]) / np.sum(Q[sel] ** 2))
            res[f'd{d}'][key] = dict(ratio_displacement_t5=ratio, ratio_end=float(g[-1] / Q[-1]), g_end=float(g[-1]), Q_end=float(Q[-1]))
            axs[0].plot(t, g, ls=ls, color=col, label=f'd={d}, measured ({"baryon" if key[0] == "b" else "energy"} centroid)')
            if key == 'bcentroids':
                axs[0].plot(t, Q, ':', color=col, lw=2, label=f'd={d}, Newtonian prediction')
                axs[1].plot(t, Dn - Dn[0], ls, color=col, label=f'd={d}, control (no gravity)')
                axs[1].plot(t, Dc - Dc[0], color=col, alpha=.5, label=f'd={d}, coupled')
    axs[0].set(xlabel='time', ylabel='change in separation from gravity', title='Gravitational pull: coupled − control'); axs[0].legend(fontsize=7); axs[0].grid(alpha=.2)
    axs[1].set(xlabel='time', ylabel='separation change (baryon centroid)', title='Raw separations (direct Skyrme push included)'); axs[1].legend(fontsize=7); axs[1].grid(alpha=.2)
    fig.suptitle('EX-GC26F — two Skyrmions attracting through their own Nordström field (exploratory)'); fig.tight_layout()
    fig.savefig(ROOT / 'paper/figures/exgc26f_attraction.png', dpi=150)
    (OUT / 'recompute.json').write_text(json.dumps(res, indent=1) + '\n'); print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
