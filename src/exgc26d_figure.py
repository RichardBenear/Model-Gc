"""EX-GC26D figure from committed records: centroid speed vs time (coupled and control) and field co-motion.

  python src/exgc26d_figure.py -> paper/figures/exgc26d_drag.png
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / 'results/EX-GC26D'; CHK = ROOT / 'review/EX-GC26D_check'


def main():
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for u, col in ((0.1, 'C0'), (0.3, 'C3')):
        for GN, ls in ((0.01, '-'), (0, ':')):
            r = json.loads((OUT / f'u{u:g}_GN{GN:g}.json').read_text()); t = np.array([x['t'] for x in r['rows']])
            C = np.array([x['centroid'][0] for x in r['rows']]); v = np.gradient(C, t)
            axs[0].plot(t[1:-1], v[1:-1] / u, ls=ls, color=col, label=f'u={u}, ' + ('coupled G_N=0.01' if GN else 'control (no field)'))
        f = CHK / f'field_only_u{u:g}.json'
        if f.exists():
            q = json.loads(f.read_text())
            for run, mk in zip(q['runs'], ('x', 'o')):
                d = [x['d'] for x in run['comotion']]
                axs[1].plot(d, [x['fore_over_aft'] for x in run['comotion']], marker=mk, color=col, ls='--' if run['case'] == 'taper' else '-',
                            label=f'u={u} fore/aft, ' + ('PC setup (far taper)' if run['case'] == 'taper' else 'exact far field'))
        rc = json.loads((OUT / f'u{u:g}_GN0.01.json').read_text())['comotion']
        axs[1].plot([x['d'] for x in rc], [x['fore_over_aft'] for x in rc], 's', mfc='none', color=col, ms=9, label=f'u={u} fore/aft, coupled PC run')
    axs[0].set(xlabel='time', ylabel='centroid speed / u', title='Coasting: coupled vs no-field control'); axs[0].legend(fontsize=7); axs[0].grid(alpha=.2)
    axs[1].axhline(1, color='k', lw=.7); axs[1].set(xlabel='distance d from knot', ylabel='Φ(X+d) / Φ(X−d)', title='Field co-motion at t = 8 (1 = no lag)')
    axs[1].legend(fontsize=6.5); axs[1].grid(alpha=.2)
    fig.suptitle('EX-GC26D — boosted Skyrmion with its own Nordström field (exploratory)'); fig.tight_layout()
    fig.savefig(ROOT / 'paper/figures/exgc26d_drag.png', dpi=150)


if __name__ == '__main__':
    main()
