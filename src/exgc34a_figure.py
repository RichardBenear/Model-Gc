"""EX-GC34A-r figure from committed records: results/EX-GC34Ar/summary.json -> paper/figures/exgc34ar_cone.png"""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
s = json.loads((ROOT / 'results/EX-GC34Ar/summary.json').read_text())
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
lab = {'F1': 'massless', 'sym': 'two-comp (sym)', 'F4': 'sigma model', 'F2': 'massive μ=1', 'anti': 'two-comp (anti) μ=√2'}
cols = {'F1': '#1f77b4', 'sym': '#17becf', 'F4': '#9467bd', 'F2': '#d62728', 'anti': '#ff7f0e'}
order = ['F1', 'sym', 'F4', 'F2', 'anti']
bgs = ['B4', 'B2', 'B3']; names = {'B4': 'Nordström\n(conformal)', 'B2': 'lapse only\n(PPN γ=0)', 'B3': 'GR-like\n(γ=1)'}
for j, m in enumerate(order):
    for i, b in enumerate(bgs):
        r = [x for x in s['M2'] if x['bg'] == b and x['mode'] == m][0]
        xp = i + (j - 2) * 0.14
        ax[0].bar(xp, -r['alpha'], 0.13, color=cols[m], label=lab[m] if i == 0 else None)
        ax[0].plot([xp - 0.065, xp + 0.065], [-r['alpha_ray']] * 2, 'k-', lw=1.5)
ax[0].set_xticks(range(3)); ax[0].set_xticklabels([names[b] for b in bgs])
ax[0].set_ylabel('deflection toward mass (rad)'); ax[0].set_title('Bending past a prescribed weak field\n(bars: fields; black ticks: ray prediction)')
ax[0].legend(fontsize=8)
for j, m in enumerate(['F1', 'F2', 'sym', 'anti', 'F4', 'F5']):
    rows = [x for x in s['M1'] if x['mode'] == m]
    y = [100 * x['rel'] for x in rows]
    ax[1].scatter(np.full(len(y), j) + np.linspace(-0.25, 0.25, len(y)), y, s=14, color=cols.get(m, 'k'))
ax[1].axhspan(-0.3, 0.3, color='g', alpha=0.15, label='±0.3% pass band')
ax[1].set_xticks(range(6)); ax[1].set_xticklabels(['massless', 'massive', 'two-comp\nsym', 'two-comp\nanti', 'sigma', 'control:\nown links'])
ax[1].set_ylabel('speed vs shared-structure prediction (%)'); ax[1].set_title('One cone on flat and anisotropic backgrounds\n(3 backgrounds × 3 directions)')
ax[1].legend(fontsize=8)
fig.suptitle('EX-GC34A-r (EXPLORATORY, Model Gc construction): every field uses one shared link structure', fontsize=10)
fig.tight_layout(); out = ROOT / 'paper/figures/exgc34ar_cone.png'; fig.savefig(out, dpi=130); print(out)
