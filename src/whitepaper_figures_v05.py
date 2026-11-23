"""White paper v0.5 figure for EX-GC34 (computed from committed records only).

  python src/whitepaper_figures_v05.py   -> paper/figures/wp_gc34_bending_waves.png
Inputs: results/EX-GC34Ar/summary.json, results/EX-GC34D/bend.json, results/EX-GC34D/posthoc_iso.json,
results/EX-GC34E/posthoc_modes.json, results/EX-GC34E/review_unstable.json (panel c, added in revision 2).
"""
import json, os
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'paper/figures'
INK, INK2, GRID = '#1f1f1e', '#5f5e58', '#e4e3dd'
BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'

ar = {r['bg']: r for r in json.loads((ROOT / 'results/EX-GC34Ar/summary.json').read_text())['M2'] if r['mode'] == 'F1'}
bd = json.loads((ROOT / 'results/EX-GC34D/bend.json').read_text())
iso = json.loads((ROOT / 'results/EX-GC34D/posthoc_iso.json').read_text())
modes = json.loads((ROOT / 'results/EX-GC34E/posthoc_modes.json').read_text())['modes']
unst = json.loads((ROOT / 'results/EX-GC34E/review_unstable.json').read_text())['imag_axis']

rows = [
    ('Nordström (prescribed, 34A-r)', ar['B4']['alpha'] / ar['B2']['alpha'], ar['B4']['alpha_ray'] / ar['B2']['alpha_ray']),
    ('lapse only (prescribed, 34A-r)', 1.0, 1.0),
    ('GR-like (prescribed, 34A-r)', ar['B3']['alpha'] / ar['B2']['alpha'], ar['B3']['alpha_ray'] / ar['B2']['alpha_ray']),
    ('Regge, solved (34D, registered)', bd['ratio_field'], bd['ratio_ray']),
    ('Regge, isotropic coords (34D, post hoc)', iso['ratio_field'], iso['ratio_ray']),
]

plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'text.color': INK})
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.2), gridspec_kw=dict(width_ratios=[1.25, 1, 0.8]))

y = np.arange(len(rows))[::-1]; hgt = 0.36
a1.barh(y + hgt / 2 + 0.01, [r[1] for r in rows], height=hgt, color=BLUE, label='wave field (packet)')
a1.barh(y - hgt / 2 - 0.01, [r[2] for r in rows], height=hgt, color=ORANGE, label='ray prediction')
for yy, r in zip(y, rows):
    a1.text(max(r[1], r[2]) + 0.04, yy, f'{abs(r[1]) if abs(r[1]) < 5e-4 else r[1]:.3f} / {abs(r[2]) if abs(r[2]) < 5e-4 else r[2]:.3f}', va='center', fontsize=8, color=INK2)
a1.set_yticks(y); a1.set_yticklabels([r[0] for r in rows])
for xv in (0, 1, 2):
    a1.axvline(xv, color=GRID, lw=1, zorder=0)
a1.set_xlim(-0.05, 2.9); a1.set_xlabel('light deflection ÷ lapse-only deflection  (GR: 2)')
a1.set_title('(a) Light bending on the shared lattice', loc='left', fontsize=10)
a1.legend(loc='upper right', frameon=False, fontsize=8)
for s in ('top', 'right'):
    a1.spines[s].set_visible(False)

cols = {'100': BLUE, '110': ORANGE, '111': AQUA}
mk = {'100': 'o', '110': 's', '111': '^'}
for dn in ('100', '110', '111'):
    xs, ys = [], []
    for m in modes:
        if m['dir'] == dn:
            for sp in [f['speed'] for f in m['fine'] if f['sv1'] < 1e-5]:   # deep zeros (off-axis floor ~1e-6)
                xs.append(m['kh']); ys.append(sp)
    a2.plot(xs, ys, mk[dn], color=cols[dn], ms=7, mec='white', mew=1.2, label=f'gravitational wave, k along ({dn})')
kh = np.array(sorted({m['kh'] for m in modes}))
cl = [next(m['c_light'] for m in modes if m['kh'] == k and m['dir'] == '100') for k in kh]
kk = np.linspace(0.01, 0.21, 100)
a2.plot(kh, cl, '-', color=INK2, lw=2, label='light (34A massless field)')
a2.axhline(1, color=GRID, lw=1, zorder=0)
a2.set_xlabel('kh (wavenumber × lattice spacing)'); a2.set_ylabel('speed ω/k')
a2.set_title('(b) Two oscillating wave modes (34E, post hoc)', loc='left', fontsize=10)
a2.legend(frameon=False, fontsize=8, loc='lower left')
for s in ('top', 'right'):
    a2.spines[s].set_visible(False)
a2.grid(axis='y', color=GRID, lw=0.8)
ucol = {'110': ORANGE, '111': AQUA, 'generic': INK2, 'near100': BLUE}
ulab = {'110': '(110)', '111': '(111)', 'generic': '(0.8, 0.5, 0.33)', 'near100': '(1, 0.05, 0)'}
for dn in ('110', 'near100', 'generic', '111'):
    rr = sorted((u['kh'], u['growth_Im_w_over_k']) for u in unst if u['dir'] == dn)
    a3.plot([r[0] for r in rr], [r[1] for r in rr], '-o', color=ucol[dn], ms=5, lw=1.5, label=f'k along {ulab[dn]}')
a3.set_xscale('log'); a3.set_ylim(0, 2.5)
a3.set_xticks([0.025, 0.05, 0.1, 0.2]); a3.set_xticklabels(['0.025', '0.05', '0.1', '0.2']); a3.minorticks_off()
a3.set_xlabel('kh (log scale)'); a3.set_ylabel('growth rate Im ω / k')
a3.set_title('(c) Growing re-slicing mode (34E, review)', loc='left', fontsize=10)
a3.legend(frameon=False, fontsize=8, loc='lower right')
for s_ in ('top', 'right'):
    a3.spines[s_].set_visible(False)
a3.grid(axis='y', color=GRID, lw=0.8)
fig.tight_layout()
os.makedirs(FIG, exist_ok=True)
fig.savefig(FIG / 'wp_gc34_bending_waves.png', dpi=150)
print('wrote', FIG / 'wp_gc34_bending_waves.png')
