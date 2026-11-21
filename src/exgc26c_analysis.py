"""EX-GC26C analysis and figure from the committed compact records only.

  python src/exgc26c_analysis.py -> results/EX-GC26C/summary.json, paper/figures/exgc26c_field.png
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / 'results/EX-GC26C'
RUNS = [('B1_GN0.01', 1, .01), ('B2_GN0.01', 2, .01), ('B1_GN0.02', 1, .02)]


def main():
    res = {}; fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for tag, B, GN in RUNS:
        p = OUT / f'{tag}.json'
        if not p.exists(): continue
        r = json.loads(p.read_text()); rows = r['rows']; T = r['T']
        prof = r['profiles'][str(int(round(T)))]
        rr = np.array([x['r'] for x in prof]); ph = np.array([x['mean'] for x in prof]); sd = np.array([x['std'] for x in prof])
        sel = (rr >= 3) & (rr <= T - 2)
        A = -float(np.sum(ph[sel] / rr[sel]) / np.sum(1 / rr[sel] ** 2))          # least squares Phi = -A/r
        c2 = np.polyfit(1 / rr[sel], ph[sel], 1); A2, C2 = -c2[0], c2[1]           # with a constant, as a check
        Mg = np.array([x['M_g'] for x in rows]); t = np.array([x['t'] for x in rows])
        Htot = np.array([x['H_m'] + x['H_phi'] for x in rows]); C = np.array([x['centroid'] for x in rows])
        rphi = -ph[sel] * rr[sel] / GN
        res[tag] = dict(B=B, GN=GN, T=T, E=r['E'], two_E2=2 * r['E2'], A_over_GN=A / GN, A_over_GN_with_const=A2 / GN, const=C2,
                        M_g_source_mean=float(Mg.mean()), M_g_source_range=[float(Mg.min()), float(Mg.max())],
                        A_over_GN_div_source=A / GN / float(Mg.mean()), A_over_GN_div_2E2=A / GN / (2 * r['E2']), A_over_GN_div_E=A / GN / r['E'],
                        rphi_spread=float((rphi.max() - rphi.min()) / rphi.mean()), max_angular_rel_std=float(np.max(sd[sel] / np.abs(ph[sel]))),
                        energy_drift=float(np.max(np.abs(Htot / Htot[0] - 1))), H_m_change=float(rows[-1]['H_m'] - rows[0]['H_m']),
                        H_phi_end=rows[-1]['H_phi'], phi_centre_end=rows[-1]['phi_centre'], centroid_drift=float(np.max(np.abs(C - C[0]))),
                        degree_initial=r['degree_initial']['degrees'], degree_final=r['degree_final']['degrees'], cpu_s=r['cpu_s'])
        for tt, prf in r['profiles'].items():
            if tt == '0': continue
            x = np.array([q['r'] for q in prf]); y = -np.array([q['rphi'] for q in prf]) / (GN * Mg.mean())
            if GN == .01: axs[0].plot(x, y, label=f'B={B}, t={tt}', ls='-' if B == 1 else '--')
        axs[1].plot(t, (np.array([x['H_m'] for x in rows]) - rows[0]['H_m']) / r['E'], label=f'B={B} G_N={GN}: matter')
        axs[1].plot(t, np.array([x['H_phi'] for x in rows]) / r['E'], ls=':', label=f'B={B} G_N={GN}: field')
    if 'B1_GN0.01' in res and 'B2_GN0.01' in res:
        a1, a2 = res['B1_GN0.01'], res['B2_GN0.01']
        res['active_mass_ratio_B2_B1'] = dict(measured=a2['A_over_GN'] / a1['A_over_GN'], source_ratio=a2['M_g_source_mean'] / a1['M_g_source_mean'],
                                              lattice_2E2_ratio=a2['two_E2'] / a1['two_E2'], energy_ratio=a2['E'] / a1['E'], no_binding=2.0)
    if 'B1_GN0.02' in res and 'B1_GN0.01' in res:
        res['linearity_A_per_GN_ratio'] = res['B1_GN0.02']['A_over_GN'] / res['B1_GN0.01']['A_over_GN']
    (OUT / 'summary.json').write_text(json.dumps(res, indent=1) + '\n')
    axs[0].axhline(1, color='k', lw=.7, ls=':'); axs[0].set(xlabel='r', ylabel='−rΦ / (G_N M_g)', title='Self-generated field: r·Φ (1 = exact 1/r)', ylim=(0, 1.2)); axs[0].legend(fontsize=7); axs[0].grid(alpha=.2)
    axs[1].set(xlabel='time', ylabel='energy / E', title='Energy exchange (total conserved)'); axs[1].legend(fontsize=7); axs[1].grid(alpha=.2)
    fig.suptitle('EX-GC26C — Skyrmion sourcing its own Nordström field (exploratory)'); fig.tight_layout(); fig.savefig(ROOT / 'paper/figures/exgc26c_field.png', dpi=150)
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
