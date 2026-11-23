"""EX-GC34G (exploratory): candidate lattice kinetic terms for the shift-link geometry. Plan: spec/EX-GC34G_plan.md.

  python src/exgc34g.py spec    -> results/EX-GC34G/candidates.json   (G1, G2, G3, G5 + reported)
  python src/exgc34g.py evolve  -> results/EX-GC34G/evolve.json       (G4)
  python src/exgc34g.py check   -> results/EX-GC34G/check.json        (G7, real-space recomputation)
Reuses src/exgc34e.py and src/exgc34f.py unchanged. Convention e^{-i w t}.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import json, platform, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp
import scipy.linalg as sl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import exgc34d as D
import exgc34e as E
import exgc34f as F

OUT = ROOT / 'results/EX-GC34G'
H = E.H; ET = E.ET; NV = 11; G1 = F.ghat(1.0)
IN, IB, IE, RED = F.IN, F.IB, F.IE, F.RED
SEED = 3434 + 7
DUMP = dict(indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))

# ---------------- edge groups ----------------
TET_GROUPS = [[(tuple(o[i]), t) for (i, j, t) in D.tet_edges(o)] for o in D.TETS]
CUBE_EDGES = F.CUBE_EDGES
VERT_EDGES = [((0, 0, 0), t) for t in range(7)] + [(tuple(-ET[t]), t) for t in range(7)]


def fit_rows(edges):
    """W (ne x 6) and, for given q, w: S (ne x 11) with Fourier phases relative to the cell base vertex."""
    Wm = np.zeros((len(edges), 6))
    for e, (oi, t) in enumerate(edges):
        Dv = ET[t] * H
        for c, (a, b) in enumerate(E.SYM):
            Wm[e, c] = Dv[a] * Dv[b] * (1 if a == b else 2)
    return Wm


def strain_S(edges, q, w):
    S = np.zeros((len(edges), NV), complex)
    for e, (oi, t) in enumerate(edges):
        oi = np.array(oi, float); oj = oi + ET[t]
        Dv = ET[t] * H; l0 = np.linalg.norm(Dv); Dh = Dv / l0
        pi = np.exp(1j * (q @ oi)); pj = np.exp(1j * (q @ oj))
        S[e, 4 + t] += l0 * (-1j * w) * pi
        for a in range(3):
            S[e, 1 + a] += -l0 * Dh[a] * (pj - pi)
    return S


PINV = {'cube': np.linalg.pinv(fit_rows(CUBE_EDGES)), 'vert': np.linalg.pinv(fit_rows(VERT_EDGES)),
        'tet': [np.linalg.inv(fit_rows(g)) for g in TET_GROUPS]}


def K_tets(q, w):
    return [Pi @ strain_S(g, q, w) for Pi, g in zip(PINV['tet'], TET_GROUPS)]


def kin_T(q, w):
    return sum(H ** 3 / 6 * A.conj().T @ G1 @ A for A in K_tets(q, w))


def kin_A(q, w):
    A = PINV['cube'] @ strain_S(CUBE_EDGES, q, w); return H ** 3 * A.conj().T @ G1 @ A


def kin_B(q, w):
    A = sum(K_tets(q, w)) / 6; return H ** 3 * A.conj().T @ G1 @ A


def kin_C(q, w):
    A = PINV['vert'] @ strain_S(VERT_EDGES, q, w); return H ** 3 * A.conj().T @ G1 @ A


def blend(s):
    return lambda q, w: s * kin_T(q, w) + (1 - s) * kin_A(q, w)


CANDS = {'T': kin_T, 'A': kin_A, 'B': kin_B, 'C': kin_C, 'D0.1': blend(0.1), 'D0.25': blend(0.25), 'D0.5': blend(0.5)}


class KernelG(F.Kernel):
    def __init__(self, st, kin):
        super().__init__(st); self.kfun = kin

    def T(self, q, w):
        return self.kfun(q, w)


def unit(d):
    return F.unit(d)


# ---------------- spectrum (G1, G2, G3, G5) ----------------
def spec():
    t0 = time.process_time()
    mp.mp.dps = 50
    st = E.tet_static(); rng = np.random.default_rng(SEED)
    reg34f = json.loads((ROOT / 'results/EX-GC34F/spectrum.json').read_text())
    dirs = dict(reg34f['directions']); appr = dict(reg34f['approach'])
    cases_def = [(n, kh) for kh in (0.0125, 0.025, 0.05, 0.1, 0.2) for n in dirs] + [(n, 0.05) for n in appr]
    for n, d in (('100', (1, 0, 0)), ('d01+', (1, 1, 0)), ('b1++', (1, 1, 1))):
        for kh in (0.5, 1.0, 2.0, 3.0):
            cases_def.append((n, kh))
        cases_def.append((n, float(np.pi / np.abs(unit(d)).max())))
    alld = dict(dirs); alld.update(appr)
    K34 = D.kernel(D.build_terms()[0], np.array([[0.3, -0.7, 1.1]]))[0]
    sel = [0] + list(range(4, 11))
    res = dict(convention='e^{-i w t}; zero-shift gauge; roots confirmed by argument principle (34F analyse_posthoc)', candidates={})
    for cn, kf in CANDS.items():
        KV = KernelG(st, kf)
        out = {}
        qz = np.array([0.3, -0.7, 1.1])
        out['G5_static_vs_34D'] = float(np.abs(KV.M(qz, 0.0)[np.ix_(sel, sel)] - K34).max() / np.abs(K34).max())
        out['G5_gauge_null'] = float(max(np.linalg.norm(KV.M(qz, w) @ E.gauge_modes(qz[None], w)[0]) / np.linalg.norm(KV.M(qz, w))
                                         for w in (0.0, 0.7, -1.3)))
        out['hermitian'] = float(max(np.abs(KV.M(qz, w) - KV.M(qz, w).conj().T).max() for w in (0.7, -1.3)))
        cases = []
        for name, kh in cases_def:
            q = kh * np.array(alld[name]) if name in alld else kh * unit({'100': (1, 0, 0), 'd01+': (1, 1, 0), 'b1++': (1, 1, 1)}[name])
            P = KV.pencil(q)
            a = F.analyse_posthoc(P, q, rng)
            sz = F.static_zero_modes(P, q)
            n11 = F.normal_rank_deficit(P, kh / H, rows=list(range(NV)), cols=list(range(NV)))
            roots = [r for r in a['roots'] if r['confirmed']]
            cases.append(dict(dir=name, kh=kh, nk=a['nk'], nullity_full11=n11, n_exact_reslicing=n11 - 3,
                              static=sz, light_speed=float(E.light_speed(q)), n_unconfirmed=sum(not r['confirmed'] for r in a['roots']),
                              zero_winding=a['zero_winding'],
                              roots=[dict(cls=r['cls'], w_over_k=r['w_over_k'], wh=r['wh'], resid11=r['resid11'],
                                          tt=r['tt_frac_of_nongauge'], ng=r['non_gauge_frac']) for r in roots]))
        # G1
        grow = [(c['dir'], c['kh'], r['w_over_k']) for c in cases for r in c['roots'] if r['cls'] == 'growing']
        out['G1_n_growing'] = len(grow); out['G1_growing_examples'] = grow[:20]
        out['G1_max_growth_over_k'] = max([g[2][1] for g in grow] or [0.0])
        out['G1'] = len(grow) == 0
        # G2
        g2 = []
        for c in cases:
            if c['kh'] <= 0.05 + 1e-12:
                low = [r for r in c['roots'] if r['cls'] == 'propagating' and 0 < r['wh'][0] < 0.5]
                other_low = [r for r in c['roots'] if r['cls'] not in ('propagating',) and abs(complex(*r['wh'])) < 0.5
                             and not (r['cls'] == 'propagating_neg')]
                dc = max([abs(r['w_over_k'][0] - c['light_speed']) for r in low] or [np.inf])
                g2.append(dict(dir=c['dir'], kh=c['kh'], n_low=len(low), max_dc=dc, n_other_low=len(other_low),
                               low_speeds=sorted(r['w_over_k'][0] for r in low)))
        out['G2_cases'] = g2
        out['G2_n_bad'] = sum(not (x['n_low'] == 2 and x['max_dc'] <= 1e-3 and x['n_other_low'] == 0) for x in g2)
        out['G2'] = out['G2_n_bad'] == 0
        # G3
        g3bad = [(c['dir'], c['kh'], c['static'].get('n_chains', 0), c['n_exact_reslicing']) for c in cases
                 if c['static'].get('n_chains', 0) != c['n_exact_reslicing']]
        out['G3_mismatches'] = g3bad; out['G3'] = len(g3bad) == 0
        # reported: E2 breaking, third modes
        brk = []
        for name, d in (('100', (1, 0, 0)), ('b1++', (1, 1, 1)), ('generic', (0.8, 0.5, 0.33))):
            for kh in (0.05, 0.1):
                q = kh * unit(d); k = kh / H; P = KV.pencil(q); vals = []
                for w in (0.0, 0.5 * k):
                    Mm = F.Mw(P, w); Gm = E.gauge_modes(q[None], w)[0]
                    Qc = sl.null_space(Gm.conj().T); sv = np.linalg.svd(Qc.conj().T @ Mm @ Qc, compute_uv=False)
                    vals.append(float(sv[-1] / sv[0]))
                brk.append(dict(dir=name, kh=kh, smin=vals))
        out['reslicing_breaking'] = brk
        # G5 E4
        orig = E.full_kernel
        E.full_kernel = lambda st_, qq, w, KV=KV: KV.M(qq[0], w)[None]
        try:
            e4 = []
            for qdir in ((1, 0, 0), (1, 1, 1), (0.8, 0.5, 0.33)):
                qd = unit(qdir)
                perp = np.cross(qd, [0, 0, 1.0]) if abs(qd[2]) < 0.9 else np.cross(qd, [1.0, 0, 0])
                for vname, vdir in (('par', qd), ('perp', perp), ('obl', qd + 1.3 * perp)):
                    r = E.e4_case(st, qdir, 0.05, vdir)
                    e4.append(dict(qdir=list(qdir), v=vname, err_Ov=r['err_Ov'], err_static=r['err_static'], mom=r['momentum_residual']))
        finally:
            E.full_kernel = orig
        out['G5_E4'] = e4; out['G5_E4_max'] = max(max(r['err_Ov'], r['err_static']) for r in e4)
        out['G5'] = bool(out['G5_static_vs_34D'] <= 1e-12 and out['G5_gauge_null'] <= 1e-10 and out['G5_E4_max'] <= 1e-2)
        # wave speed table
        out['wave_table'] = [dict(dir=c['dir'], kh=c['kh'], speeds=sorted(r['w_over_k'][0] for r in c['roots'] if r['cls'] == 'propagating' and r['wh'][0] < 0.5),
                                  c_light=c['light_speed'])
                             for c in cases if c['dir'] in ('100', 'd01+', 'b1++') and c['kh'] in (0.05, 0.1, 0.2)]
        out['n_unconfirmed_total'] = sum(c['n_unconfirmed'] for c in cases)
        out['cases'] = cases
        res['candidates'][cn] = out
        print(cn, 'G1', out['G1'], out['G1_n_growing'], 'max %.3g' % out['G1_max_growth_over_k'], '| G2', out['G2'], out['G2_n_bad'],
              '| G3', out['G3'], len(g3bad), '| G5', out['G5'], 'static %.1e gauge %.1e E4 %.2e' % (out['G5_static_vs_34D'], out['G5_gauge_null'], out['G5_E4_max']),
              '| unconf', out['n_unconfirmed_total'], '| cpu %.0f' % (time.process_time() - t0), flush=True)
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'candidates.json').write_text(json.dumps(res, **DUMP))
    print('cpu', res['cpu_s'])


# ---------------- G4 time evolution ----------------
def evolve(override=None, outname='evolve.json', label='registered G4', min_seed=0.0):
    t0 = time.process_time()
    mp.mp.dps = 50
    st = E.tet_static(); rng = np.random.default_rng(SEED + 1)
    cj = json.loads((OUT / 'candidates.json').read_text())['candidates']
    todo = override or [c for c in CANDS if cj[c]['G1'] and cj[c]['G2'] and cj[c]['G3'] and cj[c]['G5']]
    res = dict(label=label, scored=todo, runs={})
    for cn in todo:
        KV = KernelG(st, CANDS[cn]); rr = []
        for name, d in (('generic', (0.8, 0.5, 0.33)), ('d01+', (1, 1, 0))):
            kh = 0.05; q = kh * unit(d); k = kh / H; P = KV.pencil(q)
            a = F.analyse_posthoc(P, q, rng, keep_vectors=True)
            modes = [r for r in a['roots'] if r['confirmed']]
            wave = min([r for r in modes if r['cls'] == 'propagating'], key=lambda r: abs(r['w_over_k'][0] - 1))
            others = [r for r in modes if r is not wave and abs(complex(*r['w_over_k'])) > min_seed]
            n8 = len(RED)
            Cr = P['C'][np.ix_(RED, RED)]; Ar = P['A'][np.ix_(RED, RED)]; Br = P['B'][np.ix_(RED, RED)]
            Em = np.block([[np.eye(n8), np.zeros((n8, n8))], [np.zeros((n8, n8)), Ar]])
            Fm = np.block([[np.zeros((n8, n8)), np.eye(n8)], [Cr, 1j * Br]])
            zof = lambda r: np.r_[r['_x'][RED], -1j * r['_w'] * r['_x'][RED]]
            zw = zof(wave); Z = np.array([zof(r) for r in others]).T
            om = max(abs(r['_w']) for r in modes)
            Tend = 40 * 2 * np.pi / wave['_w'].real
            out = dict(dir=name, kh=kh, wave_w_over_k=wave['w_over_k'], n_seeded=len(others),
                       seeded_w_over_k=[r['w_over_k'] for r in others])
            for tag, fac in (('dt', 0.02), ('dt2', 0.01)):
                dt0 = fac / om; ns = int(np.ceil(Tend / dt0)); dt = Tend / ns
                Ms = np.linalg.solve(Em - 0.5 * dt * Fm, Em + 0.5 * dt * Fm)
                z2 = zw.copy(); z3 = zw + 1e-8 * np.linalg.norm(zw) * (Z / np.linalg.norm(Z, axis=0)).sum(1)
                c0 = np.linalg.lstsq(Z, z3 - z2, rcond=None)[0]
                lap, mom, amp, grow = [], [], [], []
                nC0 = np.linalg.norm(P['C'][IN, RED]); nCb = np.linalg.norm(P['C'][np.ix_(IB, RED)], 2); nBb = np.linalg.norm(P['B'][np.ix_(IB, RED)], 2)
                for i in range(ns + 1):
                    if i % 10 == 0 or i == ns:
                        for z in (z2, z3):
                            x = z[:n8]; v = z[n8:]
                            lap.append(abs(P['C'][IN, RED] @ x) / (nC0 * np.linalg.norm(x)))
                            mr = P['C'][np.ix_(IB, RED)] @ x + 1j * P['B'][np.ix_(IB, RED)] @ v
                            mom.append(np.linalg.norm(mr) / (nCb * np.linalg.norm(x) + nBb * np.linalg.norm(v)))
                        amp.append(abs(np.linalg.norm(z2) / np.linalg.norm(zw) - 1))
                        c = np.linalg.lstsq(Z, z3 - z2, rcond=None)[0]
                        grow.append(float(np.max(np.abs(c) / np.abs(c0))))
                    z2 = Ms @ z2; z3 = Ms @ z3
                out[tag] = dict(dt=dt, steps=ns, max_lapse=float(max(lap)), max_mom=float(max(mom)), wave_amp_change=float(max(amp)),
                                max_seed_growth_factor=float(max(grow)))
            out['G4'] = all(out[t]['max_lapse'] <= 1e-10 and out[t]['max_mom'] <= 1e-10 and out[t]['wave_amp_change'] <= 1e-8
                            and out[t]['max_seed_growth_factor'] <= 1 + 1e-6 for t in ('dt', 'dt2'))
            rr.append(out)
            print(cn, name, 'G4', out['G4'], {t: {kk: ('%.3e' % v if isinstance(v, float) else v) for kk, v in out[t].items()} for t in ('dt', 'dt2')}, flush=True)
        res['runs'][cn] = rr
    res['cpu_s'] = time.process_time() - t0
    (OUT / outname).write_text(json.dumps(res, **DUMP))
    print('cpu', res['cpu_s'])


# ---------------- G7 real-space recomputation ----------------
def check():
    N = 8; rng = np.random.default_rng(SEED + 2)
    beta = rng.standard_normal((3, N, N, N)); ldot = rng.standard_normal((7, N, N, N)) * H
    sh = lambda f, off: np.roll(f, shift=tuple(-np.array(off)), axis=(-3, -2, -1))
    def rates(edges):
        s = []
        for (oi, t) in edges:
            oj = tuple(np.array(oi) + ET[t]); Dv = ET[t] * H; l0 = np.linalg.norm(Dv); Dh = Dv / l0
            db = sh(beta, oj) - sh(beta, oi)
            s.append(l0 * (sh(ldot[t], oi) - np.einsum('a,a...->...', Dh, db)))
        return np.stack(s)
    def energy(Kf, V):
        return (V / 2 * np.einsum('c...,cd,d...->...', Kf, G1, Kf)).sum()
    Kt = [np.einsum('cd,d...->c...', np.linalg.pinv(fit_rows(g)), rates(g)) for g in TET_GROUPS]
    direct = {'T': sum(energy(Kk, H ** 3 / 6) for Kk in Kt),
              'A': energy(np.einsum('cd,d...->c...', np.linalg.pinv(fit_rows(CUBE_EDGES)), rates(CUBE_EDGES)), H ** 3),
              'B': energy(sum(Kt) / 6, H ** 3),
              'C': energy(np.einsum('cd,d...->c...', np.linalg.pinv(fit_rows(VERT_EDGES)), rates(VERT_EDGES)), H ** 3)}
    for s in (0.1, 0.25, 0.5):
        direct['D%g' % s] = s * direct['T'] + (1 - s) * direct['A']
    x = np.zeros((11, N, N, N), complex); x[1:4] = beta; x[4:] = ldot / (-1j)
    xh = np.fft.fftn(x, axes=(1, 2, 3)); kk = 2 * np.pi * np.fft.fftfreq(N)
    out = {}
    for cn, kf in CANDS.items():
        four = 0.0
        for i in range(N):
            for j in range(N):
                for l in range(N):
                    q = np.array([kk[i], kk[j], kk[l]]); v = xh[:, i, j, l]
                    four += 0.5 * (v.conj() @ kf(q, 1.0) @ v).real
        four /= N ** 3
        out[cn] = dict(direct=float(direct[cn]), fourier=float(four), rel=float(abs(direct[cn] / four - 1)))
        print(cn, out[cn], flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'check.json').write_text(json.dumps(out, **DUMP))


if __name__ == '__main__':
    {'spec': spec, 'evolve': evolve, 'check': check,
     'evolve_posthoc_B': lambda: evolve(['B'], 'evolve_posthoc_B.json', 'POST HOC: candidate B failed G1-G3 as registered; time-domain look only'),
     'evolve_posthoc_B2': lambda: evolve(['B'], 'evolve_posthoc_B_noclusters.json', 'POST HOC: B, seeds exclude split zero-cluster roots |w| < 0.1 k', 0.1)}[sys.argv[1]]()
