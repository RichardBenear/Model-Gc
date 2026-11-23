"""EX-GC34K (exploratory): diagnose the first-order E4 error of the re-slicing-exact kinetic term; gradient family;
full scoring with constraint-compatible synchronous time domain. Plan: spec/EX-GC34K_plan.md.

  python src/exgc34k.py diag    -> results/EX-GC34K/diag.json                 (A1-A3)
  python src/exgc34k.py family  -> results/EX-GC34K/family.json, derive.json, kernels.npz   (B)
  python src/exgc34k.py spec    -> results/EX-GC34K/candidates.json           (C: G1-G3, G5)
  python src/exgc34k.py e4      -> results/EX-GC34K/e4.json                   (C: E4, order, option a)
  python src/exgc34k.py sync    -> results/EX-GC34K/sync.json                 (C: G4'')
  python src/exgc34k.py check   -> results/EX-GC34K/check.json                (C: G7)
Reuses exgc34e/f/g/i/j unchanged.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import itertools, json, platform, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp
import scipy.linalg as sl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import exgc34e as E
import exgc34f as F
import exgc34g as G
import exgc34i as I

OUT = ROOT / 'results/EX-GC34K'
H = E.H; NV = 11; IN, IB, IE = F.IN, F.IB, F.IE
MU = -1.0
DUMP = G.DUMP
SEED = 3434 + 11
ET = E.ET.astype(float)
DLS = np.linalg.solve(ET.T @ ET, ET.T)                       # least-squares gradient from the 7 forward edges


def g_gamma(gam):
    return lambda q: ((1 + gam) * np.exp(1j * q) - (1 + 2 * gam) + gam * np.exp(-1j * q)) / H


def g_kuhn(q):
    return DLS @ (np.exp(1j * (ET @ q)) - 1) / H


def sigma(gf, q):
    return -I.L0 * ((I.DH * (np.exp(1j * (ET @ q)) - 1)[:, None]) @ gf(q))


def half_grid(N):
    ks = 2 * np.pi * np.arange(N) / N; out = []
    for m in itertools.product(range(N), repeat=3):
        mm = tuple((-np.array(m)) % N)
        if m <= mm:
            out.append(np.array([ks[i] for i in m]))
    return out


def solve_space(gf, R, st):
    B = I.basis(R); n = len(B); N = 2 * R + 6
    rows = []; rhs = []
    for q in half_grid(N):
        s0 = sigma(gf, q)
        Ren = E.static_kernel(st, q[None])[0][IE, IN] / I.L0 / MU
        Mq = np.zeros((7, n), complex)
        for k, (o, i, j, w) in enumerate(B):
            ph = np.exp(1j * (q @ np.array(o, float)))
            if o == (0, 0, 0):
                Mq[i, k] += s0[j]
                if i != j:
                    Mq[j, k] += s0[i]
            else:
                Mq[i, k] += ph * s0[j]
                Mq[j, k] += np.conj(ph) * s0[i]
        rows.append(Mq); rhs.append(Ren)
    A = np.vstack(rows); b = np.concatenate(rhs)
    A = np.vstack([A.real, A.imag]); b = np.concatenate([b.real, b.imag])
    W = np.array([w for (_, _, _, w) in B]); Aw = A / W
    del A
    U, s, Vh = np.linalg.svd(Aw, full_matrices=False)
    r = int((s > s[0] * 1e-11).sum())
    xp = Vh[:r].T @ ((U[:, :r].T @ b) / s[:r])
    cons = float(np.linalg.norm(Aw @ xp - b) / np.linalg.norm(b))
    Nsp = Vh[r:].T
    cT, _ = I.extract_T(R); xT = W * I.vec_from_coef(cT, R)
    x = xp + Nsp @ (Nsp.T @ (xT - xp))
    return dict(R=R, n=n, rank=r, dim=Nsp.shape[1], consistent_rel=cons, exists=bool(cons <= 1e-10),
                rel_change_vs_T=float(np.linalg.norm(x - xT) / np.linalg.norm(xT))), x, Nsp, W, xT


def I0g(coef, gf, st, rng):
    iden = []; herm = []; nul = []
    for _ in range(200):
        q = rng.uniform(-np.pi, np.pi, 3); K = I.K_from(coef, None, q)
        r_ = E.static_kernel(st, q[None])[0][IE, IN] / I.L0
        iden.append(np.linalg.norm(K @ (MU * sigma(gf, q)) - r_) / np.linalg.norm(r_))
        herm.append(np.abs(K - K.conj().T).max() / np.abs(K).max())
    kin = I.make_kin(coef)
    for _ in range(50):
        q = rng.uniform(-np.pi, np.pi, 3); w = rng.uniform(-3, 3) / H
        M = kin(q, w) + E.static_kernel(st, q[None])[0]
        v = np.zeros(NV, complex); v[IN] = -1j * w; v[1:4] = MU * gf(q)
        nul.append(np.linalg.norm(M @ v) / (np.linalg.norm(M, 2) * np.linalg.norm(v)))
    out = dict(identity_rel=float(max(iden)), hermitian=float(max(herm)), null_rel=float(max(nul)))
    out['I0'] = bool(out['identity_rel'] <= 1e-12 and out['hermitian'] <= 1e-14 and out['null_rel'] <= 1e-12)
    return out


# ---------------- E4 per case and component ----------------
CASES = []
for qdir in ((1, 0, 0), (1, 1, 1), (0.8, 0.5, 0.33)):
    qd = G.unit(qdir)
    perp = np.cross(qd, [0, 0, 1.0]) if abs(qd[2]) < 0.9 else np.cross(qd, [1.0, 0, 0])
    for vn, vd in (('par', qd), ('perp', perp), ('obl', qd + 1.3 * perp)):
        CASES.append((qdir, vn, vd))


def e4_detail(kin, st, kh, src_ph=0.0, read_ph=0.0, vmag=1e-3):
    """Per case: O(v) Riemann error vector (lattice - GR) / max|GR|, and static error."""
    out = []
    for qdir, vn, vd in CASES:
        qd = G.unit(qdir); q = kh * qd; k = q / H
        v = vmag * np.array(vd) / np.linalg.norm(vd)
        R = {}
        for sgn in (1, -1):
            w = k @ (sgn * v)
            Mm = kin(q, w) + E.static_kernel(st, q[None])[0]; Gm = E.gauge_modes(q[None], w)[0]
            f = np.zeros(NV, complex); f[0] = E.KAPPA; f[1:4] = -E.KAPPA * sgn * v * np.exp(src_ph * 1j * q)
            x, *_ = np.linalg.lstsq(Mm + Gm @ Gm.conj().T, f, rcond=1e-14)
            x2 = x.copy(); x2[1:4] = x[1:4] * np.exp(read_ph * 1j * q)
            h = E.metric_from_x(x2, q); kl = np.array([-w, *k])
            Rl = E.riemann(h, kl) * H ** 3
            hg, klg, _ = E.gr_solution(k, sgn * v); Rg = E.riemann(hg, klg)
            R[sgn] = (Rl, Rg)
        R1l = (R[1][0] - R[-1][0]) / 2; R1g = (R[1][1] - R[-1][1]) / 2
        R0l = (R[1][0] + R[-1][0]) / 2; R0g = (R[1][1] + R[-1][1]) / 2
        d = (R1l - R1g) / np.abs(R1g).max()
        out.append(dict(qdir=list(qdir), v=vn, comp_err=d, err_Ov=float(np.abs(d).max()),
                        err_static=float(np.abs(R0l - R0g).max() / np.abs(R0g).max())))
    return out


def worst(rows):
    i = int(np.argmax([r['err_Ov'] for r in rows])); c = int(np.argmax(np.abs(rows[i]['comp_err'])))
    return i, c


def load_J1():
    z = np.load(ROOT / 'results/EX-GC34J/kernels.npz')
    return {tuple(int(v) for v in k.split('|')[1].split(',')): z[k] for k in z.files if k.startswith('J1|')}


def diag():
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); rng = np.random.default_rng(SEED); res = {}
    J1 = load_J1(); kin = I.make_kin(J1)
    # A1
    tab = {kh: e4_detail(kin, st, kh) for kh in (0.0125, 0.025, 0.05, 0.1)}
    i, c = worst(tab[0.05])
    comps = [(a, b, cc, d) for a in range(4) for b in range(a + 1, 4) for cc in range(4) for d in range(cc + 1, 4)]
    res['A1'] = dict(worst_case=dict(qdir=CASES[i][0], v=CASES[i][1]), worst_component=comps[c],
                     err_by_kh={str(kh): [r['err_Ov'] for r in tab[kh]] for kh in tab},
                     worst_comp_by_kh={str(kh): float(abs(tab[kh][i]['comp_err'][c])) for kh in tab},
                     per_case_c1=[float(np.abs(tab[0.0125][j]['comp_err']).max() / 0.0125) for j in range(9)],
                     case_labels=[[list(cs[0]), cs[1]] for cs in CASES])
    c1J1 = float(abs(tab[0.0125][i]['comp_err'][c]) / 0.0125)
    print('A1', res['A1']['worst_case'], comps[c], 'c1', c1J1, flush=True)

    # A2
    def c1_of(coef):
        rows = [e4_detail(I.make_kin(coef), st, kh) for kh in (0.0125, 0.025)]
        cw = abs(rows[0][i]['comp_err'][c]) / 0.0125
        cmax = max(r['err_Ov'] for r in rows[0]) / 0.0125
        return float(cw), float(cmax), float(abs(rows[1][i]['comp_err'][c]) / abs(rows[0][i]['comp_err'][c]))
    samples = []
    info1, x1, N1, W1, xT1 = solve_space(g_gamma(0.0), 1, st)
    base = x1; dist = np.linalg.norm(x1 - xT1)
    for scale in (0.1, 0.3):
        for _ in range(10):
            z = N1 @ rng.standard_normal(N1.shape[1]); z *= scale * dist / np.linalg.norm(z)
            coef = I.coef_from_vec((base + z) / W1, 1)
            cw, cm, rat = c1_of(coef); samples.append(dict(R=1, scale=scale, c1_worst=cw, c1_max=cm, ratio=rat))
    info2, x2, N2, W2, xT2 = solve_space(g_gamma(0.0), 2, st)
    dist2 = np.linalg.norm(x2 - xT2)
    for _ in range(10):
        z = N2 @ rng.standard_normal(N2.shape[1]); z *= 0.1 * dist2 / np.linalg.norm(z)
        coef = I.coef_from_vec((x2 + z) / W2, 2)
        cw, cm, rat = c1_of(coef); samples.append(dict(R=2, scale=0.1, c1_worst=cw, c1_max=cm, ratio=rat))
    zJ = np.load(ROOT / 'results/EX-GC34J/kernels.npz')
    for nm in ('J1', 'J2', 'J3'):
        coef = {tuple(int(v) for v in k.split('|')[1].split(',')): zJ[k] for k in zJ.files if k.startswith(nm + '|')}
        cw, cm, rat = c1_of(coef); samples.append(dict(R=nm, scale=0.0, c1_worst=cw, c1_max=cm, ratio=rat))
    cws = np.array([s['c1_worst'] for s in samples]); cms = np.array([s['c1_max'] for s in samples])
    spread = float((cws.max() - cws.min()) / cws.mean())
    res['A2'] = dict(samples=samples, spread_worst=spread, spread_max=float((cms.max() - cms.min()) / cms.mean()),
                     verdict='A-fixed' if spread < 0.01 else 'A-variable', c1_min=float(cws.min()), c1_J1=c1J1,
                     space_R1=info1, space_R2=info2)
    print('A2', res['A2']['verdict'], spread, res['A2']['spread_max'], flush=True)
    if spread >= 0.01:
        jmin = int(np.argmin(cws)); res['A2']['min_sample'] = samples[jmin]
    # A3
    a3 = []
    for sp in (-0.5, 0.0, 0.5):
        for rp in (-0.5, 0.0, 0.5):
            r = e4_detail(kin, st, 0.0125, sp, rp); r2 = e4_detail(kin, st, 0.025, sp, rp)
            a3.append(dict(src=sp, read=rp, c1_max=max(x['err_Ov'] for x in r) / 0.0125,
                           ratio=max(x['err_Ov'] for x in r2) / max(x['err_Ov'] for x in r)))
    res['A3'] = a3
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'diag.json').write_text(json.dumps(res, **DUMP))
    print('A3', [(a['src'], a['read'], round(a['c1_max'], 4), round(a['ratio'], 3)) for a in a3]); print('cpu', res['cpu_s'])



def diag_g1():
    """A2 follow-up required by the plan (omitted from the first diag code): G1 on the reduced direction set
    for the smallest-c1 kernel. Regenerates that kernel by replaying the same random stream as diag()."""
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); rng = np.random.default_rng(SEED)
    d = json.loads((OUT / 'diag.json').read_text()); samp = d['A2']['samples']
    jmin = int(np.argmin([x['c1_worst'] for x in samp]))
    assert samp[jmin]['R'] == 1, 'min sample expected at range 1'
    info1, x1, N1, W1, xT1 = solve_space(g_gamma(0.0), 1, st); dist = np.linalg.norm(x1 - xT1)
    k = 0
    for scale in (0.1, 0.3):
        for _ in range(10):
            z = N1 @ rng.standard_normal(N1.shape[1]); z *= scale * dist / np.linalg.norm(z)
            if k == jmin:
                coef = I.coef_from_vec((x1 + z) / W1, 1)
            k += 1
    kin = I.make_kin(coef)
    # check replay: c1 of the regenerated kernel
    r = e4_detail(kin, st, 0.0125)
    c1_replay = max(x['err_Ov'] for x in r) / 0.0125
    reg = json.loads((ROOT / 'results/EX-GC34F/spectrum.json').read_text())['directions']
    names = ['100', '010', '001', 'd01+', 'd01-', 'd02+', 'd02-', 'd12+', 'd12-', 'b1++', 'b1+-', 'b1-+', 'b1--', 'generic']
    KV = G.KernelG(st, kin); rg = np.random.default_rng(SEED + 7); grow = []; nunc = 0
    for nm in names:
        for kh in (0.05, 0.2, 1.0, 3.0):
            q = kh * np.array(reg[nm]); a = F.analyse_posthoc(KV.pencil(q), q, rg)
            nunc += sum(not rr['confirmed'] for rr in a['roots'])
            grow += [(nm, kh, rr['w_over_k']) for rr in a['roots'] if rr['confirmed'] and rr['cls'] == 'growing']
    d['A2']['min_kernel_G1'] = dict(sample_index=jmin, c1_max_replay=float(c1_replay), c1_max_recorded=samp[jmin]['c1_max'],
                                    dirs=names, kh=[0.05, 0.2, 1.0, 3.0], n_growing=len(grow), growing=grow[:20],
                                    n_unconfirmed=nunc, G1=len(grow) == 0, cpu_s=time.process_time() - t0,
                                    note='plan A2 follow-up, run after diag (code omission, not a change of rule)')
    (OUT / 'diag.json').write_text(json.dumps(d, **DUMP))
    print('min-kernel G1', len(grow) == 0, len(grow), 'unconf', nunc, 'c1 replay', c1_replay, samp[jmin]['c1_max'], 'cpu', time.process_time() - t0)

# ---------------- Part B ----------------
FAMILY = [('g0', g_gamma(0.0)), ('gm025', g_gamma(-0.25)), ('gm05', g_gamma(-0.5)), ('gm075', g_gamma(-0.75)),
          ('gm1', g_gamma(-1.0)), ('gkuhn', g_kuhn)]


def family():
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); rng = np.random.default_rng(SEED + 1)
    res = dict(members={}); save = {}; sols = {}
    J1 = load_J1()
    sols['J1'] = dict(exists=True, I0=True, note='34J J1 (forward gradient, range 1, mu = -1, closest to T)')
    for o, Ko in J1.items():
        save['J1|%d,%d,%d' % o] = Ko
    for name, gf in FAMILY:
        for R in (1, 2, 3):
            if R == 3 and name not in ('gm05', 'gm025'):
                continue
            if R > 1 and res['members'].get('%s_R%d' % (name, R - 1), {}).get('exists'):
                continue
            t1 = time.process_time()
            info, x, Nsp, W, xT = solve_space(gf, R, st)
            info['cpu_s'] = time.process_time() - t1
            key = '%s_R%d' % (name, R)
            if info['exists']:
                coef = I.coef_from_vec(x / W, R)
                info.update(I0g(coef, gf, st, rng))
                kin = I.make_kin(coef)
                e = {kh: max(r['err_Ov'] for r in e4_detail(kin, st, kh)) for kh in (0.025, 0.05)}
                info['E4'] = {str(k): v for k, v in e.items()}; info['E4_ratio'] = e[0.05] / e[0.025]
                info['qualifies_for_C'] = bool(info['I0'] and info['E4_ratio'] >= 3)
                if info['qualifies_for_C']:
                    sols[key] = dict(exists=True, I0=True)
                    for o, Ko in coef.items():
                        save['%s|%d,%d,%d' % ((key,) + o)] = Ko
            res['members'][key] = info
            print(key, {k: v for k, v in info.items()}, flush=True)
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'family.json').write_text(json.dumps(res, **DUMP))
    (OUT / 'derive.json').write_text(json.dumps(dict(solutions=sols), **DUMP))
    np.savez_compressed(OUT / 'kernels.npz', **save)
    print('cpu', res['cpu_s'])


# ---------------- Part C ----------------
def use_here():
    I.OUT = OUT
    return I.load_cands()


def spec():
    cands = use_here()
    G.CANDS = {n: I.make_kin(c) for n, (c, s) in cands.items()}
    G.OUT = OUT
    G.spec()


def check():
    use_here(); I.check()


def e4():
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); cands = use_here(); res = {}
    for name, (coef, s) in cands.items():
        kin = I.make_kin(coef)
        rows = {kh: max(max(r['err_Ov'], r['err_static']) for r in e4_detail(kin, st, kh)) for kh in (0.0125, 0.025, 0.03, 0.04, 0.05, 0.1)}
        ok = [kh for kh, v in rows.items() if v <= 1e-2]
        res[name] = dict(rows={str(k): v for k, v in rows.items()}, ratio=rows[0.05] / rows[0.025],
                         E4_pass=bool(rows[0.05] <= 1e-2), order_pass=bool(rows[0.05] / rows[0.025] >= 3),
                         largest_kh_le_1e2=max(ok) if ok else None)
        print(name, res[name], flush=True)
    res['cpu_s'] = time.process_time() - t0
    (OUT / 'e4.json').write_text(json.dumps(res, **DUMP))


def sync():
    """G4'': synchronous gauge, constraint-compatible seeds from the outset."""
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); cands = use_here(); res = {}
    for name, (coef, s) in cands.items():
        KV = G.KernelG(st, I.make_kin(coef)); runs = []
        for dn, d in (('generic', (0.8, 0.5, 0.33)), ('d01+', (1, 1, 0))):
            q = 0.05 * G.unit(d); k = 0.05 / H; P = KV.pencil(q)
            Aee = P['A'][np.ix_(IE, IE)]; Cee = P['C'][np.ix_(IE, IE)]; Fm = np.linalg.solve(Aee, Cee)
            L = np.block([[np.zeros((7, 7)), np.eye(7)], [Fm, np.zeros((7, 7))]])
            CnE = P['C'][IN, IE]; CbE = P['C'][np.ix_(IB, IE)]; BbE = P['B'][np.ix_(IB, IE)]; AbE = P['A'][np.ix_(IB, IE)]
            C0 = np.vstack([np.r_[CnE, np.zeros(7)][None], np.hstack([CbE - AbE @ Fm, 1j * BbE])])
            Cmap = np.vstack([C0, C0 @ L]); nC = np.linalg.norm(Cmap, 2)
            mu_, V = np.linalg.eig(L); zc = np.abs(mu_) <= 0.01 * k
            waves = [j for j in range(14) if not zc[j] and abs(abs(mu_[j]) / k - 1) < 0.01 and mu_[j].imag < 0]
            jw = min(waves, key=lambda j: abs(mu_[j])); wv = abs(mu_[jw]); zw = V[:, jw] / np.linalg.norm(V[:, jw])
            nonz = [V[:, j] / np.linalg.norm(V[:, j]) for j in range(14) if not zc[j] and j != jw]
            nz_ok = [v for v in nonz if np.linalg.norm(Cmap @ v) <= 1e-8 * nC]
            _, Zs, sdim = sl.schur(L.astype(complex), output='complex', sort=lambda z: abs(z) <= 0.01 * k)
            Z = Zs[:, :sdim]; _, sz, Vz = np.linalg.svd(Cmap @ Z)
            rz = int((sz > sz[0] * 1e-10).sum()) if sz.size and sz[0] > 0 else 0
            Zk = Z @ Vz[rz:].conj().T
            sb = sum(nz_ok) if nz_ok else np.zeros(14); sc_ = Zk.sum(1) if Zk.shape[1] else np.zeros(14)
            wmax = np.abs(mu_).max(); T = 40 * 2 * np.pi / wv

            def tid(z):
                x = np.zeros(NV, complex); x[IE] = Fm @ z[:7]; return 0.5 * E.metric_from_x(x, q)[1:, 1:]

            def crel(z):
                return float(np.linalg.norm(C0 @ z) / (np.linalg.norm(C0, 2) * np.linalg.norm(z)))
            out = dict(dir=dn, wave_w_over_k=float(wv / k), n_nonzero_seeds=len(nz_ok), n_nonzero_total=len(nonz),
                       zero_cluster_dim=int(sdim), zero_kept=int(Zk.shape[1]))
            for tag, fac in (('dt', 0.02), ('dt2', 0.01)):
                ns = int(np.ceil(T / (fac / wmax))); dt = T / ns
                I14 = np.eye(14); Ms = np.linalg.solve(I14 - 0.5 * dt * L, I14 + 0.5 * dt * L)
                every = max(1, ns // 4000); Mev = np.linalg.matrix_power(Ms, every)
                za = zw.copy(); zb = zw + 1e-8 * sb; zz = zw + 1e-8 * sc_
                Ea0 = np.linalg.norm(tid(za)); Eam = 0.0; amp = 0.0; db = 0.0; dc = 0.0; cm = 0.0
                for i in range(0, ns + 1, every):
                    Ea = tid(za); nE = np.linalg.norm(Ea); Eam = max(Eam, nE); amp = max(amp, abs(nE / Ea0 - 1))
                    db = max(db, np.linalg.norm(tid(zb) - Ea)); dc = max(dc, np.linalg.norm(tid(zz) - Ea))
                    cm = max(cm, crel(za), crel(zb), crel(zz))
                    za = Mev @ za; zb = Mev @ zb; zz = Mev @ zz
                out[tag] = dict(dt=dt, steps=ns, sample_every=every, max_constraint=cm, wave_tidal_amp_change=float(amp),
                                seeds_nonzero_rel=float(db / Eam), seeds_zero_rel=float(dc / Eam))
            a, b = out['dt'], out['dt2']
            out['G4pp'] = bool(all(t['max_constraint'] <= 1e-9 and t['wave_tidal_amp_change'] <= 1e-6 and t['seeds_nonzero_rel'] <= 1e-6
                                   and t['seeds_zero_rel'] <= 1e-6 for t in (a, b)))
            runs.append(out); print(name, out, flush=True)
        res[name] = dict(runs=runs, G4pp=all(r['G4pp'] for r in runs))
    res['cpu_s'] = time.process_time() - t0
    (OUT / 'sync.json').write_text(json.dumps(res, **DUMP))


def posthoc_zb():
    """POST HOC (labelled): pointwise obstruction for the centred gradient. Wherever the gradient vanishes,
    sigma_xi = 0, so K sigma_xi = C_en/(mu l0) needs C_en = 0 there for ANY bounded K (any range)."""
    st = E.tet_static(); rows = []
    for name, gam in (('gm025', -0.25), ('gm05', -0.5), ('gm075', -0.75)):
        for q in itertools.product((0.0, np.pi), repeat=3):
            q = np.array(q)
            if not q.any():
                continue
            ren = E.static_kernel(st, q[None])[0][IE, IN] / I.L0
            rows.append(dict(member=name, q_over_pi=list(q / np.pi), max_abs_g=float(np.abs(g_gamma(gam)(q)).max()),
                             max_abs_sigma=float(np.abs(sigma(g_gamma(gam), q)).max()), max_abs_rhs=float(np.abs(ren).max())))
    (OUT / 'posthoc_zb.json').write_text(json.dumps(dict(note='POST HOC', rows=rows), **DUMP))
    for r in rows:
        print(r)


if __name__ == '__main__':
    dict(diag=diag, diag_g1=diag_g1, family=family, spec=spec, e4=e4, sync=sync, check=check, posthoc_zb=posthoc_zb)[sys.argv[1]]()
