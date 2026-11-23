"""EX-GC34J (exploratory): re-slicing-exact kinetic term, mu = -1, first-order cure, synchronous-gauge time domain.
Plan: spec/EX-GC34J_plan.md.

  python src/exgc34j.py derive -> results/EX-GC34J/derive.json + kernels.npz  (J1-J3, I0)
  python src/exgc34j.py spec   -> results/EX-GC34J/candidates.json            (34G G1-G3, G5)
  python src/exgc34j.py e4     -> results/EX-GC34J/e4.json                    (E4 at four spacings, order)
  python src/exgc34j.py sync   -> results/EX-GC34J/sync.json                  (G4')
  python src/exgc34j.py check  -> results/EX-GC34J/check.json                 (G7)
Reuses exgc34e/f/g/i unchanged.
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

OUT = ROOT / 'results/EX-GC34J'
H = E.H; NV = 11; IN, IB, IE = F.IN, F.IB, F.IE
MU = -1.0
DUMP = G.DUMP
SEED = 3434 + 10


def system(R, st):
    """Identity A c = b (real, unweighted coefficients c in I.basis(R)); K sigma0 = R_en/(mu l0)."""
    B = I.basis(R); n = len(B)
    N = 2 * R + 6; ks = 2 * np.pi * np.arange(N) / N
    rows = []; rhs = []
    for a, b, c in itertools.product(range(N), repeat=3):
        q = np.array([ks[a], ks[b], ks[c]])
        s0 = I.sigma0('f', q)
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
    return np.vstack([A.real, A.imag]), np.concatenate([b.real, b.imag]), B, N


def first_order_map(R):
    """C1 c = vec(sum_o o_k K_o), k = 0..2 (3 x 49 rows)."""
    B = I.basis(R); n = len(B)
    C1 = np.zeros((3 * 49, n))
    for kk in range(n):
        e = np.zeros(n); e[kk] = 1.0
        coef = I.coef_from_vec(e, R)
        m = sum(np.array(o, float)[:, None, None] * Ko[None] for o, Ko in coef.items())
        C1[:, kk] = m.reshape(-1)
    return C1


def derive_cand(R, st, match_first):
    A, b, B, N = system(R, st)
    W = np.array([w for (_, _, _, w) in B])
    Aw = A / W
    U, s, Vh = np.linalg.svd(Aw, full_matrices=False)
    r = int((s > s[0] * 1e-11).sum())
    xp = Vh[:r].T @ ((U[:, :r].T @ b) / s[:r])                    # particular (min-norm) solution
    consistent = float(np.linalg.norm(Aw @ xp - b) / np.linalg.norm(b))
    Nsp = Vh[r:].T                                                 # orthonormal null space (weighted variables)
    cT, _ = I.extract_T(R); xT = W * I.vec_from_coef(cT, R)
    out = dict(R=R, match_first=match_first, n_unknowns=len(B), rank=r, solution_space_dim=Nsp.shape[1],
               consistent_rel=consistent, exists=bool(consistent <= 1e-10), mu=MU, grid=N)
    C1w = first_order_map(R) / W; d1 = C1w @ xT
    if match_first:
        M1 = C1w @ Nsp
        z1, *_ = np.linalg.lstsq(M1, d1 - C1w @ xp, rcond=1e-12)
        U1, s1, V1h = np.linalg.svd(M1, full_matrices=True)
        r1 = int((s1 > s1[0] * 1e-11).sum())
        N2 = Nsp @ V1h[r1:].T
        x0 = xp + Nsp @ z1
        x = x0 + N2 @ (N2.T @ (xT - x0))
        out['first_order_space_dim'] = int(N2.shape[1])
    else:
        x = xp + Nsp @ (Nsp.T @ (xT - xp))
    out['first_order_mismatch_rel'] = float(np.linalg.norm(C1w @ x - d1) / np.linalg.norm(d1))
    out['rel_change_vs_T'] = float(np.linalg.norm(x - xT) / np.linalg.norm(xT))
    coef = I.coef_from_vec(x / W, R)
    return out, coef


def derive():
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); rng = np.random.default_rng(SEED)
    res = dict(solutions={}); save = {}
    for name, R, mf in (('J1', 1, False), ('J2', 1, True), ('J3', 2, True)):
        out, coef = derive_cand(R, st, mf)
        out.update(I.I0(coef, R, 'f', MU, st, rng))
        res['solutions'][name] = out
        for o, Ko in coef.items():
            save['%s|%d,%d,%d' % ((name,) + o)] = Ko
        print(name, out, flush=True)
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'derive.json').write_text(json.dumps(res, **DUMP))
    np.savez_compressed(OUT / 'kernels.npz', **save)
    print('cpu', res['cpu_s'])


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
    st = E.tet_static(); cands = use_here()
    orig = E.full_kernel; res = {}
    for name, (coef, s) in cands.items():
        kin = I.make_kin(coef)
        E.full_kernel = lambda st_, qq, w, kin=kin: (kin(qq[0], w) + E.static_kernel(st, qq)[0])[None]
        rows = []
        for kh in (0.0125, 0.025, 0.05, 0.1):
            e = []
            for qdir in ((1, 0, 0), (1, 1, 1), (0.8, 0.5, 0.33)):
                qd = G.unit(qdir)
                perp = np.cross(qd, [0, 0, 1.0]) if abs(qd[2]) < 0.9 else np.cross(qd, [1.0, 0, 0])
                for vd in (qd, perp, qd + 1.3 * perp):
                    r = E.e4_case(st, qdir, kh, vd); e.append(max(r['err_Ov'], r['err_static']))
            rows.append(dict(kh=kh, max_err=float(max(e))))
        E.full_kernel = orig
        m = {r['kh']: r['max_err'] for r in rows}
        res[name] = dict(rows=rows, ratio_005_0025=m[0.05] / m[0.025], E4_le_1e2=bool(m[0.05] <= 1e-2),
                         order_pass=bool(m[0.05] / m[0.025] >= 3))
        print(name, res[name], flush=True)
    res['cpu_s'] = time.process_time() - t0
    (OUT / 'e4.json').write_text(json.dumps(res, **DUMP))


def sync():
    """G4': synchronous gauge (n = beta = 0): A_ee e'' = C_ee e. Runs (a) wave, (b) + non-zero seeds, (c) + zero-cluster seeds."""
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); cands = use_here(); res = {}
    for name, (coef, s) in cands.items():
        KV = G.KernelG(st, I.make_kin(coef)); runs = []
        for dn, d in (('generic', (0.8, 0.5, 0.33)), ('d01+', (1, 1, 0))):
            q = 0.05 * G.unit(d); k = 0.05 / H; P = KV.pencil(q)
            Aee = P['A'][np.ix_(IE, IE)]; Cee = P['C'][np.ix_(IE, IE)]; Fm = np.linalg.solve(Aee, Cee)
            L = np.block([[np.zeros((7, 7)), np.eye(7)], [Fm, np.zeros((7, 7))]])
            mu_, V = np.linalg.eig(L)
            zc = np.abs(mu_) <= 0.01 * k
            waves = [j for j in range(14) if not zc[j] and abs(abs(mu_[j]) / k - 1) < 0.01 and mu_[j].imag < 0]
            jw = min(waves, key=lambda j: abs(mu_[j])); wv = abs(mu_[jw])
            zw = V[:, jw] / np.linalg.norm(V[:, jw])
            seeds_b = sum(V[:, j] / np.linalg.norm(V[:, j]) for j in range(14) if not zc[j] and j != jw)
            _, Zs, sdim = sl.schur(L.astype(complex), output='complex', sort=lambda z: abs(z) <= 0.01 * k)
            seeds_c = Zs[:, :sdim].sum(1)
            wmax = np.abs(mu_).max(); T = 40 * 2 * np.pi / wv
            CnE = P['C'][IN, IE]; CbE = P['C'][np.ix_(IB, IE)]; BbE = P['B'][np.ix_(IB, IE)]; AbE = P['A'][np.ix_(IB, IE)]

            def tid(z):
                x = np.zeros(NV, complex); x[IE] = Fm @ z[:7]
                return 0.5 * E.metric_from_x(x, q)[1:, 1:]

            def cons(z):
                e, ed = z[:7], z[7:]
                lap = abs(CnE @ e) / (np.linalg.norm(CnE) * np.linalg.norm(e))
                mr = CbE @ e + 1j * BbE @ ed - AbE @ (Fm @ e)
                mom = np.linalg.norm(mr) / (np.linalg.norm(CbE, 2) * np.linalg.norm(e) + np.linalg.norm(BbE, 2) * np.linalg.norm(ed)
                                            + np.linalg.norm(AbE, 2) * np.linalg.norm(Fm @ e) + 1e-300)
                return lap, mom
            out = dict(dir=dn, wave_w_over_k=float(wv / k), n_zero_cluster=int(sdim), wmax_over_k=float(wmax / k))
            for tag, fac in (('dt', 0.02), ('dt2', 0.01)):
                ns = int(np.ceil(T / (fac / wmax))); dt = T / ns
                I14 = np.eye(14); Ms = np.linalg.solve(I14 - 0.5 * dt * L, I14 + 0.5 * dt * L)
                every = max(1, ns // 4000); Mev = np.linalg.matrix_power(Ms, every)
                za = zw.copy(); zb = zw + 1e-8 * seeds_b; zcc = zw + 1e-8 * seeds_c
                Ea0 = np.linalg.norm(tid(za)); amp = 0.0; db = 0.0; dc = 0.0; Eamax = 0.0; lapm = 0.0; momm = 0.0
                for i in range(0, ns + 1, every):
                    Ea = tid(za); nEa = np.linalg.norm(Ea); Eamax = max(Eamax, nEa)
                    amp = max(amp, abs(nEa / Ea0 - 1))
                    db = max(db, np.linalg.norm(tid(zb) - Ea)); dc = max(dc, np.linalg.norm(tid(zcc) - Ea))
                    for z in (za, zb, zcc):
                        l_, m_ = cons(z); lapm = max(lapm, l_); momm = max(momm, m_)
                    za = Mev @ za; zb = Mev @ zb; zcc = Mev @ zcc
                out[tag] = dict(dt=dt, steps=ns, sample_every=every, max_lapse=float(lapm), max_mom=float(momm),
                                wave_tidal_amp_change=float(amp), seeds_nonzero_rel=float(db / Eamax), seeds_zero_rel=float(dc / Eamax))
            out['G4p'] = all(out[t]['max_lapse'] <= 1e-9 and out[t]['max_mom'] <= 1e-9 and out[t]['wave_tidal_amp_change'] <= 1e-6
                             and out[t]['seeds_nonzero_rel'] <= 1e-6 and out[t]['seeds_zero_rel'] <= 1e-6 for t in ('dt', 'dt2'))
            runs.append(out); print(name, out, flush=True)
        res[name] = dict(runs=runs, G4p=all(r['G4p'] for r in runs))
    res['cpu_s'] = time.process_time() - t0
    (OUT / 'sync.json').write_text(json.dumps(res, **DUMP))


if __name__ == '__main__' and sys.argv[1] != 'sync_posthoc':
    dict(derive=derive, spec=spec, e4=e4, sync=sync, check=check)[sys.argv[1]]()


def sync_posthoc():
    """POST HOC: G4' with seeds restricted to the constraint-satisfying subspace, constraint residuals per run."""
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); cands = use_here(); res = dict(label='POST HOC: constraint-filtered seeds; per-run constraint residuals')
    for name, (coef, s) in cands.items():
        KV = G.KernelG(st, I.make_kin(coef)); runs = []
        for dn, d in (('generic', (0.8, 0.5, 0.33)), ('d01+', (1, 1, 0))):
            q = 0.05 * G.unit(d); k = 0.05 / H; P = KV.pencil(q)
            Aee = P['A'][np.ix_(IE, IE)]; Cee = P['C'][np.ix_(IE, IE)]; Fm = np.linalg.solve(Aee, Cee)
            L = np.block([[np.zeros((7, 7)), np.eye(7)], [Fm, np.zeros((7, 7))]])
            CnE = P['C'][IN, IE]; CbE = P['C'][np.ix_(IB, IE)]; BbE = P['B'][np.ix_(IB, IE)]; AbE = P['A'][np.ix_(IB, IE)]
            Cmap = np.vstack([np.r_[CnE, np.zeros(7)][None], np.hstack([CbE - AbE @ Fm, 1j * BbE])])   # 4 x 14 on z = (e, e')
            Cmap = np.vstack([Cmap, Cmap @ L])                                                       # and their time derivatives
            Uc, sc, Vc = np.linalg.svd(Cmap); rc = int((sc > sc[0] * 1e-10).sum()); Pc = Vc[rc:].conj().T   # constraint subspace
            mu_, V = np.linalg.eig(L); zc = np.abs(mu_) <= 0.01 * k
            waves = [j for j in range(14) if not zc[j] and abs(abs(mu_[j]) / k - 1) < 0.01 and mu_[j].imag < 0]
            jw = min(waves, key=lambda j: abs(mu_[j])); wv = abs(mu_[jw]); zw = V[:, jw] / np.linalg.norm(V[:, jw])
            _, Zs, sdim = sl.schur(L.astype(complex), output='complex', sort=lambda z: abs(z) <= 0.01 * k)
            Z = Zs[:, :sdim]
            # zero-cluster seeds inside the constraint subspace
            Mz = Cmap @ Z; Uz, sz, Vz = np.linalg.svd(Mz); rz = int((sz > max(sz[0], 1e-300) * 1e-10).sum()) if sz.size else 0
            Zk = Z @ Vz[rz:].conj().T
            nonz = [V[:, j] / np.linalg.norm(V[:, j]) for j in range(14) if not zc[j] and j != jw]
            nz_ok = [v for v in nonz if np.linalg.norm(Cmap @ v) <= 1e-8 * np.linalg.norm(Cmap, 2)]
            seeds_b = sum(nz_ok) if nz_ok else np.zeros(14); seeds_c = Zk.sum(1) if Zk.shape[1] else np.zeros(14)
            wmax = np.abs(mu_).max(); T = 40 * 2 * np.pi / wv
            def tid(z):
                x = np.zeros(NV, complex); x[IE] = Fm @ z[:7]; return 0.5 * E.metric_from_x(x, q)[1:, 1:]
            def crel(z):
                return float(np.linalg.norm(Cmap[:4] @ z) / (np.linalg.norm(Cmap[:4], 2) * np.linalg.norm(z)))
            out = dict(dir=dn, wave_constraint_rel=crel(zw), n_nonzero_seeds=len(nz_ok), n_nonzero_total=len(nonz),
                       nonzero_violating_wk=[float(abs(mu_[j]) / k) for j in range(14) if not zc[j] and j != jw and np.linalg.norm(Cmap @ (V[:, j] / np.linalg.norm(V[:, j]))) > 1e-8 * np.linalg.norm(Cmap, 2)],
                       zero_cluster_dim=int(sdim), zero_kept=int(Zk.shape[1]))
            for tag, fac in (('dt', 0.02), ('dt2', 0.01)):
                ns = int(np.ceil(T / (fac / wmax))); dt = T / ns
                I14 = np.eye(14); Ms = np.linalg.solve(I14 - 0.5 * dt * L, I14 + 0.5 * dt * L)
                every = max(1, ns // 4000); Mev = np.linalg.matrix_power(Ms, every)
                za = zw.copy(); zb = zw + 1e-8 * seeds_b; zz = zw + 1e-8 * seeds_c
                Eam = 0.0; db = 0.0; dc = 0.0; ca = cb = cc = 0.0
                for i in range(0, ns + 1, every):
                    Ea = tid(za); Eam = max(Eam, np.linalg.norm(Ea))
                    db = max(db, np.linalg.norm(tid(zb) - Ea)); dc = max(dc, np.linalg.norm(tid(zz) - Ea))
                    ca = max(ca, crel(za)); cb = max(cb, crel(zb)); cc = max(cc, crel(zz))
                    za = Mev @ za; zb = Mev @ zb; zz = Mev @ zz
                out[tag] = dict(seeds_nonzero_rel=float(db / Eam), seeds_zero_rel=float(dc / Eam), cons_wave=ca, cons_b=cb, cons_c=cc)
            runs.append(out); print(name, out, flush=True)
        res[name] = runs
    res['cpu_s'] = time.process_time() - t0
    (OUT / 'sync_posthoc.json').write_text(json.dumps(res, **DUMP))


if __name__ == '__main__' and sys.argv[1] == 'sync_posthoc':
    sync_posthoc()
