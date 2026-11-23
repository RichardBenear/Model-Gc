"""EX-GC34H (exploratory): constraint-reduced analysis of the 34G cube-mean kinetic term B. Plan: spec/EX-GC34H_plan.md.

  python src/exgc34h.py reduced -> results/EX-GC34H/reduced.json  (H1-H4, H6, C1, S1)
  python src/exgc34h.py evolve  -> results/EX-GC34H/evolve.json   (H5)
Reuses src/exgc34e.py, exgc34f.py, exgc34g.py unchanged. Convention e^{-i w t}; zero-shift gauge.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import json, platform, sys, time
from pathlib import Path
import numpy as np
import scipy.linalg as sl
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import exgc34d as D
import exgc34e as E
import exgc34f as F
import exgc34g as G

OUT = ROOT / 'results/EX-GC34H'
H = E.H; NV = 11
IN, IB, IE, RED = F.IN, F.IB, F.IE, F.RED
DUMP = G.DUMP


def cases_def():
    reg = json.loads((ROOT / 'results/EX-GC34F/spectrum.json').read_text())
    dirs = dict(reg['directions']); appr = dict(reg['approach'])
    alld = dict(dirs); alld.update(appr)
    alld.update({'100': (1, 0, 0), 'd01+': (1, 1, 0), 'b1++': (1, 1, 1)} if '100' not in alld else {})
    out = [(n, kh) for kh in (0.0125, 0.025, 0.05, 0.1, 0.2) for n in dirs] + [(n, 0.05) for n in appr]
    for n in ('100', 'd01+', 'b1++'):
        for kh in (0.5, 1.0, 2.0, 3.0):
            out.append((n, kh))
    return out, alld


def qvec(alld, name, kh):
    return kh * G.unit(alld[name])


# ---------------- reduction ----------------
def reduce_B(P, q):
    C, B, A = P['C'], P['B'], P['A']
    k = np.linalg.norm(q) / H
    Aee = A[np.ix_(IE, IE)]
    U, s, Vh = np.linalg.svd(Aee)
    u = Vh[-1].conj(); Q = Vh[:-1].conj().T                        # null vector, complement basis (7x6)
    facts = dict(A_sv_rel=[float(v / s[0]) for v in s], B_red_max=float(np.abs(B[np.ix_(RED, RED)]).max() / np.abs(B).max()),
                 A_nrow_max=float(np.abs(A[IN]).max()), lapse_u_coupling=float(abs(C[IN, IE] @ u) / np.linalg.norm(C[IN, IE])))
    Mc = np.array([[C[IN, IN], C[IN, IE] @ u], [u.conj() @ C[np.ix_(IE, [IN])][:, 0], u.conj() @ C[np.ix_(IE, IE)] @ u]])
    Rc = np.vstack([C[IN, IE] @ Q, u.conj() @ C[np.ix_(IE, IE)] @ Q])
    sv = np.linalg.svd(Mc, compute_uv=False)
    facts['block_sv_ratio'] = float(sv[-1] / sv[0])
    if sv[-1] / sv[0] < 1e-10:
        return None, facts
    S = -np.linalg.solve(Mc, Rc)                                     # (n, a) = S y
    X = np.zeros((NV, 6), complex); X[IN] = S[0]; X[np.ix_(IE, range(6))] = Q + np.outer(u, S[1])
    Fm = np.linalg.solve(Q.conj().T @ Aee @ Q, Q.conj().T @ C[IE] @ X)
    return dict(F=Fm, X=X, u=u, Q=Q, S=S, k=k), facts


def full_resid(P, w, x):
    M = F.Mw(P, w)
    return float(np.linalg.norm(M @ x) / (np.linalg.norm(M, 2) * np.linalg.norm(x)))


def eig_reduced(P, R):
    lam, Y = np.linalg.eig(R['F'])
    nF = np.linalg.norm(R['F'], 2); k = R['k']
    zero = np.abs(lam) <= 1e-6 * nF
    roots = []
    for j in range(6):
        if zero[j]:
            continue
        w = np.sqrt(-lam[j] + 0j)
        if w.real < 0 or (abs(w.real) < 1e-14 and w.imag < 0):
            w = -w
        x = R['X'] @ Y[:, j]
        res = full_resid(P, w, x)
        res_red = float(np.linalg.norm((F.Mw(P, w) @ x)[RED]) / (np.linalg.norm(F.Mw(P, w), 2) * np.linalg.norm(x)))
        e = x[IE]
        roots.append(dict(w_over_k=[float(w.real / k), float(w.imag / k)], wh=float(abs(w) * H), resid11=res, resid_red=res_red,
                          physical=res <= 1e-8, growing=bool(abs(w.imag) > 1e-6 * abs(w)),
                          u_overlap=float(abs(R['u'].conj() @ e) / np.linalg.norm(e)), _y=Y[:, j], _w=w))
    return roots, int(zero.sum()), nF


# ---------------- Riemann with time dependence ----------------
def riemann_t(h, hd, hdd, q):
    """Linearised Riemann (a<b, c<d) of h(t,x) = h(t) e^{iq.x/H}; h, hd, hdd: h and its first two t-derivatives."""
    kk = 1j * np.r_[0.0, q / H]

    def dd(a, b, m, n):
        if a == 0 and b == 0:
            return hdd[m, n]
        if a == 0 or b == 0:
            c = b if a == 0 else a
            return kk[c] * hd[m, n]
        return kk[a] * kk[b] * h[m, n]
    out = []
    for a in range(4):
        for b in range(a + 1, 4):
            for c in range(4):
                for d in range(c + 1, 4):
                    out.append(0.5 * (dd(b, c, a, d) + dd(a, d, b, c) - dd(a, c, b, d) - dd(b, d, a, c)))
    return np.array(out)


def met(x, q):
    return E.metric_from_x(x, q)


def zero_sector(P, R, q):
    """H4: invariant subspace of L for the zero cluster, momentum-constraint filter, Riemann measure."""
    Fm = R['F']; k = R['k']; nF = np.linalg.norm(Fm, 2)
    L = np.block([[np.zeros((6, 6)), np.eye(6)], [Fm, np.zeros((6, 6))]])
    tol = np.sqrt(1e-6 * nF)
    T, Zs, sdim = sl.schur(L.astype(complex), output='complex', sort=lambda z: abs(z) <= tol)
    Z = Zs[:, :sdim]
    out = dict(dim_cluster_L=int(sdim))
    if sdim == 0:
        out.update(n_kept=0, r=[], r_max=0.0); return out
    N = Z.conj().T @ L @ Z
    out['N_norm_over_k'] = float(np.linalg.norm(N, 2) / k)
    ts = (0.0, 1.0 / k, 2.0 / k)
    X = R['X']; CB = P['C'][IB]; BB = P['B'][IB]; AB = P['A'][IB]
    blocks = []
    zt = {}
    for t in ts:
        Et = Z @ sl.expm(N * t)                                      # z(t) = Et c
        y = Et[:6]; yd = Et[6:]; ydd = Fm @ y
        zt[t] = (y, yd, ydd)
        blocks.append(CB @ X @ y + 1j * BB @ X @ yd - AB @ X @ ydd)
    Mm = np.vstack(blocks)
    scale = np.linalg.norm(CB @ X, 2) + k * np.linalg.norm(BB @ X, 2) + k * k * np.linalg.norm(AB @ X, 2)
    Uu, sv, Vh = np.linalg.svd(Mm)
    sv_full = np.zeros(sdim); sv_full[:len(sv)] = sv
    keep = Vh.conj().T[:, sv_full <= 1e-8 * scale]
    out['mom_sv_rel'] = [float(v / scale) for v in sv_full]
    out['n_kept'] = int(keep.shape[1])
    rs = []; bodies = []
    for j in range(keep.shape[1]):
        c = keep[:, j]
        Rm = []; hm = []; bd = []
        for t in ts[:2]:
            y, yd, ydd = [m @ c for m in zt[t]]
            x, xd, xdd = X @ y, X @ yd, X @ ydd
            h, hd_, hdd_ = met(x, q), met(xd, q), met(xdd, q)
            Rm.append(np.linalg.norm(riemann_t(h, hd_, hdd_, q))); hm.append(np.linalg.norm(h))
            x8 = np.concatenate([x[:1], x[4:]])
            _, hh, resid = D.metric_from_edges(x8, q)
            bd.append(abs(resid) / max(np.linalg.norm(hh), 1e-300))
        rs.append(float(max(Rm) / (k * k * max(hm)))); bodies.append(float(max(bd)))
    out['r'] = rs; out['r_max'] = max(rs) if rs else 0.0; out['body_resid_rel'] = bodies
    return out


def r_disp(q):
    k = np.linalg.norm(q) / H
    G0 = E.gauge_modes(q[None], 0.0)[0]
    vals = []
    for a in range(3):
        h = met(G0[:, a], q); z = np.zeros((4, 4), complex)
        vals.append(float(np.linalg.norm(riemann_t(h, z, z, q)) / (k * k * np.linalg.norm(h))))
    return max(vals)


def r_wave(R, rt, q):
    k = R['k']; w = rt['_w']; x = R['X'] @ rt['_y']
    h = met(x, q)
    return float(np.linalg.norm(riemann_t(h, -1j * w * h, -w * w * h, q)) / (k * k * np.linalg.norm(h)))


# ---------------- driver: H1-H4, H6, C1, S1 ----------------
def reduced():
    t0 = time.process_time()
    mp.mp.dps = 50                                     # 34D finite-difference Hessian needs high precision
    st = E.tet_static()
    KB = G.KernelG(st, G.kin_B); KT = G.KernelG(st, G.kin_T)
    cdef, alld = cases_def()
    g34 = json.loads((ROOT / 'results/EX-GC34G/candidates.json').read_text())['candidates']
    gB = {(c['dir'], c['kh']): c for c in g34['B']['cases']}
    cases = []
    for name, kh in cdef:
        q = qvec(alld, name, kh)
        P = KB.pencil(q)
        R, facts = reduce_B(P, q)
        c = dict(dir=name, kh=kh, light_speed=float(E.light_speed(q)), **facts)
        sA = facts['A_sv_rel']
        c['H1'] = bool(sA[-1] < 1e-12 and sA[-2] > 1e-6)
        c['H2'] = R is not None
        if R is None:
            cases.append(c); continue
        roots, nzero, nF = eig_reduced(P, R)
        c['n_zero_cluster'] = nzero
        c['roots'] = [{kk: v for kk, v in r.items() if not kk.startswith('_')} for r in roots]
        grow = [r for r in roots if r['growing']]
        phys = [r for r in roots if r['physical']]
        low = [r for r in phys if not r['growing'] and r['wh'] < 0.5]
        other_bad = [r for r in phys if r['wh'] < 0.5 and r['growing']]
        c['n_growing'] = len(grow); c['n_physical'] = len(phys); c['n_low_physical'] = len(low)
        c['max_resid11_physical'] = max([r['resid11'] for r in phys] or [0.0])
        c['H3'] = (len(grow) == 0) and (kh > 0.05 or (len(low) == 2 and not other_bad))
        c['H3_internal'] = all(r['resid11'] <= 1e-10 for r in phys)
        zs = zero_sector(P, R, q)
        c['zero'] = zs
        c['r_disp'] = r_disp(q)
        c['r_wave'] = [r_wave(R, r, q) for r in low]
        c['H4'] = all(v <= max(10 * c['r_disp'], 1e-10) for v in zs['r'])
        c['wave_speed_err'] = [r['w_over_k'][0] - c['light_speed'] for r in low]
        c['wave_u_overlap'] = [r['u_overlap'] for r in low]
        if kh <= 0.05:
            c['H6'] = len(low) == 2 and all(abs(v) <= 1e-3 for v in c['wave_speed_err'])
        # u at small kh (phase fixed: largest component real positive)
        uu = R['u'] * np.exp(-1j * np.angle(R['u'][np.argmax(np.abs(R['u']))]))
        c['u'] = [[float(v.real), float(v.imag)] for v in uu]
        # S1: match to 34G confirmed B roots
        if (name, kh) in gB:
            ref = [r['w_over_k'] for r in gB[(name, kh)]['roots']]
            mism = []
            for r in phys:
                wk = complex(*r['w_over_k'])
                d = min(abs(wk - complex(*rr)) for rr in ref) / abs(wk) if ref else np.inf
                mism.append(float(d))
            c['S1_rel'] = mism
            c['S1'] = all(m <= 1e-6 for m in mism)
        cases.append(c)
    res = dict(convention='e^{-i w t}; zero-shift gauge; (n, a) eliminated algebraically', cases=cases)
    def allp(key, sub=None):
        cc = [c for c in cases if key in c and (sub is None or sub(c))]
        return dict(n=len(cc), n_fail=sum(not c[key] for c in cc), fail=[(c['dir'], c['kh']) for c in cc if not c[key]][:40])
    res['H1'] = allp('H1'); res['H2'] = allp('H2'); res['H3'] = allp('H3'); res['H3_internal'] = allp('H3_internal')
    res['H4'] = allp('H4'); res['H6'] = allp('H6'); res['S1'] = allp('S1')
    # C1 control: T in the 34F reduction, generic kh = 0.05
    q = qvec(alld, 'generic', 0.05); PT = KT.pencil(q)
    rsT = F.reduced_system(PT, q)
    lam = np.linalg.eigvals(rsT['F']); k = np.linalg.norm(q) / H
    gT = max(np.sqrt(-lam + 0j).imag.max(), (-np.sqrt(-lam + 0j)).imag.max()) / k
    refT = [r['w_over_k'][1] for r in next(c for c in g34['T']['cases'] if c['dir'] == 'generic' and c['kh'] == 0.05)['roots']
            if r['cls'] == 'growing']
    res['C1'] = dict(growth_over_k=float(gT), ref=refT, rel=float(abs(gT - refT[0]) / refT[0]), pass_=bool(abs(gT - refT[0]) / refT[0] <= 1e-6))
    # wave speed table and fit
    tab = {}
    for c in cases:
        if 'wave_speed_err' in c and c['dir'] in ('100', 'd01+', 'b1++', 'generic', 'rand00', 'rand01'):
            tab.setdefault(c['dir'], []).append((c['kh'], c['wave_speed_err']))
    res['wave_speed_table'] = tab
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, scipy=__import__('scipy').__version__, system=platform.system())
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'reduced.json').write_text(json.dumps(res, **DUMP))
    for key in ('H1', 'H2', 'H3', 'H3_internal', 'H4', 'H6', 'S1'):
        print(key, res[key]['n'], 'fail', res[key]['n_fail'], res[key]['fail'][:6])
    print('C1', res['C1']); print('cpu', res['cpu_s'])


# ---------------- H5 ----------------
def evolve():
    t0 = time.process_time()
    mp.mp.dps = 50
    st = E.tet_static(); KB = G.KernelG(st, G.kin_B)
    _, alld = cases_def()
    out = dict(runs=[])
    for name in ('generic', 'd01+'):
        q = qvec(alld, name, 0.05); P = KB.pencil(q)
        R, _ = reduce_B(P, q); k = R['k']; Fm = R['F']; X = R['X']
        L = np.block([[np.zeros((6, 6)), np.eye(6)], [Fm, np.zeros((6, 6))]])
        roots, _, nF = eig_reduced(P, R)
        low = sorted([r for r in roots if r['physical'] and not r['growing'] and r['wh'] < 0.5], key=lambda r: r['_w'].real)
        wv = low[0]; w = wv['_w']
        zw = np.r_[wv['_y'], -1j * w * wv['_y']]; zw /= np.linalg.norm(zw)
        # seeds: every eigenvector of L outside the zero cluster except the wave itself, plus the zero-cluster
        # invariant-subspace basis (ordered Schur), each unit-normalised and scaled by 1e-8
        tol = np.sqrt(1e-6 * nF)
        _, Zs, sdim = sl.schur(L.astype(complex), output='complex', sort=lambda z: abs(z) <= tol)
        mu, V = np.linalg.eig(L)
        vecs = [Zs[:, j] for j in range(sdim)]
        for j in range(12):
            if abs(mu[j]) > tol and abs(abs(zw.conj() @ V[:, j]) / np.linalg.norm(V[:, j]) - 1) > 1e-8:
                vecs.append(V[:, j] / np.linalg.norm(V[:, j]))
        seeds = 1e-8 * np.sum(vecs, axis=0)
        wmax = np.abs(mu).max()
        out.setdefault('n_seeds', {})[name] = len(vecs)
        T = 40 * 2 * np.pi / w.real
        CB = P['C'][IB]; BB = P['B'][IB]; AB = P['A'][IB]
        scale = np.linalg.norm(CB @ X, 2) + k * np.linalg.norm(BB @ X, 2) + k * k * np.linalg.norm(AB @ X, 2)

        def tid(z):
            y = z[:6]; ydd = Fm @ y
            x = X @ y; xdd = X @ ydd
            h00 = met(x, q)[0, 0]; hdd = met(xdd, q)[1:, 1:]
            return 0.5 * (hdd - np.outer(q, q) / H ** 2 * h00)

        def mom(z):
            y = z[:6]; yd = z[6:]
            r = CB @ X @ y + 1j * BB @ X @ yd - AB @ X @ (Fm @ y)
            return float(np.linalg.norm(r) / (scale * np.linalg.norm(z)))
        run = dict(dir=name, w_over_k=float(w.real / k), wmax=float(wmax), T=float(T), steps={})
        for lab, dt in (('dt', 0.02 / wmax), ('dt2', 0.01 / wmax)):
            I = np.eye(12); Ms = np.linalg.solve(I - 0.5 * dt * L, I + 0.5 * dt * L)
            n = int(np.ceil(T / dt)); every = max(1, n // 4000)
            z1 = zw.astype(complex).copy(); z2 = (zw + seeds).astype(complex)
            Mev = np.linalg.matrix_power(Ms, every)
            Ew = []; Ed = []; mm = 0.0; amp = []
            for i in range(0, n + 1, every):
                e1 = tid(z1); e2 = tid(z2)
                Ew.append(np.linalg.norm(e1)); Ed.append(np.linalg.norm(e2 - e1)); mm = max(mm, mom(z2))
                amp.append(abs(zw.conj() @ z1))
                z1 = Mev @ z1; z2 = Mev @ z2
            meas = float(max(Ed) / max(Ew))
            run['steps'][lab] = dict(dt=dt, n=n, sample_every=every, measure=meas, mom_max=mm,
                                     wave_amp_change=float(max(abs(a - amp[0]) for a in amp) / amp[0]),
                                     Ed_first_last=[float(Ed[0]), float(Ed[-1])], Ew_max=float(max(Ew)))
        a, b = run['steps']['dt']['measure'], run['steps']['dt2']['measure']
        run['H5'] = bool(all(s['mom_max'] <= 1e-8 and s['measure'] <= 1e-6 for s in run['steps'].values())
                         and max(a, b) <= 2 * max(min(a, b), 1e-300))
        out['runs'].append(run)
        print(name, {l: (s['measure'], s['mom_max'], s['wave_amp_change']) for l, s in run['steps'].items()}, 'H5', run['H5'])
    out['H5'] = all(r['H5'] for r in out['runs'])
    out['cpu_s'] = time.process_time() - t0
    out['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    (OUT / 'evolve.json').write_text(json.dumps(out, **DUMP))
    print('cpu', out['cpu_s'])


if __name__ == '__main__' and sys.argv[1] in ('reduced', 'evolve'):
    dict(reduced=reduced, evolve=evolve)[sys.argv[1]]()


# ---------------- POST HOC (labelled; not part of the registered scoring) ----------------
def reduce_B_posthoc(P, q):
    """No block threshold. Off-axis: (n, a) from the 2x2 block as registered. Where the lapse does not couple to u
    (|C_n,e u| <= 1e-8 |C_n,e|): a from the u-row, n from the twice-differentiated lapse constraint (34F style)."""
    C, A = P['C'], P['A']; k = np.linalg.norm(q) / H
    Aee = A[np.ix_(IE, IE)]
    _, s, Vh = np.linalg.svd(Aee)
    u = Vh[-1].conj(); Q = Vh[:-1].conj().T
    Cne = C[IN, IE]; Cen = C[np.ix_(IE, [IN])][:, 0]; Cee = C[np.ix_(IE, IE)]
    coup = abs(Cne @ u) / np.linalg.norm(Cne)
    M6 = Q.conj().T @ Aee @ Q
    if coup > 1e-8:
        Mc = np.array([[C[IN, IN], Cne @ u], [u.conj() @ Cen, u.conj() @ Cee @ u]])
        Rc = np.vstack([Cne @ Q, u.conj() @ Cee @ Q])
        S = -np.linalg.solve(Mc, Rc)
        X = np.zeros((NV, 6), complex); X[IN] = S[0]; X[np.ix_(IE, range(6))] = Q + np.outer(u, S[1])
        Fm = np.linalg.solve(M6, Q.conj().T @ C[IE] @ X)
        mode = 'block'
    else:
        uCu = u.conj() @ Cee @ u
        a_n = -(u.conj() @ Cen) / uCu; a_y = -(u.conj() @ Cee @ Q) / uCu       # a = a_n n + a_y y
        Ey = Q + np.outer(u, a_y); En = Cen * 0 + u * a_n                          # edges = Ey y + En n
        Fy = np.linalg.solve(M6, Q.conj().T @ Cee @ Ey)
        fn = np.linalg.solve(M6, Q.conj().T @ (Cee @ En + Cen))
        cq = Cne @ Q
        Pi = -(cq @ Fy) / (cq @ fn)
        X = np.zeros((NV, 6), complex); X[IN] = Pi; X[np.ix_(IE, range(6))] = Ey + np.outer(En, Pi)
        Fm = Fy + np.outer(fn, Pi)
        mode = 'differentiated_lapse'
    return dict(F=Fm, X=X, u=u, Q=Q, k=k), dict(mode=mode, lapse_u_coupling=float(coup))


def eig_posthoc(P, R):
    """Post hoc: zero cluster |lam| <= 1e-2 k^2; a root whose vector is still annihilated at 1.01 w is a direction of
    the singular pencil (normal-rank null space), not a root."""
    lam, Y = np.linalg.eig(R['F']); k = R['k']
    zero = np.abs(lam) <= 1e-2 * k * k
    roots = []
    for j in range(6):
        if zero[j]:
            continue
        w = np.sqrt(-lam[j] + 0j)
        if w.real < 0 or (abs(w.real) < 1e-14 and w.imag < 0):
            w = -w
        x = R['X'] @ Y[:, j]; e = x[IE]
        res = full_resid(P, w, x); res_pert = full_resid(P, 1.01 * w + 0.01j * abs(w), x)
        roots.append(dict(w_over_k=[float(w.real / k), float(w.imag / k)], wh=float(abs(w) * H), resid11=res, resid_perturbed=res_pert,
                          pencil_null=bool(res_pert <= 1e-8), physical=bool(res <= 1e-8 and res_pert > 1e-8),
                          growing=bool(abs(w.imag) > 1e-6 * abs(w)), u_overlap=float(abs(R['u'].conj() @ e) / np.linalg.norm(e)),
                          _y=Y[:, j], _w=w))
    return roots, int(zero.sum())


def zero_sector_posthoc(P, R, q):
    """Zero cluster; filter on ALL 11 rows of C x + iB x' - A x'' at t = 0, 1/k, 2/k; classify static vs drifting."""
    Fm = R['F']; k = R['k']; nF = np.linalg.norm(Fm, 2)
    L = np.block([[np.zeros((6, 6)), np.eye(6)], [Fm, np.zeros((6, 6))]])
    tol = 0.1 * k
    _, Zs, sdim = sl.schur(L.astype(complex), output='complex', sort=lambda z: abs(z) <= tol)
    Z = Zs[:, :sdim]; N = Z.conj().T @ L @ Z
    ts = (0.0, 1.0 / k, 2.0 / k); X = R['X']; C, B, A = P['C'], P['B'], P['A']
    zt = {}; blocks = []
    for t in ts:
        Et = Z @ sl.expm(N * t); y = Et[:6]; yd = Et[6:]; ydd = Fm @ y
        zt[t] = (y, yd, ydd)
        blocks.append(C @ X @ y + 1j * B @ X @ yd - A @ X @ ydd)
    Mm = np.vstack(blocks)
    scale = np.linalg.norm(C @ X, 2) + k * np.linalg.norm(B @ X, 2) + k * k * np.linalg.norm(A @ X, 2)
    _, sv, Vh = np.linalg.svd(Mm)
    svf = np.zeros(sdim); svf[:len(sv)] = sv
    keep = Vh.conj().T[:, svf <= 1e-8 * scale]
    sols = []
    G0 = E.gauge_modes(q[None], 0.0)[0]
    for j in range(keep.shape[1]):
        c = keep[:, j]; Rm = []; hm = []; bd = []; xs = []
        for t in ts[:2]:
            y, yd, ydd = [m @ c for m in zt[t]]
            x, xd, xdd = X @ y, X @ yd, X @ ydd; xs.append(x)
            h = met(x, q); Rm.append(np.linalg.norm(riemann_t(h, met(xd, q), met(xdd, q), q))); hm.append(np.linalg.norm(h))
            _, hh, resid = D.metric_from_edges(np.concatenate([x[:1], x[4:]]), q)
            bd.append(abs(resid) / max(np.linalg.norm(hh), 1e-300))
        y0, yd0, _ = [m @ c for m in zt[0.0]]
        x0 = X @ y0
        # overlap with exact displacement gauge (static)
        cdisp, *_ = np.linalg.lstsq(G0, x0, rcond=None)
        sols.append(dict(r=float(max(Rm) / (k * k * max(hm))), body_resid_rel=float(max(bd)),
                         drift=float(np.linalg.norm(X @ yd0) / (k * np.linalg.norm(x0))),
                         disp_frac=float(np.linalg.norm(G0 @ cdisp) / np.linalg.norm(x0)),
                         growth_2_over_k=float(np.linalg.norm(zt[2.0 / k][0] @ c) / np.linalg.norm(y0))))
    return dict(dim_cluster_L=int(sdim), n_kept=len(sols), sols=sols)


def posthoc():
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); KB = G.KernelG(st, G.kin_B)
    cdef, alld = cases_def()
    g34 = json.loads((ROOT / 'results/EX-GC34G/candidates.json').read_text())['candidates']
    gB = {(c['dir'], c['kh']): c for c in g34['B']['cases']}
    cases = []
    for name, kh in cdef:
        q = qvec(alld, name, kh); P = KB.pencil(q)
        R, facts = reduce_B_posthoc(P, q)
        roots, nzero = eig_posthoc(P, R)
        phys = [r for r in roots if r['physical']]
        low = [r for r in phys if not r['growing'] and r['wh'] < 0.5]
        c = dict(dir=name, kh=kh, light_speed=float(E.light_speed(q)), **facts, n_zero_cluster=nzero,
                 n_pencil_null=sum(r['pencil_null'] for r in roots),
                 n_growing=sum(r['growing'] for r in roots if not r['pencil_null']), n_physical=len(phys), n_low_physical=len(low),
                 max_resid11_physical=max([r['resid11'] for r in phys] or [0.0]),
                 roots=[{kk: v for kk, v in r.items() if not kk.startswith('_')} for r in roots])
        c['H3'] = c['n_growing'] == 0 and (kh > 0.05 or len(low) == 2)
        c['wave_speed_err'] = [r['w_over_k'][0] - c['light_speed'] for r in low]
        c['wave_u_overlap'] = [r['u_overlap'] for r in low]
        c['zero'] = zero_sector_posthoc(P, R, q); c['r_disp'] = r_disp(q)
        if (name, kh) in gB:
            ref = [complex(*r['w_over_k']) for r in gB[(name, kh)]['roots']]
            c['S1_rel'] = [float(min([abs(complex(*r['w_over_k']) - z) for z in ref] or [np.inf]) / abs(complex(*r['w_over_k']))) for r in phys]
            # 34G low propagating roots that the reduction does not reproduce
            mine = [complex(*r['w_over_k']) for r in phys]
            c['g34_low_unmatched'] = [[z.real, z.imag] for z in ref if 0 < z.real and abs(z.imag) < 1e-6 * abs(z)
                                      and z.real * kh < 0.5 and min([abs(z - m) for m in mine] or [np.inf]) > 1e-6 * abs(z)]
        cases.append(c)
    res = dict(label='POST HOC (after the registered run): no block threshold; differentiated lapse where the lapse does not couple to u; zero cluster |lam| <= 1e-2 k^2 (Schur |mu| <= 0.1 k); singular-pencil directions excluded; zero-sector filter on all 11 rows',
               cases=cases)
    res['summary'] = dict(
        modes={m: sum(c['mode'] == m for c in cases) for m in ('block', 'differentiated_lapse')},
        H3_fail=[(c['dir'], c['kh']) for c in cases if not c['H3']],
        n_growing_total=sum(c['n_growing'] for c in cases),
        max_resid11=max(c['max_resid11_physical'] for c in cases),
        max_abs_speed_err_kh_le_005=max(abs(v) for c in cases if c['kh'] <= 0.05 for v in c['wave_speed_err']),
        speed_fail_kh_le_005=[(c['dir'], c['kh'], c['wave_speed_err']) for c in cases if c['kh'] <= 0.05 and any(abs(v) > 1e-3 for v in c['wave_speed_err'])],
        S1_max=max(v for c in cases for v in c.get('S1_rel', [0.0]) if c['kh'] <= 0.1),
        g34_low_unmatched=[(c['dir'], c['kh'], c['g34_low_unmatched']) for c in cases if c.get('g34_low_unmatched')])
    res['cpu_s'] = time.process_time() - t0
    (OUT / 'reduced_posthoc.json').write_text(json.dumps(res, **DUMP))
    s = res['summary']
    for kk, v in s.items():
        print(kk, v if not isinstance(v, list) else (len(v), v[:8]))


if __name__ == '__main__' and sys.argv[1] == 'posthoc':
    posthoc()
