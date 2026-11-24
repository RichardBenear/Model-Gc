"""EX-GC34N (exploratory): Y1, a lattice staggered by index parity (n, phi, xi, h_aa at sites; beta, u on links;
h_ab on plaquettes). Vacuum = linearised Einstein-Hilbert (Fierz-Pauli) with staggered differences, k -> K =
2 sin(q/2)/h. Plan: spec/EX-GC34N_plan.md.

  python src/exgc34n.py build   -> results/EX-GC34N/build.json   (A: I0, G7)
  python src/exgc34n.py spec    -> results/EX-GC34N/spec.json    (B: G1, G2, G3)
  python src/exgc34n.py static  -> results/EX-GC34N/static.json  (C: gamma, E4)
  python src/exgc34n.py sync    -> results/EX-GC34N/sync.json    (D: G4'')
  python src/exgc34n.py matter  -> results/EX-GC34N/matter.json  (E: M0, M1, M2)
Reuses exgc34e.gr_solution / riemann (continuum GR), exgc34m.packet, the 34F direction set.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import itertools, json, platform, sys, time
from pathlib import Path
import numpy as np
import scipy.linalg as sl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import exgc34e as E
import exgc34m as Mm

OUT = ROOT / 'results/EX-GC34N'
H = E.H; KAP = E.KAPPA; ETA = np.diag([-1.0, 1, 1, 1])
NX = 10
PAIRS = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)]      # h components after n, beta
DUMP = dict(indent=1, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else (o.item() if hasattr(o, 'item') else str(o)))
OFFS = [np.zeros(3)] + [np.eye(3)[a] / 2 for a in range(3)] + \
       [(np.eye(3)[a] + np.eye(3)[b]) / 2 if a != b else np.zeros(3) for a, b in PAIRS]


def Kvec(q):
    return 2 * np.sin(np.asarray(q, float) / 2) / H


def Pmat():
    """x (10) -> h_{mu nu} flattened (16)."""
    P = np.zeros((16, NX))
    P[0, 0] = -2
    for a in range(3):
        P[0 * 4 + 1 + a, 1 + a] = P[(1 + a) * 4 + 0, 1 + a] = 1
    for c, (a, b) in enumerate(PAIRS):
        P[(1 + a) * 4 + 1 + b, 4 + c] = 1; P[(1 + b) * 4 + 1 + a, 4 + c] = 1
    return P


PM = Pmat()


def G_up(h, kl):
    """Linearised Einstein tensor with indices raised, for lower wavevector kl (4)."""
    ku = ETA @ kl; k2 = kl @ ku; tr = np.trace(ETA @ h); kh = h @ ku; khk = ku @ h @ ku
    G = 0.5 * (k2 * h - np.outer(kl, kh) - np.outer(kh, kl) + np.outer(kl, kl) * tr + ETA * khk - ETA * k2 * tr)
    return ETA @ G @ ETA


def kernel(q, w):
    """10x10 kernel, amplitudes referred to each field's own position: M x = kappa P^T T^{up}."""
    kl = np.array([-w, *Kvec(q)])
    Gm = np.zeros((16, NX), complex)
    for j in range(NX):
        Gm[:, j] = G_up((PM[:, j]).reshape(4, 4).astype(complex), kl).reshape(16)
    return PM.T @ Gm


def gauge(q, w):
    K = Kvec(q); W = np.zeros((NX, 4), complex)
    W[0, 0] = -1j * w; W[1:4, 0] = -1j * K                    # re-slicing
    for a in range(3):
        W[1 + a, 1 + a] = -1j * w
        for c, (b, d) in enumerate(PAIRS):
            if b == d == a:
                W[4 + c, 1 + a] = 2j * K[a]
            elif b != d and a in (b, d):
                W[4 + c, 1 + a] = 1j * K[d if a == b else b]
    return W


def site_kernel(q, w):
    ph = np.exp(1j * np.array([o @ np.asarray(q, float) for o in OFFS]))
    return (ph.conj()[:, None]) * kernel(q, w) * ph[None, :]


# ---------------- A ----------------
def build():
    t0 = time.process_time(); rng = np.random.default_rng(3434 + 14); res = {}
    i0 = []
    for _ in range(20):
        q = rng.uniform(-np.pi, np.pi, 3); w = rng.uniform(-3, 3) / H
        Mq = kernel(q, w); W = gauge(q, w)
        i0.append(dict(null=float(np.linalg.norm(Mq @ W) / (np.linalg.norm(Mq) * np.linalg.norm(W))),
                       herm=float(np.abs(Mq - Mq.conj().T).max() / np.abs(Mq).max())))
    res['I0_cases'] = i0
    res['I0'] = bool(max(c['null'] for c in i0) <= 1e-12 and max(c['herm'] for c in i0) <= 1e-14)
    # G7: periodicity and stencil range
    per = []
    for _ in range(10):
        q = rng.uniform(-np.pi, np.pi, 3); w = rng.uniform(-3, 3) / H
        for a in range(3):
            q2 = q.copy(); q2[a] += 2 * np.pi
            per.append(float(np.abs(site_kernel(q2, w) - site_kernel(q, w)).max() / np.abs(site_kernel(q, w)).max()))
    N = 16; ks = 2 * np.pi * np.arange(N) / N; w = 0.37 / H
    S = np.zeros((N, N, N, NX, NX), complex)
    for i, j, k in itertools.product(range(N), repeat=3):
        S[i, j, k] = site_kernel(np.array([ks[i], ks[j], ks[k]]), w)
    R = np.fft.ifftn(S, axes=(0, 1, 2))
    idx = np.fft.fftfreq(N, 1 / N).astype(int)
    rng_max = 0; outside = 0.0; tot = np.abs(R).max()
    for i, j, k in itertools.product(range(N), repeat=3):
        r = max(abs(idx[i]), abs(idx[j]), abs(idx[k])); m = np.abs(R[i, j, k]).max()
        if m > 1e-13 * tot:
            rng_max = max(rng_max, r)
        if r > 2:
            outside = max(outside, m / tot)
    res['G7'] = dict(periodicity_max=max(per), stencil_range=int(rng_max), outside_range2=float(outside),
                     G7=bool(max(per) <= 1e-12 and outside <= 1e-13))
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'build.json').write_text(json.dumps(res, **DUMP))
    print('I0', res['I0'], max(c['null'] for c in i0), '| G7', res['G7'])


# ---------------- B ----------------
IH = list(range(4, 10)); ICON = [0, 1, 2, 3]


def abc(q):
    M0 = kernel(q, 0.0); Mp = kernel(q, 1.0); Mn = kernel(q, -1.0)
    return M0, (Mp - Mn) / 2, (Mp + Mn) / 2 - M0          # C, B, A


def sync_system(q):
    C, B, A = abc(q)
    Ahh = A[np.ix_(IH, IH)]; Chh = C[np.ix_(IH, IH)]; F = np.linalg.solve(Ahh, Chh)
    L = np.block([[np.zeros((6, 6)), np.eye(6)], [F, np.zeros((6, 6))]])
    Bh = B[np.ix_(IH, IH)]
    C0 = np.hstack([C[np.ix_(ICON, IH)] - A[np.ix_(ICON, IH)] @ F, 1j * B[np.ix_(ICON, IH)]])
    return L, F, C0, float(np.abs(Bh).max())


def light34a(q):
    return E.light_speed(np.asarray(q, float))


def spec():
    t0 = time.process_time(); res = dict(cases=[])
    reg = json.loads((ROOT / 'results/EX-GC34F/spectrum.json').read_text())
    dirs = dict(reg['directions']); dirs.update(reg['approach'])
    cases = [(n, kh) for kh in (0.0125, 0.025, 0.05, 0.1, 0.2) for n in dirs]
    for n, d in (('100', (1, 0, 0)), ('d01+', (1, 1, 0)), ('b1++', (1, 1, 1))):
        dd = np.array(d, float) / np.linalg.norm(d)
        for kh in (0.5, 1.0, 2.0, 3.0, float(np.pi / np.abs(dd).max())):
            cases.append((n, kh))
    for n, kh in cases:
        d = np.array(dirs.get(n, {'100': (1, 0, 0), 'd01+': (1, 1, 0), 'b1++': (1, 1, 1)}.get(n)), float)
        d = d / np.linalg.norm(d); q = kh * d; K = np.linalg.norm(Kvec(q)); k = kh / H
        L, F, C0, Bhh = sync_system(q)
        lam = np.linalg.eigvals(L)
        grow = float(np.abs(lam.real).max() / K)
        waves = sorted(abs(l.imag) / k for l in lam if abs(abs(l) / K - 1) < 1e-6 and abs(l.real) < 1e-8 * K)
        zero = int(sum(abs(l) <= 1e-6 * K for l in lam))
        other = [complex(l) / K for l in lam if not (abs(abs(l) / K - 1) < 1e-6) and abs(l) > 1e-6 * K]
        nul = int(sum(s <= 1e-10 * np.linalg.svd(kernel(q, 0.0), compute_uv=False)[0]
                      for s in np.linalg.svd(kernel(q, 0.0), compute_uv=False)))
        res['cases'].append(dict(dir=n, kh=kh, max_re_over_K=grow, n_wave_roots=len(waves), wave_speed=waves[:1],
                                 light34a=float(light34a(q)), n_zero=zero, other=[[o.real, o.imag] for o in other],
                                 static_nullity=nul, Bhh=Bhh))
    C = res['cases']
    res['G1'] = bool(max(c['max_re_over_K'] for c in C) <= 1e-8)
    g2 = [c for c in C if c['kh'] <= 0.05 + 1e-12]
    res['G2_max_dc'] = max(abs(c['wave_speed'][0] - c['light34a']) for c in g2 if c['wave_speed'])
    res['G2'] = bool(all(c['n_wave_roots'] == 4 and not c['other'] for c in C) and res['G2_max_dc'] <= 1e-3)
    res['G3'] = bool(all(c['static_nullity'] == 4 for c in C))
    res['n_cases'] = len(C); res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'spec.json').write_text(json.dumps(res, **DUMP))
    print('G1', res['G1'], max(c['max_re_over_K'] for c in C), '| G2', res['G2'], res['G2_max_dc'],
          'waves', set(c['n_wave_roots'] for c in C), 'other', max(len(c['other']) for c in C), 'zero', set(c['n_zero'] for c in C),
          '| G3', res['G3'], set(c['static_nullity'] for c in C), '| cases', len(C), 'cpu', res['cpu_s'])


# ---------------- C ----------------
def source(q, v, rho=1.0):
    k = np.asarray(q, float) / H; K = Kvec(q)
    rat = np.array([1.0 if abs(ka) < 1e-14 else ka / Ka for ka, Ka in zip(k, K)])
    p = rho * v * rat
    Tup = np.zeros((4, 4)); Tup[0, 0] = rho; Tup[0, 1:] = Tup[1:, 0] = p
    Tup[1:, 1:] = rho * np.outer(v * rat, v * rat)
    return KAP * PM.T @ Tup.reshape(16)


def lattice_riemann(q, v):
    w = (np.asarray(q, float) / H) @ v
    Mq = kernel(q, w); W = gauge(q, w); f = source(q, v)
    x, *_ = np.linalg.lstsq(Mq + W @ W.conj().T, f, rcond=1e-14)
    kl = np.array([-w, *Kvec(q)])
    return E.riemann((PM @ x).reshape(4, 4), kl), x, Mq, W, f, kl


def static():
    t0 = time.process_time(); res = {}
    gam = []
    for kh in (0.05, 0.5, 2.0):
        for d in ((1, 0, 0), (1, 1, 1), (0.8, 0.5, 0.33)):
            q = kh * np.array(d, float) / np.linalg.norm(d)
            Mq = kernel(q, 0.0); W = gauge(q, 0.0); f = source(q, np.zeros(3))
            x, *_ = np.linalg.lstsq(Mq + W @ W.conj().T, f, rcond=1e-14)
            K = Kvec(q); Kh = K / np.linalg.norm(K); Pt = np.eye(3) - np.outer(Kh, Kh)
            hm = (PM @ x).reshape(4, 4)[1:, 1:]
            A = 0.5 * np.einsum('ij,ij->', Pt, hm); g = -A / (2 * x[0])
            gam.append(dict(kh=kh, dir=list(d), gamma=[g.real, g.imag], resid=float(np.linalg.norm(Mq @ x - f) / np.linalg.norm(f))))
    res['gamma'] = gam
    res['G5_static'] = bool(max(abs(complex(*c['gamma']) - 1) for c in gam) <= 1e-10)
    cases = []
    for qdir in ((1, 0, 0), (1, 1, 1), (0.8, 0.5, 0.33)):
        qd = np.array(qdir, float) / np.linalg.norm(qdir)
        perp = np.cross(qd, [0, 0, 1.0]) if abs(qd[2]) < 0.9 else np.cross(qd, [1.0, 0, 0])
        for vn, vd in (('par', qd), ('perp', perp), ('obl', qd + 1.3 * perp)):
            cases.append((qdir, vn, np.array(vd) / np.linalg.norm(vd)))
    e4 = {}; inv = []; resid = []
    for kh in (0.0125, 0.025, 0.05, 0.1):
        errs = []
        for qdir, vn, vd in cases:
            qd = np.array(qdir, float) / np.linalg.norm(qdir); q = kh * qd; k = q / H; v = 1e-3 * vd
            R = {}
            for s in (1, -1):
                Rl, x, Mq, W, f, kl = lattice_riemann(q, s * v)
                hg, klg, _ = E.gr_solution(k, s * v); R[s] = (Rl, E.riemann(hg, klg))
                resid.append(float(np.linalg.norm(Mq @ x - f) / np.linalg.norm(f)))
                if s == 1 and kh == 0.05:
                    for j in range(4):
                        inv.append(float(np.abs(E.riemann((PM @ (x + W[:, j] * np.linalg.norm(x) / np.linalg.norm(W[:, j]))).reshape(4, 4), kl) - Rl).max()
                                         / np.abs(Rl).max()))
            R1l = (R[1][0] - R[-1][0]) / 2; R1g = (R[1][1] - R[-1][1]) / 2
            R0l = (R[1][0] + R[-1][0]) / 2; R0g = (R[1][1] + R[-1][1]) / 2
            errs.append(dict(qdir=list(qdir), v=vn, err_Ov=float(np.abs(R1l - R1g).max() / np.abs(R1g).max()),
                             err_static=float(np.abs(R0l - R0g).max() / np.abs(R0g).max())))
        e4[str(kh)] = dict(cases=errs, max=max(max(c['err_Ov'], c['err_static']) for c in errs))
        print('kh', kh, 'E4 %.3e' % e4[str(kh)]['max'], flush=True)
    res['E4'] = e4; res['E4_ratio'] = e4['0.05']['max'] / e4['0.025']['max']
    res['riemann_gauge_change_max'] = max(inv); res['forced_resid_max'] = max(resid)
    res['G5_E4'] = bool(e4['0.05']['max'] <= 1e-2 and res['E4_ratio'] >= 3)
    res['cpu_s'] = time.process_time() - t0
    (OUT / 'static.json').write_text(json.dumps(res, **DUMP))
    print('gamma', res['G5_static'], max(abs(complex(*c['gamma']) - 1) for c in gam), '| E4', res['G5_E4'],
          e4['0.05']['max'], 'ratio', res['E4_ratio'], '| gauge change', res['riemann_gauge_change_max'], 'resid', res['forced_resid_max'])


# ---------------- D ----------------
def sync():
    t0 = time.process_time(); res = dict(runs=[])
    for dn, d in (('generic', (0.8, 0.5, 0.33)), ('d01+', (1, 1, 0))):
        q = 0.05 * np.array(d) / np.linalg.norm(d); K = np.linalg.norm(Kvec(q))
        L, F, C0, _ = sync_system(q)
        Cmap = np.vstack([C0, C0 @ L]); nC = np.linalg.norm(Cmap, 2)
        mu_, V = np.linalg.eig(L); zc = np.abs(mu_) <= 0.01 * K
        waves = [j for j in range(12) if not zc[j] and abs(abs(mu_[j]) / K - 1) < 0.01 and mu_[j].imag < 0]
        jw = waves[0]; wv = abs(mu_[jw]); zw = V[:, jw] / np.linalg.norm(V[:, jw])
        nonz = [V[:, j] / np.linalg.norm(V[:, j]) for j in range(12) if not zc[j] and j != jw]
        nz_ok = [v for v in nonz if np.linalg.norm(Cmap @ v) <= 1e-8 * nC]
        _, Zs, sdim = sl.schur(L.astype(complex), output='complex', sort=lambda z: abs(z) <= 0.01 * K)
        Z = Zs[:, :sdim]; _, sz, Vz = np.linalg.svd(Cmap @ Z)
        rz = int((sz > sz[0] * 1e-10).sum()) if sz.size and sz[0] > 0 else 0
        Zk = Z @ Vz[rz:].conj().T
        sb = sum(nz_ok) if nz_ok else np.zeros(12); sc_ = Zk.sum(1) if Zk.shape[1] else np.zeros(12)
        wmax = np.abs(mu_).max(); T = 40 * 2 * np.pi / wv

        def tid(z):
            hh = F @ z[:6]; m = np.zeros((3, 3), complex)
            for c, (a, b) in enumerate(PAIRS):
                m[a, b] = m[b, a] = 0.5 * hh[c]
            return m

        def crel(z):
            return float(np.linalg.norm(C0 @ z) / (np.linalg.norm(C0, 2) * np.linalg.norm(z)))
        out = dict(dir=dn, wave_w_over_K=float(wv / K), n_nonzero_seeds=len(nz_ok), n_nonzero_total=len(nonz),
                   zero_cluster_dim=int(sdim), zero_kept=int(Zk.shape[1]))
        for tag, fac in (('dt', 0.02), ('dt2', 0.01)):
            ns = int(np.ceil(T / (fac / wmax))); dt = T / ns
            I12 = np.eye(12); Ms = np.linalg.solve(I12 - 0.5 * dt * L, I12 + 0.5 * dt * L)
            every = max(1, ns // 4000); Mev = np.linalg.matrix_power(Ms, every)
            za = zw.copy(); zb = zw + 1e-8 * sb; zz = zw + 1e-8 * sc_
            Ea0 = np.linalg.norm(tid(za)); Eam = 0.0; amp = 0.0; db = 0.0; dc = 0.0; cm = 0.0
            for i in range(0, ns + 1, every):
                Ea = tid(za); nE = np.linalg.norm(Ea); Eam = max(Eam, nE); amp = max(amp, abs(nE / Ea0 - 1))
                db = max(db, np.linalg.norm(tid(zb) - Ea)); dc = max(dc, np.linalg.norm(tid(zz) - Ea))
                cm = max(cm, crel(za), crel(zb), crel(zz))
                za = Mev @ za; zb = Mev @ zb; zz = Mev @ zz
            out[tag] = dict(dt=dt, steps=ns, max_constraint=cm, wave_tidal_amp_change=float(amp),
                            seeds_nonzero_rel=float(db / Eam), seeds_zero_rel=float(dc / Eam))
        a, b = out['dt'], out['dt2']
        out['G4pp'] = bool(all(t['max_constraint'] <= 1e-9 and t['wave_tidal_amp_change'] <= 1e-6 and t['seeds_nonzero_rel'] <= 1e-6
                               and t['seeds_zero_rel'] <= 1e-6 for t in (a, b)))
        res['runs'].append(out); print(out, flush=True)
    res['G4pp'] = all(r['G4pp'] for r in res['runs']); res['cpu_s'] = time.process_time() - t0
    (OUT / 'sync.json').write_text(json.dumps(res, **DUMP))
    print('G4pp', res['G4pp'])


# ---------------- E: matter on Y1 ----------------
MASS = 1.0


def fw(a, ax, s=1):
    return np.roll(a, -s, ax)          # a(v + s e_ax)


def grads(phi, h):
    return [(fw(phi, a) - phi) / h for a in range(3)]


def msources(phi, phid, h):
    """s_n (site), s_b[a] (link a, index v -> v+e_a/2), s_d[a] = dS/dh_aa (site), s_o[c] = dS/dh_ab (plaquette, PAIRS[3+c])."""
    m2 = MASS ** 2; V = h ** 3; G = grads(phi, h)
    s_n = V * (-0.5 * phid ** 2 - 0.5 * m2 * phi ** 2)
    s_d = [V * (0.25 * phid ** 2 - 0.25 * m2 * phi ** 2) for _ in range(3)]
    s_b = []
    for a in range(3):
        g2 = G[a] ** 2; g2s = g2 + fw(g2, a, -1)            # links v and v - e_a, both touching site v
        s_n = s_n - V * 0.25 * g2s
        for c in range(3):
            s_d[c] = s_d[c] - V * 0.125 * g2s             # -1/2 * 1/2 * T-bar
        s_d[a] = s_d[a] + V * 0.25 * g2s                  # + 1/2 h-bar_aa G_a^2
        s_b.append(-V * 0.5 * (phid + fw(phid, a)) * G[a])
    s_o = []
    for (a, b) in PAIRS[3:]:
        Ga = 0.5 * (G[a] + fw(G[a], b)); Gb = 0.5 * (G[b] + fw(G[b], a))
        s_o.append(V * Ga * Gb)
    return s_n, np.array(s_b), np.array(s_d), np.array(s_o)


def mlagrangian(phi, phid, n, b, hd, ho, h):
    m2 = MASS ** 2; V = h ** 3; G = grads(phi, h); T = hd.sum(0)
    L = np.sum(0.5 * (1 - n + 0.5 * T) * phid ** 2 - 0.5 * (1 + n + 0.5 * T) * m2 * phi ** 2)
    for a in range(3):
        nb = 0.5 * (n + fw(n, a)); Tb = 0.5 * (T + fw(T, a)); hb = 0.5 * (hd[a] + fw(hd[a], a))
        L += np.sum(-0.5 * (1 + nb + 0.5 * Tb - hb) * G[a] ** 2 - 0.5 * (phid + fw(phid, a)) * b[a] * G[a])
    for c, (a, bb) in enumerate(PAIRS[3:]):
        Ga = 0.5 * (G[a] + fw(G[a], bb)); Gb = 0.5 * (G[bb] + fw(G[bb], a))
        L += np.sum(ho[c] * Ga * Gb)
    return V * L


def lap(phi, h):
    return sum(np.roll(phi, 1, a) + np.roll(phi, -1, a) - 2 * phi for a in range(3)) / h ** 2


def msources_dot(phi, phid, h):
    phidd = -MASS ** 2 * phi + lap(phi, h)
    sp = msources(phi + phid, phid + phidd, h); sm = msources(phi - phid, phid - phidd, h)
    return [(x - y) / 2 for x, y in zip(sp, sm)]


def midentities(s, sd, h):
    s_n, s_b, s_d, s_o = s; sn_d, sb_d, _, _ = sd
    IE = -sn_d + sum(s_b[a] - fw(s_b[a], a, -1) for a in range(3)) / h
    div = np.zeros_like(s_b)
    for a in range(3):
        div[a] += (2 / h) * (s_d[a] - fw(s_d[a], a))
    for c, (a, b) in enumerate(PAIRS[3:]):
        div[a] += (fw(s_o[c], b, -1) - s_o[c]) / h          # d h_ab / d u_a: Delta_b u_a
        div[b] += (fw(s_o[c], a, -1) - s_o[c]) / h          # Delta_a u_b
    IP = -sb_d + div
    return IE, IP, div


def matter():
    t0 = time.process_time(); rng = np.random.default_rng(3434 + 15); res = {}
    m0 = []
    for h in (1.0, 0.5):
        N = 8; phi = rng.standard_normal((N,) * 3); phid = rng.standard_normal((N,) * 3)
        n = rng.standard_normal((N,) * 3); b = rng.standard_normal((3,) + (N,) * 3)
        hd = rng.standard_normal((3,) + (N,) * 3); ho = rng.standard_normal((3,) + (N,) * 3)
        s_n, s_b, s_d, s_o = msources(phi, phid, h)
        pred = np.sum(s_n * n) + np.sum(s_b * b) + np.sum(s_d * hd) + np.sum(s_o * ho)

        def cd(eps):
            return (mlagrangian(phi, phid, eps * n, eps * b, eps * hd, eps * ho, h)
                    - mlagrangian(phi, phid, -eps * n, -eps * b, -eps * hd, -eps * ho, h)) / (2 * eps)
        rich = (4 * cd(1e-4) - cd(2e-4)) / 3
        m0.append(dict(h=h, rel=float(abs(rich - pred) / abs(pred)), rel_plain=float(abs(cd(1e-5) - pred) / abs(pred))))
    res['M0_rows'] = m0; res['M0'] = bool(max(r['rel'] for r in m0) <= 1e-8)
    print('M0', res['M0'], m0, flush=True)
    rows = []
    for dname in Mm.DIRS:
        for h in (1.0, 0.5, 0.25):
            N = int(round(Mm.LBOX / h))
            for t in (0.0, 6.0, 12.0):
                phi, phid, KK = Mm.packet(N, h, dname, t)
                s = msources(phi, phid, h); sd = msources_dot(phi, phid, h)
                IE, IP, div = midentities(s, sd, h)
                rE = float(np.abs(IE).max() / np.abs(sd[0]).max())
                rP = float(np.abs(IP).max() / max(np.abs(sd[1]).max(), np.abs(div).max()))
                glob = float(np.abs(IP.reshape(3, -1).sum(1)).max() / np.abs(sd[1]).max())
                rows.append(dict(dir=dname, h=h, t=t, rE=rE, rP=rP, global_rel=glob))
                print(dname, h, t, 'rE %.2e rP %.3e glob %.1e cpu %.0f' % (rE, rP, glob, time.process_time() - t0), flush=True)
    res['rows'] = rows
    res['M1'] = bool(max(r['rE'] for r in rows) <= 1e-10)
    ratios = []
    for dname in Mm.DIRS:
        for t in (0.0, 6.0, 12.0):
            r = [next(x['rP'] for x in rows if x['dir'] == dname and x['h'] == h and x['t'] == t) for h in (1.0, 0.5, 0.25)]
            ratios.append(dict(dir=dname, t=t, r1=r[0] / r[1], r2=r[1] / r[2]))
    res['M2_ratios'] = ratios
    res['M2'] = 'EXACT' if max(r['rP'] for r in rows) <= 1e-10 else (
        'SECOND' if all(min(x['r1'], x['r2']) >= 3.5 for x in ratios) else 'FAIL')
    res['cpu_s'] = time.process_time() - t0
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    (OUT / 'matter.json').write_text(json.dumps(res, **DUMP))
    print('M1', res['M1'], 'M2', res['M2'], [(x['dir'], x['t'], round(x['r1'], 2), round(x['r2'], 2)) for x in ratios])


# ---------------- B (high precision roots) ----------------
def mp_sync_L(q, dps=50):
    """Synchronous-gauge companion matrix built in mpmath arithmetic (the zero cluster is defective, so float64
    eigenvalues split by eps^(1/p); exact roots need the matrix itself in high precision)."""
    import mpmath as mp
    mp.mp.dps = dps
    Hm = mp.mpf(H) if False else mp.mpf('0.2')
    K = [2 * mp.sin(mp.mpf(x) / 2) / Hm for x in q]
    eta = mp.diag([-1, 1, 1, 1])

    def Gup(hm, w):
        kl = mp.matrix([[-w], [K[0]], [K[1]], [K[2]]])
        ku = eta * kl; k2 = (kl.T * ku)[0]; tr = sum(eta[i, i] * hm[i, i] for i in range(4))
        kh = hm * ku; khk = (ku.T * hm * ku)[0]
        G = mp.matrix(4, 4)
        for i in range(4):
            for j in range(4):
                G[i, j] = 0.5 * (k2 * hm[i, j] - kl[i] * kh[j] - kh[i] * kl[j] + kl[i] * kl[j] * tr
                                 + eta[i, j] * khk - eta[i, j] * k2 * tr)
        return eta * G * eta

    def kern(w):
        M = mp.matrix(NX, NX)
        for j in range(NX):
            hm = mp.matrix(4, 4)
            for r in range(16):
                hm[r // 4, r % 4] = PM[r, j]
            Gm = Gup(hm, w)
            for i in range(NX):
                M[i, j] = sum(PM[r, i] * Gm[r // 4, r % 4] for r in range(16))
        return M
    M0 = kern(mp.mpf(0)); Mp = kern(mp.mpf(1)); Mn = kern(mp.mpf(-1))
    A = (Mp + Mn) / 2 - M0
    Ahh = mp.matrix([[A[i, j] for j in IH] for i in IH]); Chh = mp.matrix([[M0[i, j] for j in IH] for i in IH])
    F = mp.lu_solve(Ahh, Chh) if False else mp.inverse(Ahh) * Chh
    L = mp.zeros(12, 12)
    for i in range(6):
        L[i, 6 + i] = 1
        for j in range(6):
            L[6 + i, j] = F[i, j]
    Kn = mp.sqrt(sum(x ** 2 for x in K))
    return L, Kn


def spec_mp():
    """Roots of the synchronous-gauge system in 50-digit arithmetic for every spec case (G1, G2 counting)."""
    import mpmath as mp
    t0 = time.process_time(); r = json.loads((OUT / 'spec.json').read_text())
    reg = json.loads((ROOT / 'results/EX-GC34F/spectrum.json').read_text())
    dirs = dict(reg['directions']); dirs.update(reg['approach'])
    rows = []
    for c in r['cases']:
        d = np.array(dirs.get(c['dir'], {'100': (1, 0, 0), 'd01+': (1, 1, 0), 'b1++': (1, 1, 1)}.get(c['dir'])), float)
        d = d / np.linalg.norm(d); q = [mp.mpf(repr(float(x))) for x in c['kh'] * d]
        L, Kn = mp_sync_L(q)
        ev = mp.eig(L, left=False, right=False)
        lam = [complex(x) for x in ev]; K = float(Kn)
        grow = max(abs(l.real) for l in lam) / K
        nwave = sum(abs(abs(l) / K - 1) < 1e-10 for l in lam)
        nzero = sum(abs(l) <= 1e-8 * K for l in lam)
        rows.append(dict(dir=c['dir'], kh=c['kh'], max_re_over_K=grow, n_wave=nwave, n_zero=nzero,
                         max_zero=max([abs(l) / K for l in lam if abs(l) <= 1e-3 * K] or [0.0])))
    res = dict(rows=rows, G1=bool(max(x['max_re_over_K'] for x in rows) <= 1e-8),
               G2_counts=bool(all(x['n_wave'] == 4 and x['n_zero'] == 8 for x in rows)), cpu_s=time.process_time() - t0,
               note='50-digit roots; float64 eig splits the defective 8-fold zero cluster by eps^(1/4)')
    (OUT / 'spec_mp.json').write_text(json.dumps(res, **DUMP))
    print('mp G1', res['G1'], max(x['max_re_over_K'] for x in rows), 'counts', res['G2_counts'],
          set((x['n_wave'], x['n_zero']) for x in rows), 'max zero', max(x['max_zero'] for x in rows), 'cpu', res['cpu_s'])


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    dict(build=build, spec=spec, spec_mp=spec_mp, static=static, sync=sync, matter=matter)[sys.argv[1]]()
