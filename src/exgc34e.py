"""EX-GC34E (exploratory): lapse x Regge (34D) + shift link variable with lattice DeWitt kinetic energy.

Per vertex x = (n, beta_x, beta_y, beta_z, dl_1..dl_7) (11). Fourier e^{i q.v - i w t}. Edge strain rate
s_e = l0 (-i w dl_e - Dhat.(beta_end - beta_start)) = D^T K D fixes K_ij per tetrahedron; kinetic
(1/2kappa) sum V (K:K - trK^2). Potential: 34D Regge Hessian. Matter: -N M + beta.P.
L2 = 1/2 x^+ M(w,q) x / kappa - Re(n* M) + Re(beta*.P)  ->  M x = kappa (M e_n - P e_beta).

  python src/exgc34e.py      -> results/EX-GC34E/analysis.json
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

OUT = ROOT / 'results/EX-GC34E'
H = D.H; KAPPA = D.KAPPA; ET = D.ET
NV = 11
SYM = ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))
GHAT = np.diag([1.0, 1, 1, 2, 2, 2]) - np.outer([1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0])   # K:K - (trK)^2
VTET = H ** 3 / 6


def tet_static():
    """34D static kernel terms embedded in 11 components."""
    terms, _ = D.build_terms()
    out = {}
    for (ca, cb, off), c in terms.items():
        a = 0 if ca == 0 else 3 + ca; b = 0 if cb == 0 else 3 + cb
        out[(a, b, off)] = out.get((a, b, off), 0.0) + c
    return out


def static_kernel(terms, q):
    C = {}
    for (a, b, off), c in terms.items():
        C.setdefault(off, np.zeros((NV, NV)))[a, b] += c
    K = np.zeros(q.shape[:-1] + (NV, NV), complex)
    for off, Cm in C.items():
        K += np.exp(1j * (q @ np.array(off, float)))[..., None, None] * Cm
    return K


def strain_map(o, q, w):
    """A (..., 6, 11): K_ij (SYM order) of tetrahedron with vertex offsets o, from the Fourier amplitudes at its base
    vertex (phases e^{i q.o} applied)."""
    es = D.tet_edges(o)
    Wm = np.zeros((6, 6))
    S = np.zeros(q.shape[:-1] + (6, NV), complex)
    for e, (i, j, t) in enumerate(es):
        Dv = ET[t] * H; l0 = np.linalg.norm(Dv); Dh = Dv / l0
        for c, (a, b) in enumerate(SYM):
            Wm[e, c] = Dv[a] * Dv[b] * (1 if a == b else 2)
        pi = np.exp(1j * (q @ o[i].astype(float))); pj = np.exp(1j * (q @ o[j].astype(float)))
        S[..., e, 4 + t] += l0 * (-1j * w) * pi
        for a in range(3):
            S[..., e, 1 + a] += -l0 * Dh[a] * (pj - pi)
    return np.einsum('cd,...dv->...cv', np.linalg.inv(Wm), S)


def kinetic_kernel(q, w):
    T = np.zeros(q.shape[:-1] + (NV, NV), complex)
    for o in D.TETS:
        A = strain_map(o, q, w)
        T += VTET * np.einsum('...cv,cd,...du->...vu', np.conj(A), GHAT, A)
    return T


def full_kernel(st, q, w):
    return kinetic_kernel(q, w) + static_kernel(st, q)


def gauge_modes(q, w):
    """Vertex displacement u: beta += du/dt = -i w u, dl_e += Dhat.(u(v+d) - u(v))."""
    G = np.zeros(q.shape[:-1] + (NV, 3), complex)
    for a in range(3):
        G[..., 1 + a, a] = -1j * w
    for t in range(7):
        d = ET[t].astype(float)
        G[..., 4 + t, :] = (d / np.linalg.norm(d)) * (np.exp(1j * (q @ d)) - 1)[..., None]
    return G


def metric_from_x(x, q):
    """Vertex metric: h00 = -2n, h0i = beta_i, h_ij from edges (34D map)."""
    x8 = np.concatenate([x[..., :1], x[..., 4:]], -1)
    n, hm, _ = D.metric_from_edges(x8, q)
    h = np.zeros(x.shape[:-1] + (4, 4), complex)
    h[..., 0, 0] = -2 * n
    for a in range(3):
        h[..., 0, 1 + a] = h[..., 1 + a, 0] = x[..., 1 + a]
    h[..., 1:, 1:] = hm
    return h


def riemann(h, kl):
    """Linearised Riemann R_abcd, a<b, c<d (36 comps); kl lower 4-wavevector (-w, k)."""
    out = []
    for a in range(4):
        for b in range(a + 1, 4):
            for c in range(4):
                for d in range(c + 1, 4):
                    out.append(0.5 * (kl[b] * kl[c] * h[a, d] + kl[a] * kl[d] * h[b, c]
                                      - kl[a] * kl[c] * h[b, d] - kl[b] * kl[d] * h[a, c]))
    return np.array(out)


# ---------------- continuum linearised GR (independent solve) ----------------
ETA = np.diag([-1.0, 1, 1, 1])


def gr_solution(k, v, rho=1.0):
    """Solve G_{mu nu}(h) = kappa T_{mu nu} for dust (T_00 = rho, T_0i = -rho v_i, T_ij = rho v_i v_j)."""
    w = k @ v
    kl = np.array([-w, *k]); ku = ETA @ kl; k2 = kl @ ku
    idx = [(m, n) for m in range(4) for n in range(m, 4)]
    def G_of(hm):
        tr = np.trace(ETA @ hm); kh = hm @ ku; khk = ku @ hm @ ku
        return 0.5 * (k2 * hm - np.outer(kl, kh) - np.outer(kh, kl) + np.outer(kl, kl) * tr + ETA * khk - ETA * k2 * tr)
    Amat = np.zeros((10, 10), complex)
    for c, (m, n) in enumerate(idx):
        hm = np.zeros((4, 4), complex); hm[m, n] = hm[n, m] = 1
        Gm = G_of(hm)
        Amat[:, c] = [Gm[a, b] for a, b in idx]
    T = np.zeros((4, 4)); T[0, 0] = rho
    T[0, 1:] = T[1:, 0] = -rho * v; T[1:, 1:] = rho * np.outer(v, v)
    sol, res, rk, sv = np.linalg.lstsq(Amat, KAPPA * np.array([T[a, b] for a, b in idx]), rcond=1e-12)
    hm = np.zeros((4, 4), complex)
    for c, (m, n) in enumerate(idx):
        hm[m, n] = hm[n, m] = sol[c]
    resid = np.linalg.norm(Amat @ sol - KAPPA * np.array([T[a, b] for a, b in idx]))
    return hm, kl, float(resid)


# ---------------- E4 forced solve ----------------
def lattice_solution(st, q, v_lat, M=1.0):
    """x at w = k.v with k = q/h; source M e_n - P e_beta, P = M v. Gauge penalty + lstsq."""
    k = q / H; w = k @ v_lat
    Mm = full_kernel(st, q[None], w)[0]; G = gauge_modes(q[None], w)[0]
    Mg = Mm + G @ G.conj().T
    f = np.zeros(NV, complex); f[0] = KAPPA * M; f[1:4] = -KAPPA * M * v_lat
    x, *_ = np.linalg.lstsq(Mg, f, rcond=1e-14)
    resid = np.linalg.norm(Mm @ x - f) / np.linalg.norm(f)
    mom_resid = np.linalg.norm((Mm @ x - f)[1:4]) / np.linalg.norm(f)
    sv = np.linalg.svd(Mg, compute_uv=False)
    return x, w, float(resid), float(mom_resid), float(sv[-1] / sv[0])


def e4_case(st, qdir, kh, vdir, vmag=1e-3):
    qd = np.array(qdir, float) / np.linalg.norm(qdir); q = kh * qd; k = q / H
    vd = np.array(vdir, float); vd = vd / np.linalg.norm(vd); v = vmag * vd
    R = {}
    for sgn in (1, -1):
        x, w, res, mres, cond = lattice_solution(st, q, sgn * v)
        h = metric_from_x(x, q); kl = np.array([-w, *k])
        Rl = riemann(h, kl) * (H ** 3)            # lattice M is per vertex: continuum rho_hat = M / h^3
        hg, klg, gres = gr_solution(k, sgn * v)
        Rg = riemann(hg, klg)
        R[sgn] = (Rl, Rg, res, mres, cond, gres)
    R0l = (R[1][0] + R[-1][0]) / 2; R0g = (R[1][1] + R[-1][1]) / 2
    R1l = (R[1][0] - R[-1][0]) / 2; R1g = (R[1][1] - R[-1][1]) / 2
    return dict(kh=kh, qdir=list(qdir), vdir=list(vdir),
                err_static=float(np.abs(R0l - R0g).max() / np.abs(R0g).max()),
                err_Ov=float(np.abs(R1l - R1g).max() / np.abs(R1g).max()),
                maxR_Ov=float(np.abs(R1g).max()),
                residual=max(R[1][2], R[-1][2]), momentum_residual=max(R[1][3], R[-1][3]),
                cond=R[1][4], gr_residual=max(R[1][5], R[-1][5]))


# ---------------- E3 waves ----------------
def wave_modes(st, q):
    """Finite eigenfrequencies of the vacuum pencil after dropping axis-edge columns (gauge) and axis-edge rows
    (implied by the gauge identity). M(w) = C + w B + w^2 A; quadratic eigenproblem linearised."""
    ws = (0.0, 1.0, -1.0)
    Ms = [full_kernel(st, q[None], w)[0] for w in ws]
    C = Ms[0]; A = (Ms[1] + Ms[2] - 2 * C) / 2; B = (Ms[1] - Ms[2]) / 2
    keep = [0, 1, 2, 3] + [4 + t for t in range(3, 7)]        # drop axis edges (types 0..2)
    C, B, A = [X[np.ix_(keep, keep)] for X in (C, B, A)]
    n = len(keep); I = np.eye(n); Z = np.zeros((n, n))
    L1 = np.block([[Z, I], [-C, -B]]); L2 = np.block([[I, Z], [Z, A]])
    ev = sl.eig(L1, L2, right=False)
    ev = ev[np.isfinite(ev)]
    return ev


def light_speed(q):
    """34A massless field on a flat background, isotropic split (axis 1/3, diagonals 1/6), continuous time."""
    dirs = np.array([(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)])
    wts = np.array([1 / 3] * 3 + [1 / 6] * 6)
    om2 = (wts * 4 * np.sin(dirs @ q / 2) ** 2).sum() / H ** 2
    return np.sqrt(om2) / (np.linalg.norm(q) / H)


def main():
    mp.mp.dps = 50
    t0 = time.process_time(); res = {}
    st = tet_static()
    # E1
    rng = np.random.default_rng(345)
    samp = [(rng.uniform(-np.pi, np.pi, 3), rng.uniform(-3, 3)) for _ in range(6)] + [(np.array([np.pi] * 3), 2.0)]
    kg = []
    for q, w in samp:
        Mm = full_kernel(st, q[None], w)[0]; G = gauge_modes(q[None], w)[0]
        kg.append(float(np.linalg.norm(Mm @ G) / np.linalg.norm(Mm)))
    res['E1_MG_rel_max'] = max(kg)
    q = np.array([0.3, -0.7, 1.1])
    M0 = full_kernel(st, q[None], 0.0)[0]
    K34 = D.kernel(D.build_terms()[0], q[None])[0]
    sel = [0] + list(range(4, 11))
    res['E1_static_vs_34D'] = float(np.abs(M0[np.ix_(sel, sel)] - K34).max() / np.abs(K34).max())
    res['E1_herm_max'] = float(max(np.abs(full_kernel(st, q[None], w)[0] - full_kernel(st, q[None], w)[0].conj().T).max()
                                   for w in (0.7, -1.3)))
    print('E1', res['E1_MG_rel_max'], res['E1_static_vs_34D'], res['E1_herm_max'], flush=True)
    # E2: refoliation. Continuum: beta -> beta + grad xi0, n -> n - i w xi0. Lattice gradient candidates at vertex:
    # symmetric difference sum over axes. Report how far the best (n, beta-longitudinal) combination is from null.
    e2 = []
    for kh in (0.05, 0.1, 0.2, 0.4):
        for dn, dv in (('100', (1, 0, 0)), ('111', (1, 1, 1))):
            qq = kh * np.array(dv, float) / np.linalg.norm(dv)
            for w in (0.0, 0.5 * kh / H):
                Mm = full_kernel(st, qq[None], w)[0]; G = gauge_modes(qq[None], w)[0]
                # complement of the displacement gauge, then smallest singular value
                Qc = sl.null_space(G.conj().T)
                sv = np.linalg.svd(Qc.conj().T @ Mm @ Qc, compute_uv=False)
                bb = np.linalg.svd(Mm[1:4, 1:4], compute_uv=False)
                e2.append(dict(kh=kh, dir=dn, w=w, smin_over_smax=float(sv[-1] / sv[0]),
                               sv_tail=[float(s / sv[0]) for s in sv[-3:]],
                               beta_block_smin_over_smax=float(bb[-1] / bb[0])))
    res['E2'] = e2
    for r in e2:
        print('E2', r['kh'], r['dir'], round(r['w'], 3), '%.2e' % r['smin_over_smax'], ['%.1e' % s for s in r['sv_tail']],
              '%.2e' % r['beta_block_smin_over_smax'], flush=True)
    # E3
    e3 = []
    for kh in (0.05, 0.1, 0.2):
        for dn, dv in (('100', (1, 0, 0)), ('110', (1, 1, 0)), ('111', (1, 1, 1))):
            qq = kh * np.array(dv, float) / np.linalg.norm(dv)
            ev = wave_modes(st, qq)
            real = sorted(float(e.real) for e in ev if abs(e.imag) < 1e-8 * max(1, abs(e)) and e.real > 0)
            light = [e for e in real if e * H < 0.5]
            allf = sorted([[float(e.real), float(e.imag)] for e in ev], key=lambda z: abs(complex(*z)))
            c_light = light_speed(qq)
            e3.append(dict(kh=kh, dir=dn, n_light=len(light), speeds=[e / (kh / H) for e in light],
                           c_light=float(c_light), all_finite=allf[:12], n_finite=len(ev)))
    res['E3'] = e3
    for r in e3:
        print('E3', r['kh'], r['dir'], r['n_light'], [round(s, 6) for s in r['speeds']], round(r['c_light'], 6),
              r['n_finite'], flush=True)
    # E4
    e4 = []
    for kh in (0.025, 0.05, 0.1):
        for qdir in ((1, 0, 0), (1, 1, 1), (0.8, 0.5, 0.33)):
            qd = np.array(qdir, float); qd /= np.linalg.norm(qd)
            perp = np.cross(qd, [0, 0, 1.0]) if abs(qd[2]) < 0.9 else np.cross(qd, [1.0, 0, 0])
            for vname, vdir in (('par', qd), ('perp', perp), ('obl', qd + 1.3 * perp)):
                r = e4_case(st, qdir, kh, vdir); r['v'] = vname; e4.append(r)
                print('E4', kh, qdir, vname, '%.2e %.2e res %.1e mom %.1e cond %.1e' % (
                    r['err_static'], r['err_Ov'], r['residual'], r['momentum_residual'], r['cond']), flush=True)
    res['E4'] = e4
    res['cpu_s'] = time.process_time() - t0
    import scipy
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                           mpmath=mp.__version__)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'analysis.json').write_text(json.dumps(res, indent=1) + '\n')


if __name__ == '__main__':
    main()
