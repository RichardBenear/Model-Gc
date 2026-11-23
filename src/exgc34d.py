"""EX-GC34D (exploratory): lapse x Regge gravity energy on the cubic lattice -> static weak field and PPN gamma.

Kuhn triangulation (6 tetrahedra per cube, 7 edge types per vertex). Variables per vertex: lapse n and 7 edge-length
changes dl. Gravity energy S_g = (1/kappa) sum_tet sum_{e in tet} N_e l_e (theta0_e - theta_e(l)), N_e = mean of end
lapses; matter -sum_v N_v M_v. Linearised by high-precision finite differences per tetrahedron, assembled into the
8x8 Fourier kernel K(k); vertex-displacement null modes fixed with K + G G^dagger.

  python src/exgc34d.py          -> results/EX-GC34D/analysis.json (D1, D2, D3) + local_raw/EX-GC34D/gauss_fields.npz
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import itertools, json, platform, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC34D'; RAW = ROOT / 'local_raw/EX-GC34D'
sys.path.insert(0, str(ROOT / 'src'))

H = 0.2; KAPPA = 1.0
ET = np.array([(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1), (1, 0, 1), (1, 1, 1)])
NBOX = (512, 256, 256)
MASS = 0.4 * np.pi            # far field n -> -kappa M / (8 pi r) = -0.05 / r
SIG_SRC = 1.0

# ---------------- triangulation ----------------
TETS = []
for s in itertools.permutations(range(3)):
    o = [np.zeros(3, int)]
    for a in s:
        o.append(o[-1] + np.eye(3, dtype=int)[a])
    TETS.append(np.array(o))


def tet_edges(o):
    """[(i, j, type)] for the 6 edges of a tetrahedron with vertex offsets o (4x3)."""
    out = []
    for i, j in itertools.combinations(range(4), 2):
        d = o[j] - o[i]
        out.append((i, j, int([t for t in range(7) if (ET[t] == d).all()][0])))
    return out


# ---------------- exact geometry (mpmath) ----------------
def mp_coords(L2):
    """Embed a tetrahedron from squared edge lengths L2[i][j] (vertex 0 at the origin)."""
    G = [[(L2[0][a + 1] + L2[0][b + 1] - L2[a + 1][b + 1]) / 2 for b in range(3)] for a in range(3)]
    C = mp.cholesky(mp.matrix(G))
    return [mp.matrix([0, 0, 0])] + [mp.matrix([C[a, 0], C[a, 1], C[a, 2]]) for a in range(3)]


def mp_dihedral(P, i, j):
    k, l = [m for m in range(4) if m not in (i, j)]
    e = P[j] - P[i]; e = e / mp.norm(e)
    u = P[k] - P[i]; v = P[l] - P[i]
    u = u - e * (u.T * e)[0]; v = v - e * (v.T * e)[0]
    return mp.acos((u.T * v)[0] / (mp.norm(u) * mp.norm(v)))


def mp_theta(o, l):
    es = tet_edges(o); L2 = [[mp.mpf(0)] * 4 for _ in range(4)]
    for (i, j, _), le in zip(es, l):
        L2[i][j] = L2[j][i] = le ** 2
    P = mp_coords(L2)
    return [mp_dihedral(P, i, j) for i, j, _ in es]


def l0_of(o):
    return [mp.mpf(H) * mp.sqrt(int((ET[t] ** 2).sum())) for _, _, t in tet_edges(o)]


def S_tet(o, x, th0):
    """x: 4 lapse perturbations + 6 length changes (mp)."""
    es = tet_edges(o); l0 = l0_of(o)
    l = [l0[e] + x[4 + e] for e in range(6)]
    th = mp_theta(o, l)
    return sum((1 + (x[i] + x[j]) / 2) * l[e] * (th0[e] - th[e]) for e, (i, j, _) in enumerate(es))


def tet_hessian(o, d=mp.mpf('1e-12')):
    th0 = mp_theta(o, l0_of(o))
    f = lambda x: S_tet(o, x, th0)
    Hm = np.zeros((10, 10)); I = [[mp.mpf(0)] * 10 for _ in range(10)]
    for a in range(10):
        I[a][a] = d
    add = lambda p, q, sp, sq: [sp * p[m] + sq * q[m] for m in range(10)]
    for a in range(10):
        for b in range(a, 10):
            v = (f(add(I[a], I[b], 1, 1)) - f(add(I[a], I[b], 1, -1)) - f(add(I[a], I[b], -1, 1))
                 + f(add(I[a], I[b], -1, -1))) / (4 * d * d)
            Hm[a, b] = Hm[b, a] = float(v)
    return Hm


def schlafli(o, d=mp.mpf('1e-15')):
    """max_f |sum_e l_e dtheta_e/dl_f| (relative to max |dtheta/dl| * l)."""
    l0 = l0_of(o); worst = 0.0; scale = 0.0
    for f in range(6):
        lp = list(l0); lm = list(l0); lp[f] += d; lm[f] -= d
        tp = mp_theta(o, lp); tm = mp_theta(o, lm)
        dth = [(tp[e] - tm[e]) / (2 * d) for e in range(6)]
        worst = max(worst, abs(float(sum(l0[e] * dth[e] for e in range(6)))))
        scale = max(scale, max(abs(float(l0[e] * dth[e])) for e in range(6)))
    return worst / scale


# ---------------- Fourier kernel ----------------
def build_terms():
    """Aggregate per-tet Hessians into {(ca, cb, offset): coefficient}. Component 0 = lapse, 1..7 = edge types."""
    terms = {}; hess = []
    for o in TETS:
        Hm = tet_hessian(o); hess.append(Hm)
        var = [(tuple(o[i]), 0) for i in range(4)] + [(tuple(o[i]), 1 + t) for i, j, t in tet_edges(o)]
        for a, (oa, ca) in enumerate(var):
            for b, (ob, cb) in enumerate(var):
                key = (ca, cb, tuple(np.array(ob) - np.array(oa)))
                terms[key] = terms.get(key, 0.0) + Hm[a, b]
    return terms, hess


def kernel(terms, q):
    """q (..., 3) lattice wavenumbers (k h) -> K (..., 8, 8)."""
    C = {}
    for (ca, cb, off), c in terms.items():
        C.setdefault(off, np.zeros((8, 8)))[ca, cb] += c
    K = np.zeros(q.shape[:-1] + (8, 8), complex)
    for off, Cm in C.items():
        K += np.exp(1j * (q @ np.array(off, float)))[..., None, None] * Cm
    return K


def gauge_modes(q):
    G = np.zeros(q.shape[:-1] + (8, 3), complex)
    for t in range(7):
        d = ET[t].astype(float)
        G[..., 1 + t, :] = (d / np.linalg.norm(d)) * (np.exp(1j * (q @ d)) - 1)[..., None]
    return G


def response(terms, q):
    """Unit-mass response x = kappa (K + G G^+)^{-1} e_n for each q; q = 0 entries return 0."""
    K = kernel(terms, q); G = gauge_modes(q)
    Kg = K + G @ np.conj(np.swapaxes(G, -1, -2))
    z = np.all(q == 0, axis=-1)
    Kg[z] = np.eye(8)
    f = np.zeros(q.shape[:-1] + (8,), complex); f[..., 0] = KAPPA
    x = np.linalg.solve(Kg, f[..., None])[..., 0]
    x[z] = 0
    return x


def metric_from_edges(x, q):
    """x (...,8) -> n, h (...,3,3) at the vertex (edges moved from midpoints by phase e^{-i q.d/2})."""
    e = [x[..., 1 + t] * np.exp(-0.5j * (q @ ET[t].astype(float))) for t in range(7)]
    hm = np.zeros(x.shape[:-1] + (3, 3), complex)
    for i in range(3):
        hm[..., i, i] = 2 * e[i] / H                         # d^T h d = 2 |d| dl / h, |d| = 1
    for (i, j), t in (((0, 1), 3), ((1, 2), 4), ((0, 2), 5)):
        hm[..., i, j] = hm[..., j, i] = (2 * np.sqrt(2) * e[t] / H - hm[..., i, i] - hm[..., j, j]) / 2
    d = np.ones(3); body = 2 * np.sqrt(3) * e[6] / H        # d^T h d from the body diagonal
    resid = body - np.einsum('i,...ij,j->...', d, hm, d)
    return x[..., 0], hm, resid


def gamma_of(n, hm, q):
    kh = q / np.linalg.norm(q, axis=-1, keepdims=True)
    P = np.eye(3) - kh[..., :, None] * kh[..., None, :]
    A = 0.5 * np.einsum('...ij,...ij->...', P, hm)
    return -A / (2 * n)


# ---------------- D1 separate path: direct Regge sum ----------------
def np_dihedrals(P):
    """P (..., 4, 3) -> angles (..., 6) in tet_edges order."""
    out = []
    for i, j in itertools.combinations(range(4), 2):
        k, l = [m for m in range(4) if m not in (i, j)]
        e = P[..., j, :] - P[..., i, :]; e = e / np.linalg.norm(e, axis=-1, keepdims=True)
        u = P[..., k, :] - P[..., i, :]; v = P[..., l, :] - P[..., i, :]
        u = u - e * (u * e).sum(-1, keepdims=True); v = v - e * (v * e).sum(-1, keepdims=True)
        out.append(np.arccos(np.clip((u * v).sum(-1) / np.linalg.norm(u, axis=-1) / np.linalg.norm(v, axis=-1),
                                     -1, 1)))
    return np.stack(out, -1)


def np_embed(L2):
    """L2 (..., 4, 4) squared lengths -> coordinates (..., 4, 3)."""
    G = 0.5 * (L2[..., 0, 1:, None] + L2[..., 0, None, 1:] - L2[..., 1:, 1:])
    C = np.linalg.cholesky(G)
    return np.concatenate([np.zeros(C.shape[:-2] + (1, 3)), C], -2)


def direct_sum(nf, dl):
    """nf (N,N,N) lapse perturbation, dl (7,N,N,N) length changes: Regge sum with N_e l_e (theta0 - theta)."""
    total = 0.0
    for o in TETS:
        es = tet_edges(o)
        sh = lambda f, off: np.roll(f, shift=tuple(-np.array(off)), axis=(-3, -2, -1))
        l = np.stack([H * np.linalg.norm(ET[t]) + sh(dl[t], o[i]) for i, j, t in es], -1)
        L2 = np.zeros(l.shape[:-1] + (4, 4))
        for e, (i, j, t) in enumerate(es):
            L2[..., i, j] = L2[..., j, i] = l[..., e] ** 2
        th = np_dihedrals(np_embed(L2))
        th0 = np_dihedrals(o.astype(float) * H)
        Nv = [1 + sh(nf, o[i]) for i in range(4)]
        for e, (i, j, t) in enumerate(es):
            total += ((Nv[i] + Nv[j]) / 2 * l[..., e] * (th0[e] - th[..., e])).sum()
    return total


def d1_direct(terms, N=8, eps=1e-3, seed=341):
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((8, N, N, N)) * np.array([1.0] + [H] * 7)[:, None, None, None]
    S = lambda s: direct_sum(s * x[0], s * x[1:])
    Qd = (S(eps) + S(-eps) - 2 * S(0.0)) / (2 * eps ** 2)
    xh = np.fft.fftn(x, axes=(1, 2, 3)); kk = 2 * np.pi * np.fft.fftfreq(N)
    q = np.stack(np.meshgrid(kk, kk, kk, indexing='ij'), -1)
    K = kernel(terms, q)
    Qf = 0.5 * np.einsum('axyz,xyzab,bxyz->', np.conj(xh), K, xh).real / N ** 3
    return dict(N=N, eps=eps, Q_direct=float(Qd), Q_fourier=float(Qf), rel=float(abs(Qd / Qf - 1)), S0=float(S(0.0)))


# ---------------- sources ----------------
def box_axes():
    return [(np.arange(NBOX[a]) - NBOX[a] // 2) * H for a in range(3)]


def gauss_source():
    ax = box_axes()
    r2 = ax[0][:, None, None] ** 2 + ax[1][None, :, None] ** 2 + ax[2][None, None, :] ** 2
    m = np.exp(-r2 / (2 * SIG_SRC ** 2)).astype(np.float64)
    return m * (MASS / m.sum())


def knot_density(B):
    """Per-site energy density of the 34B rest knot (h = 0.2), orientation stencil as in exgc34b."""
    import exgc34b as kb
    import exgc26a_skyrmion as ref
    e1 = kb.e1mod(0.2); G = ref.Grid(0.2, 7.0)
    p = np.load(e1.restpath(B)); parts = kb.knot_parts(G, p)
    e = None
    for m, D in parts:
        tr = np.einsum('ii...->...', D); D2 = np.einsum('ik...,kj...->ij...', D, D)
        ed = (tr + 0.5 * (tr * tr - np.einsum('ii...->...', D2))) * m / 8
        e = ed if e is None else e + ed
    return e


def embed_center(e):
    """Place a small array at the box centre, scaled to MASS."""
    out = np.zeros(NBOX); c = [NBOX[a] // 2 - e.shape[a] // 2 for a in range(3)]
    out[c[0]:c[0] + e.shape[0], c[1]:c[1] + e.shape[1], c[2]:c[2] + e.shape[2]] = e
    return out * (MASS / out.sum())


# ---------------- D2 ----------------
def d2(terms):
    res = {}
    for name, dv in (('100', (1, 0, 0)), ('110', (1, 1, 0)), ('111', (1, 1, 1))):
        dv = np.array(dv, float) / np.linalg.norm(dv); rows = []
        for qm in (0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8):
            q = (qm * dv)[None]
            x = response(terms, q); n, hm, resid = metric_from_edges(x, q)
            g = gamma_of(n, hm, q)[0]
            newton = (qm / H) ** 2 * H ** 3 * n[0] / KAPPA            # k^2 h^3 n / (kappa M) -> -1/2
            Kg = kernel(terms, q)[0] + (gauge_modes(q)[0] @ np.conj(gauge_modes(q)[0].T))
            ev = np.sort(np.abs(np.linalg.eigvalsh(Kg)))
            rows.append(dict(kh=qm, gamma=float(g.real), gamma_imag=float(abs(g.imag)),
                             newton=float(newton.real), body_resid_over_n=float(abs(resid[0] / n[0])),
                             abs_eigs=ev.tolist()))
        res[name] = rows
    return res


# ---------------- D3 / fields ----------------
def straight_ray(xn, xh, mhat_plane, bs):
    """kx = 0 plane: deflection of a ray along x at (y=b, z=0) for Phi = -n (lapse) and -n + h_xx/2 (full)."""
    out = {}
    for lab, phi in (('lapse', -xn), ('full', -xn + 0.5 * xh)):
        P = phi * mhat_plane
        ky = 2 * np.pi * np.fft.fftfreq(NBOX[1], d=H)[:, None]
        a = np.fft.ifft2(1j * ky * P).real               # d/dy of the x-projection (up to a common factor)
        out[lab] = [float(a[int(round(b / H)), 0]) for b in bs]
    out['ratio'] = [f / l for f, l in zip(out['full'], out['lapse'])]
    return out


def solve_plane(terms):
    """Response on the kx = 0 plane (full kernel, not separable)."""
    ky = 2 * np.pi * np.fft.fftfreq(NBOX[1]); kz = 2 * np.pi * np.fft.fftfreq(NBOX[2])
    q = np.stack(np.meshgrid([0.0], ky, kz, indexing='ij'), -1)[0]
    x = response(terms, q); n, hm, _ = metric_from_edges(x, q)
    return n, hm[..., 0, 0]


def mhat_plane(src):
    """FFT of the source at kx = 0, with the box centre moved to the origin."""
    s = np.fft.ifftshift(src)
    return np.fft.fft2(s.sum(0))


def solve_full(terms, src, path):
    """Full 3D response for one source: stores n_hat and h_hat (6 comps) on the rfft grid (origin-centred)."""
    sh = np.fft.rfftn(np.fft.ifftshift(src))
    kx = 2 * np.pi * np.fft.fftfreq(NBOX[0]); ky = 2 * np.pi * np.fft.fftfreq(NBOX[1])
    kz = 2 * np.pi * np.fft.rfftfreq(NBOX[2])
    nh = np.zeros(sh.shape, complex); hh = np.zeros((6,) + sh.shape, complex)
    QY, QZ = np.meshgrid(ky, kz, indexing='ij')
    for i in range(NBOX[0]):
        q = np.stack([np.full(QY.shape, kx[i]), QY, QZ], -1)
        x = response(terms, q); n, hm, _ = metric_from_edges(x, q)
        nh[i] = n * sh[i]
        for c, (a, b) in enumerate(((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))):
            hh[c, i] = hm[..., a, b] * sh[i]
    RAW.mkdir(parents=True, exist_ok=True)
    np.savez(path, nh=nh, hh=hh)
    return nh, hh


def platform_info():
    import scipy
    return {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
            'mpmath': mp.__version__, 'machine': platform.platform()}


def main():
    mp.mp.dps = 50
    t0 = time.process_time(); res = {}
    terms, hess = build_terms()
    # D1
    res['schlafli_rel'] = max(schlafli(o) for o in TETS)
    res['tet_nn_block_max'] = float(max(abs(Hm[:4, :4]).max() for Hm in hess))
    qs = np.array([[0.05, 0, 0], [0.03, 0.02, 0.01], [np.pi, 0, 0], [np.pi, np.pi, np.pi], [1.0, -2.0, 0.5]])
    K = kernel(terms, qs); G = gauge_modes(qs)
    res['K_nn_max'] = float(abs(K[..., 0, 0]).max())
    res['K_herm_max'] = float(abs(K - np.conj(np.swapaxes(K, -1, -2))).max())
    res['KG_rel'] = [float(np.linalg.norm(K[i] @ G[i]) / np.linalg.norm(K[i])) for i in range(len(qs))]
    res['D1_direct'] = d1_direct(terms)
    print('D1', {k: res[k] for k in ('schlafli_rel', 'K_nn_max', 'KG_rel')}, res['D1_direct'], flush=True)
    # D2
    res['D2'] = d2(terms)
    for k, rows in res['D2'].items():
        print('D2', k, [(r['kh'], round(r['gamma'], 6), round(r['newton'], 6)) for r in rows], flush=True)
    # D3
    bs = [6.0, 8.0, 12.0]
    xn, xh = solve_plane(terms)
    srcs = {'gauss': gauss_source(), 'knot_B1': embed_center(knot_density(1)), 'knot_B2': embed_center(knot_density(2))}
    res['D3'] = {k: straight_ray(xn, xh, mhat_plane(s), bs) for k, s in srcs.items()}
    res['D3_b'] = bs
    for k, v in res['D3'].items():
        print('D3', k, v['ratio'], flush=True)
    # full 3D Gaussian field for D4
    t1 = time.time()
    solve_full(terms, srcs['gauss'], RAW / 'gauss_fields.npz')
    res['full_solve_wall_s'] = time.time() - t1
    res['cpu_s'] = time.process_time() - t0
    res['platform'] = platform_info()
    res['params'] = dict(h=H, kappa=KAPPA, box=NBOX, mass=MASS, sigma_src=SIG_SRC)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'analysis.json').write_text(json.dumps(res, indent=1) + '\n')
    print('done cpu', res['cpu_s'])


if __name__ == '__main__':
    main()
