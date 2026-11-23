"""EX-GC34F (exploratory): diagnostics of the 34E growing (re-slicing) mode. Plan: spec/EX-GC34F_plan.md.

Uses the committed 34E kernel (src/exgc34e.py) unchanged; diagnostic kernel variants for F3 are built here.

  python src/exgc34f.py f1   -> results/EX-GC34F/spectrum.json
  python src/exgc34f.py f2   -> results/EX-GC34F/evolve.json
  python src/exgc34f.py f3   -> results/EX-GC34F/mechanism.json
Convention e^{-i w t}: +Im w grows. M(w) = C + w B + w^2 A (exactly quadratic); time domain C x + i B x' - A x'' = 0.
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
import exgc34d as D
import exgc34e as E

OUT = ROOT / 'results/EX-GC34F'
H = E.H; ET = E.ET; NV = 11
IN = 0; IB = [1, 2, 3]; IE = list(range(4, 11))
RED = [IN] + IE                                   # zero-shift gauge unknowns / rows
IU = np.triu_indices(4)
SEED = 3434


# ---------------- kernels ----------------
def ghat(lam):
    return np.diag([1.0, 1, 1, 2, 2, 2]) - lam * np.outer([1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0])


def kin_tet(q, w, lam=1.0):
    G = ghat(lam); T = np.zeros((NV, NV), complex)
    for o in D.TETS:
        A = E.strain_map(o, q[None], w)[0]
        T += E.VTET * A.conj().T @ G @ A
    return T


CUBE_EDGES = sorted({(tuple(o[i]), t) for o in D.TETS for (i, j, t) in D.tet_edges(o)})


def kin_cube(q, w, lam=1.0):
    """Per-cube least-squares strain rate from the 19 edges of the cube (diagnostic variant a)."""
    G = ghat(lam)
    Wm = np.zeros((len(CUBE_EDGES), 6)); S = np.zeros((len(CUBE_EDGES), NV), complex)
    for e, (oi, t) in enumerate(CUBE_EDGES):
        oi = np.array(oi, float); oj = oi + ET[t]
        Dv = ET[t] * H; l0 = np.linalg.norm(Dv); Dh = Dv / l0
        for c, (a, b) in enumerate(E.SYM):
            Wm[e, c] = Dv[a] * Dv[b] * (1 if a == b else 2)
        pi = np.exp(1j * (q @ oi)); pj = np.exp(1j * (q @ oj))
        S[e, 4 + t] += l0 * (-1j * w) * pi
        for a in range(3):
            S[e, 1 + a] += -l0 * Dh[a] * (pj - pi)
    A = np.linalg.pinv(Wm) @ S
    return H ** 3 * A.conj().T @ G @ A


class Kernel:
    def __init__(self, st, kin='tet', lam=1.0, lapse_regge=True):
        self.st = dict(st) if lapse_regge else {k: v for k, v in st.items() if k[0] != 0 and k[1] != 0}
        self.kin = kin; self.lam = lam; self.lapse_regge = lapse_regge

    def T(self, q, w):
        return kin_tet(q, w, self.lam) if self.kin == 'tet' else kin_cube(q, w, self.lam)

    def R(self, q):
        return E.static_kernel(self.st, q[None])[0]

    def M(self, q, w):
        return self.T(q, w) + self.R(q)

    def pencil(self, q):
        Ts = [self.T(q, w) for w in (0.0, 1.0, -1.0)]
        T0 = Ts[0]; T1 = (Ts[1] - Ts[2]) / 2; T2 = (Ts[1] + Ts[2] - 2 * T0) / 2
        R = self.R(q)
        return dict(C=T0 + R, B=T1, A=T2, T0=T0, R=R)


def Mw(P, w):
    return P['C'] + w * P['B'] + w * w * P['A']


# ---------------- F1 root finder ----------------
def normal_rank_deficit(P, k, rows=RED, cols=RED, tol=1e-15):
    defs = []
    for g in (0.37, 0.81, 1.3, 2.2):
        s = np.linalg.svd(Mw(P, g * k)[np.ix_(rows, cols)], compute_uv=False)
        defs.append(int((s / s[0] < tol).sum()))
    return max(defs)


def red_roots(P, k, rng, nk=None):
    """Finite roots in the zero-shift gauge; rank-completing perturbation, kept if stable under two draws."""
    C, B, A = [P[x][np.ix_(RED, RED)] for x in ('C', 'B', 'A')]
    n = len(RED); I = np.eye(n); Z = np.zeros((n, n))
    L1 = np.block([[Z, I], [-C, -B]]); L2 = np.block([[I, Z], [Z, A]])
    if nk is None:
        nk = normal_rank_deficit(P, k)
    sc = np.linalg.norm(L1, 2)
    draws = []
    for _ in range(2):
        if nk:
            U = rng.standard_normal((2 * n, nk)) + 1j * rng.standard_normal((2 * n, nk))
            V = rng.standard_normal((2 * n, nk)) + 1j * rng.standard_normal((2 * n, nk))
            D1 = np.diag(rng.standard_normal(nk)); D2 = np.diag(rng.standard_normal(nk))
            ev, vec = sl.eig(L1 + sc * U @ D1 @ V.conj().T, L2 + U @ D2 @ V.conj().T)
        else:
            ev, vec = sl.eig(L1, L2)
        draws.append((ev, vec))
    out = []
    ev0, vec0 = draws[0]; ev1 = draws[1][0]
    for e, v in zip(ev0, vec0.T):
        if not np.isfinite(e) or abs(e) > 1e6 * k:
            continue
        if nk and np.min(np.abs(ev1 - e)) > 1e-8 * max(abs(e), k):
            continue
        x = np.zeros(NV, complex); x[RED] = v[:n]
        nx = np.linalg.norm(x)
        if nx == 0:
            continue
        out.append((complex(e), x / nx))
    return out, nk


def validate(P, k, w, x, nk):
    M = Mw(P, w)
    res = float(np.linalg.norm(M @ x) / (np.linalg.norm(M, 2) * np.linalg.norm(x)))
    r = len(RED) - nk
    def sr(ww):
        s = np.linalg.svd(Mw(P, ww)[np.ix_(RED, RED)], compute_uv=False)
        return s[r - 1]
    d = 1e-3 * abs(w)
    nb = min(sr(w + d), sr(w - d), sr(w + 1j * d), sr(w - 1j * d))
    contrast = float(sr(w) / nb)
    mom = float(np.linalg.norm((M @ x)[IB]) / (np.linalg.norm(M, 2) * np.linalg.norm(x)))
    return res, contrast, mom


def gauge_split(x, q, w):
    h = E.metric_from_x(x, q); kl = np.r_[-w, q / H]
    Gb = []
    for m in range(4):
        hg = np.zeros((4, 4), complex); hg[m, :] += 1j * kl; hg[:, m] += 1j * kl; Gb.append(hg[IU])
    Gb = np.array(Gb).T; hv = h[IU]
    c, *_ = np.linalg.lstsq(Gb, hv, rcond=None)
    ng = hv - Gb @ c
    hn = np.zeros((4, 4), complex); hn[IU] = ng; hn = hn + hn.T - np.diag(np.diag(hn))
    qh = q / np.linalg.norm(q); Pp = np.eye(3) - np.outer(qh, qh)
    hij = h[1:, 1:]
    TT = Pp @ hij @ Pp - 0.5 * Pp * np.trace(Pp @ hij)
    return float(np.linalg.norm(ng) / np.linalg.norm(hv)), float(np.linalg.norm(TT) / max(np.linalg.norm(hn), 1e-300))


def analyse(P, q, rng, want_vectors=False):
    k = np.linalg.norm(q) / H
    roots, nk = red_roots(P, k, rng)
    rows = []
    for w, x in roots:
        if abs(w) <= 1e-8 * k:
            rows.append(dict(cls='zero', w_over_k=[w.real / k, w.imag / k]))
            continue
        res, con, mom = validate(P, k, w, x, nk)
        ngf, ttf = gauge_split(x, q, w)
        rel_im = w.imag / abs(w)
        if abs(rel_im) <= 1e-8:
            cls = 'propagating' if w.real > 0 else 'propagating_neg'
        elif w.imag > 0:
            cls = 'growing'
        else:
            cls = 'decaying'
        r = dict(cls=cls, w_over_k=[w.real / k, w.imag / k], wh=[w.real * H, w.imag * H], resid11=res,
                 contrast=con, mom_resid=mom, non_gauge_frac=ngf, tt_frac_of_nongauge=ttf,
                 abs_n=float(abs(x[IN])), abs_edges=float(np.linalg.norm(x[IE])))
        if want_vectors:
            r['_x'] = x; r['_w'] = w
        rows.append(r)
    return rows, nk


def static_zero_modes(P, q):
    G0 = E.gauge_modes(q[None], 0.0)[0]
    M0 = P['C'] + G0 @ G0.conj().T
    s, V = np.linalg.eigh(M0)
    sc = np.abs(s).max()
    N = V[:, np.abs(s) < 1e-12 * sc]
    out = dict(n_static_zero=int(N.shape[1]))
    if N.shape[1]:
        sb = np.linalg.svd(N.conj().T @ P['B'] @ N, compute_uv=False)
        nb = np.linalg.norm(P['B'], 2)
        out['NBN_sv_rel'] = [float(v / nb) for v in sb]
        out['n_chains'] = int((sb / nb <= 1e-10).sum())
    else:
        out['n_chains'] = 0
    return out


def unit(d):
    d = np.array(d, float); return d / np.linalg.norm(d)


# ---------------- F1 ----------------
SIGNED_PERMS = [(p, s) for p in itertools.permutations(range(3)) for s in itertools.product((1, -1), repeat=3)]


def apply_perm(ps, q):
    p, s = ps
    return np.array([s[i] * q[p[i]] for i in range(3)])


def match_sets(a, b):
    a = np.array(a); b = np.array(b)
    if len(a) != len(b):
        return np.inf
    used = np.zeros(len(b), bool); worst = 0.0
    for z in a:
        d = np.abs(b - z); d[used] = np.inf; j = int(np.argmin(d)); used[j] = True
        worst = max(worst, d[j] / max(abs(z), 1e-300))
    return worst


def f1():
    t0 = time.process_time()
    mp.mp.dps = 50
    st = E.tet_static(); K = Kernel(st); rng = np.random.default_rng(SEED)
    res = dict(convention='e^{-i w t}; +Im w grows', gauge='zero shift (beta = 0)')
    # 1. symmetry test
    rs = np.random.default_rng(SEED + 1)
    sym = []
    for i in range(20):
        q = unit(rs.standard_normal(3)) * rs.uniform(0.05, 0.6)
        base = [w for w, _ in red_roots(K.pencil(q), np.linalg.norm(q) / H, rng)[0] if abs(w) > 1e-8 * np.linalg.norm(q) / H]
        row = []
        for ps in SIGNED_PERMS:
            qp = apply_perm(ps, q)
            rr = [w for w, _ in red_roots(K.pencil(qp), np.linalg.norm(qp) / H, rng)[0] if abs(w) > 1e-8 * np.linalg.norm(q) / H]
            row.append(match_sets(base, rr))
        sym.append(row)
    sym = np.array(sym)
    worst = sym.max(axis=0)
    group = [i for i in range(len(SIGNED_PERMS)) if worst[i] <= 1e-10]
    res['symmetry'] = dict(n_tested=len(SIGNED_PERMS), worst_rel_mismatch=[float(v) for v in worst],
                           perms=[[list(p), list(s)] for p, s in SIGNED_PERMS],
                           group=[[list(SIGNED_PERMS[i][0]), list(SIGNED_PERMS[i][1])] for i in group],
                           group_order=len(group))
    print('symmetry group order', len(group), 'of 48; worst mismatch among non-members',
          float(min([worst[i] for i in range(48) if i not in group] or [np.nan])), flush=True)
    gperms = [SIGNED_PERMS[i] for i in group]
    def canon(d):
        imgs = [apply_perm(ps, d) for ps in gperms]
        return max(imgs, key=lambda v: tuple(np.round(v, 12)))
    # 2. directions
    dirs = {'100': (1, 0, 0), '010': (0, 1, 0), '001': (0, 0, 1)}
    for a, b in ((0, 1), (0, 2), (1, 2)):
        for sg in (1, -1):
            v = [0, 0, 0]; v[a] = 1; v[b] = sg; dirs['d%d%d%s' % (a, b, '+' if sg > 0 else '-')] = tuple(v)
    for s1 in (1, -1):
        for s2 in (1, -1):
            dirs['b1%s%s' % ('+' if s1 > 0 else '-', '+' if s2 > 0 else '-')] = (1, s1, s2)
    dirs['generic'] = (0.8, 0.5, 0.33); dirs['near100'] = (1, 0.05, 0)
    rd = np.random.default_rng(SEED + 2)
    for i in range(30):
        dirs['rand%02d' % i] = tuple(canon(unit(rd.standard_normal(3))))
    approach = {}
    for th in (1e-1, 1e-2, 1e-3, 1e-4):
        approach['A100plane_%g' % th] = (1, np.tan(th), 0)
        approach['A100off_%g' % th] = (1, np.tan(th) * np.cos(0.3), np.tan(th) * np.sin(0.3))
        v110 = unit((1, 1, 0)); approach['A110off_%g' % th] = tuple(np.cos(th) * v110 + np.sin(th) * np.array([0, 0, 1.0]))
    res['directions'] = {k_: list(map(float, unit(v))) for k_, v in dirs.items()}
    res['approach'] = {k_: list(map(float, unit(v))) for k_, v in approach.items()}
    cases = []
    def run(name, d, kh):
        q = kh * unit(d)
        P = K.pencil(q)
        rows, nk = analyse(P, q, rng)
        sz = static_zero_modes(P, q)
        nfull = normal_rank_deficit(P, kh / H, rows=list(range(NV)), cols=list(range(NV)))
        c = dict(dir=name, kh=kh, nk=nk, nullity_full11=nfull, static=sz, roots=rows,
                 light_speed=float(E.light_speed(q)))
        cases.append(c)
        return c
    for kh in (0.0125, 0.025, 0.05, 0.1, 0.2):
        for name, d in dirs.items():
            run(name, d, kh)
        print('kh', kh, 'done', time.process_time() - t0, flush=True)
    for name, d in approach.items():
        run(name, d, 0.05)
    for name, d in (('100', (1, 0, 0)), ('d01+', (1, 1, 0)), ('b1++', (1, 1, 1))):
        for kh in (0.5, 1.0, 2.0, 3.0):
            run(name, d, kh)
        qd = unit(d); khzb = np.pi / np.abs(qd).max()
        run(name, d, float(khzb))
    res['cases'] = cases
    # summaries
    def summ(c):
        g = [r for r in c['roots'] if r['cls'] == 'growing']
        p = [r for r in c['roots'] if r['cls'] == 'propagating']
        return g, p
    allroots = [r for c in cases for r in c['roots'] if r['cls'] != 'zero']
    res['max_resid11'] = max(r['resid11'] for r in allroots)
    res['max_contrast'] = max(r['contrast'] for r in allroots)
    bad = [(c['dir'], c['kh'], r['w_over_k'], r['resid11'], r['contrast']) for c in cases for r in c['roots']
           if r['cls'] != 'zero' and (r['resid11'] > 1e-10 or r['contrast'] > 1e-6)]
    res['failed_validation'] = bad
    # F1(a) vs review_unstable
    rv = json.loads((ROOT / 'results/EX-GC34E/review_unstable.json').read_text())['imag_axis']
    mapdir = {'110': 'd01+', '111': 'b1++', 'generic': 'generic', 'near100': 'near100'}
    cmp_ = []
    for a in rv:
        c = next(c for c in cases if c['dir'] == mapdir[a['dir']] and abs(c['kh'] - a['kh']) < 1e-12)
        g = [r['w_over_k'][1] for r in c['roots'] if r['cls'] == 'growing']
        cmp_.append(dict(dir=a['dir'], kh=a['kh'], review=a['growth_Im_w_over_k'], f1=g,
                         diff=float(min(abs(x - a['growth_Im_w_over_k']) for x in g)) if g else None))
    res['F1a_vs_review'] = cmp_
    # F1(d) propagating count at kh <= 0.05
    d_rows = []
    for c in cases:
        if c['kh'] <= 0.05 + 1e-12 and c['dir'] in dirs:
            p = [r for r in c['roots'] if r['cls'] == 'propagating' and 0.5 < r['w_over_k'][0] < 2]
            d_rows.append(dict(dir=c['dir'], kh=c['kh'], n_prop=len(p),
                               max_dc=float(max([abs(r['w_over_k'][0] - c['light_speed']) for r in p] or [np.nan]))))
    res['F1d'] = d_rows
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, scipy=sl.__name__ and __import__('scipy').__version__,
                           system=platform.system())
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'spectrum.json').write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))
    print('F1 cpu', res['cpu_s'])


# ---------------- F1 post hoc (added after the registered F1 run; see the note) ----------------
# The registered run classified zero roots with |w| <= 1e-8 k, but the zero-shift gauge has a 6-fold zero root
# (3 static displacements, each with a chain e = t G0 that violates only the momentum constraint), which splits
# numerically to |w| ~ 1e-6..1e-3 k; and the singular-value contrast check rejected roots with residual ~1e-15.
# Post hoc: zero-cluster multiplicity from the winding number of a rank-completed determinant around w = 0, and
# every root confirmed by the same winding count around it (argument principle), a computation separate from the
# eigensolver.
def det_fn(P, nk, rng):
    comps = []
    for _ in range(2):
        if nk:
            U = rng.standard_normal((len(RED), nk)) + 1j * rng.standard_normal((len(RED), nk))
            V = rng.standard_normal((len(RED), nk)) + 1j * rng.standard_normal((len(RED), nk))
        else:
            U = V = None
        comps.append((U, V))
    def f(w, j):
        M = Mw(P, w)[np.ix_(RED, RED)]
        if comps[j][0] is not None:
            M = M + np.linalg.norm(M, 2) * comps[j][0] @ comps[j][1].conj().T
        sgn, _ = np.linalg.slogdet(M)
        return sgn
    return f


def winding(f, j, c, r, npt=96):
    ph = np.array([f(c + r * np.exp(2j * np.pi * t / npt), j) for t in range(npt + 1)])
    d = np.angle(ph[1:] / ph[:-1])
    return int(np.round(d.sum() / (2 * np.pi))), float(np.abs(d).max())


def analyse_posthoc(P, q, rng, keep_vectors=False):
    k = np.linalg.norm(q) / H
    roots, nk = red_roots(P, k, rng)
    f = det_fn(P, nk, rng)
    r0 = 1e-2 * k
    m0 = [winding(f, j, 0.0, r0)[0] for j in (0, 1)]
    roots = sorted(roots, key=lambda r: abs(r[0]))
    nz = sum(abs(w) < r0 for w, _ in roots)
    out = dict(nk=nk, zero_winding=m0, zero_found=nz, roots=[])
    ws = np.array([w for w, _ in roots])
    for i, (w, x) in enumerate(roots):
        if abs(w) < r0:
            continue
        others = np.delete(ws, i)
        dmin = np.min(np.abs(others - w)) if len(others) else abs(w)
        rad = min(0.25 * dmin, 1e-3 * abs(w)) if dmin > 0 else 1e-3 * abs(w)
        if dmin < 1e-9 * abs(w):            # numerically coincident pair: circle around both
            rad = 1e-6 * abs(w)
        inside = int((np.abs(ws - w) < rad).sum())
        wn = [winding(f, j, w, rad) for j in (0, 1)]
        res, _, mom = validate(P, k, w, x, nk)
        ngf, ttf = gauge_split(x, q, w)
        rel_im = w.imag / abs(w)
        cls = ('propagating' if w.real > 0 else 'propagating_neg') if abs(rel_im) <= 1e-8 else ('growing' if w.imag > 0 else 'decaying')
        row = dict(cls=cls, w_over_k=[w.real / k, w.imag / k], wh=[w.real * H, w.imag * H], resid11=res, mom_resid=mom,
                   winding=[wn[0][0], wn[1][0]], roots_inside=inside, max_phase_step=max(wn[0][1], wn[1][1]),
                   confirmed=bool(res <= 1e-10 and wn[0][0] == inside and wn[1][0] == inside),
                   non_gauge_frac=ngf, tt_frac_of_nongauge=ttf, abs_n=float(abs(x[IN])), abs_edges=float(np.linalg.norm(x[IE])))
        if keep_vectors:
            row['_x'] = x; row['_w'] = w
        out['roots'].append(row)
    return out


def f1_posthoc():
    t0 = time.process_time()
    mp.mp.dps = 50
    st = E.tet_static(); K = Kernel(st); rng = np.random.default_rng(SEED + 10)
    reg = json.loads((OUT / 'spectrum.json').read_text())
    res = dict(label='POST HOC (after the registered F1 run); not a re-scoring', convention=reg['convention'])
    # symmetry test with confirmed, non-zero roots only
    rs = np.random.default_rng(SEED + 1)
    worst = np.zeros(len(SIGNED_PERMS)); nconf_fail = 0
    for i in range(20):
        q = unit(rs.standard_normal(3)) * rs.uniform(0.05, 0.6)
        a = analyse_posthoc(K.pencil(q), q, rng)
        base = [complex(*r['w_over_k']) for r in a['roots'] if r['confirmed']]
        nconf_fail += sum(not r['confirmed'] for r in a['roots'])
        for j, ps in enumerate(SIGNED_PERMS):
            qp = apply_perm(ps, q)
            b = analyse_posthoc(K.pencil(qp), qp, rng)
            worst[j] = max(worst[j], match_sets(base, [complex(*r['w_over_k']) for r in b['roots'] if r['confirmed']]))
    group = [i for i in range(48) if worst[i] <= 1e-10]
    res['symmetry'] = dict(worst_rel_mismatch=[float(v) for v in worst], perms=reg['symmetry']['perms'],
                           group=[[list(SIGNED_PERMS[i][0]), list(SIGNED_PERMS[i][1])] for i in group], group_order=len(group),
                           unconfirmed_roots_in_base_sets=nconf_fail)
    print('post hoc symmetry group order', len(group), flush=True)
    cases = []
    allq = [(c['dir'], c['kh']) for c in reg['cases']]
    dirs = dict(reg['directions']); dirs.update(reg['approach'])
    for name, kh in allq:
        q = kh * np.array(dirs[name])
        P = K.pencil(q)
        a = analyse_posthoc(P, q, rng)
        a.update(dir=name, kh=kh, light_speed=float(E.light_speed(q)), static=static_zero_modes(P, q),
                 nullity_full11=normal_rank_deficit(P, kh / H, rows=list(range(NV)), cols=list(range(NV))))
        cases.append(a)
    res['cases'] = cases
    allr = [r for c in cases for r in c['roots']]
    res['n_roots'] = len(allr); res['n_unconfirmed'] = sum(not r['confirmed'] for r in allr)
    res['unconfirmed'] = [(c['dir'], c['kh'], r['w_over_k'], r['resid11'], r['winding'], r['roots_inside'])
                          for c in cases for r in c['roots'] if not r['confirmed']]
    res['zero_cluster'] = sorted({(c['zero_winding'][0], c['zero_winding'][1], c['zero_found'], c['nk']) for c in cases})
    rv = json.loads((ROOT / 'results/EX-GC34E/review_unstable.json').read_text())['imag_axis']
    mapdir = {'110': 'd01+', '111': 'b1++', 'generic': 'generic', 'near100': 'near100'}
    cmp_ = []
    for a in rv:
        c = next(c for c in cases if c['dir'] == mapdir[a['dir']] and abs(c['kh'] - a['kh']) < 1e-12)
        g = [r['w_over_k'][1] for r in c['roots'] if r['cls'] == 'growing' and r['confirmed']]
        cmp_.append(dict(dir=a['dir'], kh=a['kh'], review=a['growth_Im_w_over_k'], f1=g,
                         diff=float(min(abs(x - a['growth_Im_w_over_k']) for x in g)) if g else None))
    res['F1a_vs_review'] = cmp_
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    res['cpu_s'] = time.process_time() - t0
    (OUT / 'spectrum_posthoc.json').write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))
    print('F1 post hoc cpu', res['cpu_s'], 'unconfirmed', res['n_unconfirmed'], 'of', res['n_roots'])


# ---------------- F2 ----------------
def reduced_system(P, q):
    A = P['A'][np.ix_(IE, IE)]; C = P['C']
    Ai = np.linalg.inv(A)
    Cee = C[np.ix_(IE, IE)]; Cen = C[IE, IN]; Cne = C[IN, IE]
    s = Cne @ Ai @ Cen
    Pi = -(Cne @ Ai @ Cee) / s                     # n = Pi . e
    F = Ai @ (Cee + np.outer(Cen, Pi))              # e'' = F e
    return dict(F=F, Pi=Pi, Cne=Cne, CbE=C[np.ix_(IB, IE)], Cbn=C[IB, IN], BbE=P['B'][np.ix_(IB, IE)])


def tidal(rs, q, e, edd):
    n = rs['Pi'] @ e
    x = np.zeros(NV, complex); x[IN] = n; x[IE] = e
    h00 = E.metric_from_x(x, q)[0, 0]
    xd = np.zeros(NV, complex); xd[IE] = edd
    hdd = E.metric_from_x(xd, q)[1:, 1:]
    qq = np.outer(q, q) / H ** 2
    return 0.5 * (hdd - qq * h00)


def evolve(rs, q, z0, dt, nsteps):
    n7 = len(IE); F = rs['F']
    L = np.block([[np.zeros((n7, n7)), np.eye(n7)], [F, np.zeros((n7, n7))]])
    I = np.eye(2 * n7)
    Mstep = np.linalg.solve(I - 0.5 * dt * L, I + 0.5 * dt * L)
    z = z0.copy(); out = dict(t=[], z=[], lapse=[], mom=[], E=[])
    nCne = np.linalg.norm(rs['Cne']); nCb = np.linalg.norm(rs['CbE'], 2); nBb = np.linalg.norm(rs['BbE'], 2)
    for i in range(nsteps + 1):
        e = z[:n7]; ed = z[n7:]; n = rs['Pi'] @ e
        out['t'].append(i * dt); out['z'].append(z.copy())
        out['lapse'].append(float(abs(rs['Cne'] @ e) / (nCne * np.linalg.norm(e))))
        mr = rs['CbE'] @ e + rs['Cbn'] * n + 1j * rs['BbE'] @ ed
        out['mom'].append(float(np.linalg.norm(mr) / (nCb * np.linalg.norm(e) + nBb * np.linalg.norm(ed))))
        out['E'].append(tidal(rs, q, e, F @ e))
        z = Mstep @ z
    return out


def slope(t, y):
    return float(np.polyfit(np.asarray(t), np.log(np.asarray(y)), 1)[0])


def f2():
    t0 = time.process_time()
    mp.mp.dps = 50
    st = E.tet_static(); K = Kernel(st)
    spec = json.loads((OUT / 'spectrum_posthoc.json').read_text())   # confirmed roots (see note)
    res = dict(convention='e^{-i w t}; implicit midpoint on z = (e, de/dt), zero-shift gauge, index-reduced lapse', runs=[])
    for name, d in (('generic', (0.8, 0.5, 0.33)), ('d01+', (1, 1, 0))):
        kh = 0.05; q = kh * unit(d); k = kh / H
        P = K.pencil(q); rs = reduced_system(P, q)
        lam, V = np.linalg.eig(-rs['F'])           # W v = w^2 v
        ig = int(np.argmin(lam.real)); gam = float(np.sqrt(-lam[ig].real))
        cand = [i for i in range(len(lam)) if lam[i].real > 0.25 * k * k and abs(lam[i].imag) < 1e-8 * k * k]
        iw = min(cand, key=lambda i: abs(np.sqrt(lam[i].real) / k - 1)); om = float(np.sqrt(lam[iw].real))
        c = next(c for c in spec['cases'] if c['dir'] == name and abs(c['kh'] - kh) < 1e-12)
        f1g = [r['w_over_k'][1] for r in c['roots'] if r['cls'] == 'growing' and r['confirmed']]
        f1p = [r['w_over_k'][0] for r in c['roots'] if r['cls'] == 'propagating' and r['confirmed']]
        vg = V[:, ig] / np.linalg.norm(V[:, ig]); vw = V[:, iw] / np.linalg.norm(V[:, iw])
        zg = np.r_[vg, gam * vg]; zw = np.r_[vw, -1j * om * vw]
        # sanity: tidal field of eigenvectors vs riemann()
        xg = np.zeros(NV, complex); xg[IN] = rs['Pi'] @ vg; xg[IE] = vg
        Rg = E.riemann(E.metric_from_x(xg, q), np.r_[-1j * gam, q / H])
        comp = [i for i, (a, b, cc, dd) in enumerate((a, b, cc, dd) for a in range(4) for b in range(a + 1, 4)
                                                    for cc in range(4) for dd in range(cc + 1, 4)) if a == 0 and cc == 0]
        Eg = tidal(rs, q, vg, rs['F'] @ vg)
        Rg0 = np.array([Rg[i] for i in comp])
        Eg0 = np.array([Eg[b - 1, dd - 1] for (a, b, cc, dd) in ((a, b, cc, dd) for a in range(4) for b in range(a + 1, 4)
                                                              for cc in range(4) for dd in range(cc + 1, 4)) if a == 0 and cc == 0])
        sanity = float(np.linalg.norm(Eg0 - Rg0) / np.linalg.norm(Rg0))
        T = 15.0 / gam
        runs = {}
        for tag, dt0 in (('dt', 0.02 / max(om, gam)), ('dt2', 0.01 / max(om, gam))):
            ns = int(np.ceil(T / dt0)); dt = T / ns
            fac = 0.02 if tag == 'dt' else 0.01
            r1 = evolve(rs, q, zg, T / int(np.ceil(T * gam / fac)), int(np.ceil(T * gam / fac)))
            r2 = evolve(rs, q, zw, dt, ns)
            seed = 1e-8
            r3 = evolve(rs, q, zw + seed * zg, dt, ns)
            # mode amplitudes by eigen-decomposition of the first-order operator
            n7 = len(IE)
            Lop = np.block([[np.zeros((n7, n7)), np.eye(n7)], [rs['F'], np.zeros((n7, n7))]])
            mu, Vf = np.linalg.eig(Lop)
            jg = int(np.argmin(np.abs(mu - gam))); jw = int(np.argmin(np.abs(mu + 1j * om)))
            amp = lambda zz, j: np.linalg.solve(Vf, zz)[j]
            a1 = [abs(amp(z, jg)) for z in r1['z']]
            E1 = [np.linalg.norm(Ei) for Ei in r1['E']]
            a2 = [abs(amp(z, jw)) for z in r2['z']]
            dE = [np.linalg.norm(Ea - Eb) for Ea, Eb in zip(r3['E'], r2['E'])]
            Ew = [np.linalg.norm(Eb) for Eb in r2['E']]
            tt = np.array(r2['t'])
            ratio = np.array(dE) / np.array(Ew)
            j = int(np.argmax(ratio >= 1.0)) if (ratio >= 1.0).any() else None
            if j:
                tx = float(tt[j - 1] + (tt[j] - tt[j - 1]) * (np.log(1) - np.log(ratio[j - 1])) / (np.log(ratio[j]) - np.log(ratio[j - 1])))
            else:
                tx = None
            txp = float(np.log(Ew[0] / dE[0]) / gam)
            runs[tag] = dict(
                dt=dt, steps=ns,
                i_growth_amp=slope(r1['t'], a1), i_growth_E=slope(r1['t'], E1),
                i_max_lapse=max(r1['lapse']), i_max_mom=max(r1['mom']),
                ii_amp_change=float(max(abs(a / a2[0] - 1) for a in a2)),
                ii_max_lapse=max(r2['lapse']), ii_max_mom=max(r2['mom']),
                iii_max_lapse=max(r3['lapse']), iii_max_mom=max(r3['mom']),
                iii_seed_growth_E=slope(tt, dE), t_cross=tx, t_cross_pred=txp,
                t_cross_rel_err=(abs(tx - txp) / txp if tx else None), E_seed0_over_Ewave=float(dE[0] / Ew[0]))
        # POST HOC (added after the registered F2 run): t_cross lies beyond the registered run length 15/Im w,
        # so runs (ii) and (iii) are extended to 40/Im w at the registered step to measure it.
        dtl = 0.02 / max(om, gam); nsl = int(np.ceil(40.0 / gam / dtl)); dtl = 40.0 / gam / nsl
        L2 = evolve(rs, q, zw, dtl, nsl); L3 = evolve(rs, q, zw + 1e-8 * zg, dtl, nsl)
        dEl = np.array([np.linalg.norm(Ea - Eb) for Ea, Eb in zip(L3['E'], L2['E'])]); Ewl = np.array([np.linalg.norm(Eb) for Eb in L2['E']])
        rl = dEl / Ewl; ttl = np.array(L2['t']); jj = int(np.argmax(rl >= 1.0)) if (rl >= 1.0).any() else 0
        txl = float(ttl[jj - 1] + (ttl[jj] - ttl[jj - 1]) * (-np.log(rl[jj - 1])) / (np.log(rl[jj]) - np.log(rl[jj - 1]))) if jj else None
        txpl = float(np.log(Ewl[0] / dEl[0]) / gam)
        posthoc = dict(label='POST HOC extended run, 40/Im w', t_cross=txl, t_cross_pred=txpl,
                       t_cross_rel_err=(abs(txl - txpl) / txpl if txl else None), E_seed0_over_Ewave=float(dEl[0] / Ewl[0]),
                       seed_growth_E=slope(ttl, dEl), max_lapse=max(L3['lapse']), max_mom=max(L3['mom']),
                       t_cross_times_k=(txl * k if txl else None), t_cross_in_wave_periods=(txl * om / (2 * np.pi) if txl else None))
        rr = dict(dir=name, kh=kh, gamma=gam, gamma_over_k=gam / k, omega_wave_over_k=om / k,
                  f1_growing=f1g, f1_propagating=f1p, tidal_sanity_vs_riemann=sanity, T=T, runs=runs, posthoc_long=posthoc)
        res['runs'].append(rr)
        print(name, 'gamma/k %.6f (F1 %s)  wave %.6f' % (gam / k, f1g, om / k), 'sanity %.1e' % sanity, flush=True)
        for tag, r in runs.items():
            print('  ', tag, {kk: ('%.4e' % v if isinstance(v, float) else v) for kk, v in r.items()}, flush=True)
        print('   posthoc', {kk: ('%.4e' % v if isinstance(v, float) else v) for kk, v in posthoc.items()}, flush=True)
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    res['cpu_s'] = time.process_time() - t0
    (OUT / 'evolve.json').write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))
    print('F2 cpu', res['cpu_s'])


# ---------------- F3 ----------------
AXG = [0, 1, 2, 3] + [4 + t for t in range(3, 7)]          # axis-edge gauge (needs all q_a != 0)
BLOCKS = {'kin_bb': (IB, IB), 'kin_be': (IB, IE), 'kin_ee': (IE, IE), 'regge_ne': ([IN], IE), 'regge_ee': (IE, IE)}


def block_parts(P):
    """Split M(w) = sum_X M_X(w); each part as (C, B, A) triple."""
    parts = {}
    T = {'C': P['T0'], 'B': P['B'], 'A': P['A']}
    R = {'C': P['R'], 'B': 0 * P['B'], 'A': 0 * P['A']}
    for name, (r, c) in BLOCKS.items():
        src = T if name.startswith('kin') else R
        m = np.zeros((NV, NV), bool); m[np.ix_(r, c)] = True; m[np.ix_(c, r)] = True
        parts[name] = {x: np.where(m, src[x], 0) for x in ('C', 'B', 'A')}
    tot = {x: sum(parts[n][x] for n in parts) for x in ('C', 'B', 'A')}
    leftover = max(float(np.abs(tot[x] - P[x]).max()) for x in ('C', 'B', 'A')) / float(np.abs(P['C']).max())
    return parts, leftover


def combine(parts, scale):
    return {x: sum(scale.get(n, 1.0) * parts[n][x] for n in parts) for x in ('C', 'B', 'A')}


def roots_idx(P, k, idx, rng):
    C, B, A = [P[x][np.ix_(idx, idx)] for x in ('C', 'B', 'A')]
    n = len(idx); I = np.eye(n); Z = np.zeros((n, n))
    ev, V = sl.eig(np.block([[Z, I], [-C, -B]]), np.block([[I, Z], [Z, A]]))
    out = []
    for e, v in zip(ev, V.T):
        if np.isfinite(e) and abs(e) < 1e6 * k:
            x = np.zeros(NV, complex); x[idx] = v[:n]; out.append((complex(e), x / np.linalg.norm(x)))
    return out


def track(P, k, idx, w_ref, rng):
    rr = roots_idx(P, k, idx, rng)
    w, x = min(rr, key=lambda r: abs(r[0] - w_ref))
    return w, x


def f3():
    t0 = time.process_time()
    mp.mp.dps = 50
    st = E.tet_static(); K = Kernel(st); rng = np.random.default_rng(SEED + 20)
    res = dict(convention='e^{-i w t}', blocks=list(BLOCKS))
    # 1. contributions and sensitivities
    mech = []
    for name, d in (('b1++', (1, 1, 1)), ('generic', (0.8, 0.5, 0.33))):
        for kh in (0.0125, 0.025, 0.05, 0.1, 0.2):
            q = kh * unit(d); k = kh / H; P = K.pencil(q)
            parts, left = block_parts(P)
            row = dict(dir=name, kh=kh, leftover=left)
            for gname, idx in (('axis_edge', AXG), ('zero_shift', RED)):
                w, x = track(P, k, idx, 1.42j * k, rng)
                x = x / x[IN]
                contrib = {n: complex(x.conj() @ Mw(parts[n], w) @ x) for n in parts}
                tot = sum(contrib.values())
                sens = {}
                for n in parts:
                    vals = []
                    for eps in (1e-6, -1e-6):
                        Ps = combine(parts, {n: 1 + eps}); Ps['T0'] = None
                        ww, _ = track(Ps, k, idx, w, rng); vals.append(np.log(-(ww * ww).real))
                    sens[n] = float((vals[0] - vals[1]) / 2e-6)
                keep = [n for n in parts if abs(sens[n]) >= 0.05]
                Pm = combine(parts, {n: (1.0 if n in keep else 0.0) for n in parts})
                try:
                    wm, _ = track(Pm, k, idx, w, rng); mrel = float(abs(wm - w) / abs(w))
                except ValueError:
                    wm, mrel = None, None
                row[gname] = dict(w_over_k=[w.real / k, w.imag / k],
                                  contrib_over_k2={n: [v.real / k ** 2, v.imag / k ** 2] for n, v in contrib.items()},
                                  contrib_sum_rel=float(abs(tot) / max(abs(v) for v in contrib.values())),
                                  sens=sens, minimal_keep=keep,
                                  minimal_w_over_k=None if wm is None else [wm.real / k, wm.imag / k], minimal_rel_err=mrel,
                                  n_over_beta_par=float(abs(x[IN]) / max(abs(unit(d) @ x[IB]), 1e-300)),
                                  edges_over_n=float(np.linalg.norm(x[IE])))
            mech.append(row)
            print(name, kh, 'AXG sens', {n: round(v, 4) for n, v in row['axis_edge']['sens'].items()},
                  'min err', row['axis_edge']['minimal_rel_err'], '| ZS sens', {n: round(v, 4) for n, v in row['zero_shift']['sens'].items()},
                  'min err', row['zero_shift']['minimal_rel_err'], flush=True)
    res['mechanism'] = mech
    # exponents of contributions (axis-edge gauge), fixed physical k: fit |c|/k^2 vs kh
    expo = {}
    for name in ('b1++', 'generic'):
        rows = [r for r in mech if r['dir'] == name]
        khs = np.array([r['kh'] for r in rows])
        for g in ('axis_edge', 'zero_shift'):
            for n in BLOCKS:
                c = np.array([abs(complex(*r[g]['contrib_over_k2'][n])) for r in rows])
                expo['%s/%s/%s' % (name, g, n)] = float(np.polyfit(np.log(khs), np.log(c + 1e-300), 1)[0]) if (c > 0).all() else None
    res['contrib_exponent_in_kh'] = expo
    # 2 + 3. constraint bookkeeping and diagnostic variants
    variants = {'base': Kernel(st), 'a_cube': Kernel(st, kin='cube'), 'b_lambda_0.5': Kernel(st, lam=0.5),
                'b_lambda_1/3': Kernel(st, lam=1 / 3), 'c_no_lapse_regge': Kernel(st, lapse_regge=False)}
    K34 = D.kernel(D.build_terms()[0], np.array([[0.3, -0.7, 1.1]]))[0]
    sel = [0] + list(range(4, 11))
    vres = {}
    for vn, KV in variants.items():
        out = dict()
        qz = np.array([0.3, -0.7, 1.1])
        out['static_vs_34D'] = float(np.abs(KV.M(qz, 0.0)[np.ix_(sel, sel)] - K34).max() / np.abs(K34).max())
        # displacement-gauge null check
        out['gauge_null_rel'] = float(max(np.linalg.norm(KV.M(qz, w) @ E.gauge_modes(qz[None], w)[0]) / np.linalg.norm(KV.M(qz, w))
                                          for w in (0.0, 0.7, -1.3)))
        sp = []
        for name, d in (('b1++', (1, 1, 1)), ('generic', (0.8, 0.5, 0.33)), ('d01+', (1, 1, 0)), ('100', (1, 0, 0))):
            for kh in (0.05, 0.1):
                q = kh * unit(d); k = kh / H; P = KV.pencil(q)
                a = analyse_posthoc(P, q, rng)
                cons = dict(nullity_full11=normal_rank_deficit(P, k, rows=list(range(NV)), cols=list(range(NV))), nk=a['nk'],
                            rank_lapse_shift_rows=int(np.linalg.matrix_rank(Mw(P, 0.81 * k)[[0] + IB, :], tol=1e-12 * np.abs(Mw(P, 0.81 * k)).max())))
                G0 = E.gauge_modes(q[None], 0.0)[0]
                brk = []
                for w in (0.0, 0.5 * k):
                    Mm = Mw(P, w); G = E.gauge_modes(q[None], w)[0]
                    Qc = sl.null_space(G.conj().T); sv = np.linalg.svd(Qc.conj().T @ Mm @ Qc, compute_uv=False)
                    brk.append(float(sv[-1] / sv[0]))
                sp.append(dict(dir=name, kh=kh, constraints=cons, reslicing_smin=brk,
                               growing=[r['w_over_k'] for r in a['roots'] if r['cls'] == 'growing' and r['confirmed']],
                               propagating=[r['w_over_k'][0] for r in a['roots'] if r['cls'] == 'propagating' and r['confirmed'] and r['w_over_k'][0] < 20],
                               other=[r['w_over_k'] for r in a['roots'] if r['cls'] in ('growing',) and not r['confirmed']],
                               n_unconfirmed=sum(not r['confirmed'] for r in a['roots'])))
        out['spectra'] = sp
        # E4 gravitomagnetism at kh = 0.05 (34E case list)
        orig = E.full_kernel
        E.full_kernel = lambda st_, qq, w, KV=KV: KV.M(qq[0], w)[None]
        try:
            e4 = []
            for qdir in ((1, 0, 0), (1, 1, 1), (0.8, 0.5, 0.33)):
                qd = unit(qdir)
                perp = np.cross(qd, [0, 0, 1.0]) if abs(qd[2]) < 0.9 else np.cross(qd, [1.0, 0, 0])
                for vname, vdir in (('par', qd), ('perp', perp), ('obl', qd + 1.3 * perp)):
                    r = E.e4_case(st, qdir, 0.05, vdir); e4.append(dict(qdir=list(qdir), v=vname, err_Ov=r['err_Ov'], err_static=r['err_static'],
                                                                         momentum_residual=r['momentum_residual']))
        finally:
            E.full_kernel = orig
        out['E4_max_err_Ov'] = max(r['err_Ov'] for r in e4); out['E4_max_err_static'] = max(r['err_static'] for r in e4); out['E4'] = e4
        vres[vn] = out
        print(vn, 'static %.1e gauge %.1e E4 %.2e/%.2e' % (out['static_vs_34D'], out['gauge_null_rel'], out['E4_max_err_static'], out['E4_max_err_Ov']),
              [(r['dir'], r['kh'], [round(g[1], 4) for g in r['growing']], r['constraints']['nullity_full11'], r['constraints']['nk']) for r in sp], flush=True)
    res['variants'] = vres
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    res['cpu_s'] = time.process_time() - t0
    (OUT / 'mechanism.json').write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))
    print('F3 cpu', res['cpu_s'])


if __name__ == '__main__':
    {'f1': f1, 'f1p': f1_posthoc, 'f2': f2, 'f3': f3}[sys.argv[1]]()
