"""EX-GC34L (exploratory): is J1's first-order E4 error in the moving source, the curvature readout, or the
geometry? Plan: spec/EX-GC34L_plan.md.

  python src/exgc34l.py audit  -> results/EX-GC34L/audit.json    (A: source overlaps, forced residual, readout gauge sensitivity)
  python src/exgc34l.py bench  -> results/EX-GC34L/bench.json    (B: S0/S1 x R0/R1/R1w/R2; registered S1+R1)
  python src/exgc34l.py fit    -> results/EX-GC34L/fit.json, kernel_C.npz   (C: only if L-FAIL)
Reuses exgc34e/f/g/i/k unchanged.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import json, platform, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp
import scipy.optimize as so

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import exgc34e as E
import exgc34f as F
import exgc34g as G
import exgc34i as I
import exgc34k as K

OUT = ROOT / 'results/EX-GC34L'
H = E.H; NV = 11; MU = K.MU; KAP = E.KAPPA
KHS = (0.0125, 0.025, 0.05, 0.1)
VMAG = 1e-3
DUMP = G.DUMP
SEED = 3434 + 12
gfwd = K.g_gamma(0.0)


# ---------------- gauge directions, sources, readouts ----------------
def Wmat(q, w):
    z = np.zeros(NV, complex); z[0] = -1j * w; z[1:4] = MU * gfwd(q)
    return np.hstack([E.gauge_modes(q[None], w)[0], z[:, None]])


def source(kind, q, v):
    f = np.zeros(NV, complex); f[0] = KAP
    if kind == 'S0':
        f[1:4] = -KAP * v
    else:   # S1: link current from exact lattice continuity
        sinc = np.array([1.0 if abs(qa) < 1e-14 else (qa / 2) / np.sin(qa / 2) for qa in q])
        f[1:4] = -KAP * v * np.exp(0.5j * q) * sinc
    return f


def L0mat(q):
    return np.stack([E.metric_from_x(np.eye(NV)[j].astype(complex), q).reshape(16) for j in range(NV)], 1)


def Gcmat(kl):
    out = np.zeros((16, 4), complex)
    for nu in range(4):
        h = np.zeros((4, 4), complex)
        for m in range(4):
            h[m, nu] += 1j * kl[m]; h[nu, m] += 1j * kl[m]
        out[:, nu] = h.reshape(16)
    return out


def readout(kind, q, w):
    """Returns a 16 x 11 matrix x -> h (flattened 4x4)."""
    L0 = L0mat(q)
    if kind == 'R0':
        return L0
    if kind == 'R2':
        P = np.eye(NV, dtype=complex); P[1:4, 1:4] = np.diag(np.exp(-0.5j * q))
        return L0 @ P
    wts = {'R1': (1, 1, 1), 'R1w_a': (2, 1, 1), 'R1w_b': (1, 1, 2)}[kind]
    Dg = np.diag([wts[0]] + [wts[1]] * 3 + [wts[2]] * 7).astype(complex)
    W = Wmat(q, w); kl = np.array([-w, *(q / H)]); Gc = Gcmat(kl)
    A = np.linalg.lstsq(Gc, L0 @ W, rcond=None)[0]
    Wp = np.linalg.solve(W.conj().T @ Dg @ W, W.conj().T @ Dg)
    return L0 + (Gc @ A - L0 @ W) @ Wp


def riem(L, x, q, w):
    return E.riemann((L @ x).reshape(4, 4), np.array([-w, *(q / H)])) * H ** 3


def solve(kin, st, q, w, f):
    Mm = kin(q, w) + E.static_kernel(st, q[None])[0]; W = Wmat(q, w)
    x, *_ = np.linalg.lstsq(Mm + W @ W.conj().T, f, rcond=1e-14)
    return x, Mm, W


def cases():
    return K.CASES


def e4_rows(kin, st, kh, src, rd, keep_x=False):
    out = []
    for qdir, vn, vd in cases():
        qd = G.unit(qdir); q = kh * qd; k = q / H
        v = VMAG * np.array(vd) / np.linalg.norm(vd)
        R = {}; extra = {}
        for sgn in (1, -1):
            w = k @ (sgn * v); f = source(src, q, sgn * v)
            x, Mm, W = solve(kin, st, q, w, f)
            Rl = riem(readout(rd, q, w), x, q, w)
            hg, klg, _ = E.gr_solution(k, sgn * v); Rg = E.riemann(hg, klg)
            R[sgn] = (Rl, Rg)
            if sgn == 1:
                extra = dict(resid=float(np.linalg.norm(Mm @ x - f) / np.linalg.norm(f)))
                if keep_x:
                    extra.update(x=x, w=w, q=q, Rg=Rg)
        R1l = (R[1][0] - R[-1][0]) / 2; R1g = (R[1][1] - R[-1][1]) / 2
        R0l = (R[1][0] + R[-1][0]) / 2; R0g = (R[1][1] + R[-1][1]) / 2
        d = (R1l - R1g) / np.abs(R1g).max()
        out.append(dict(qdir=list(qdir), v=vn, comp_err=d, err_Ov=float(np.abs(d).max()),
                        err_static=float(np.abs(R0l - R0g).max() / np.abs(R0g).max()), R1g_max=float(np.abs(R1g).max()), **extra))
    return out


def e4_max(rows):
    return max(max(r['err_Ov'], r['err_static']) for r in rows)


# ---------------- Part A ----------------
def audit():
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); kin = I.make_kin(K.load_J1()); res = dict(overlaps=[], sensitivity=[])
    for kh in KHS:
        for src in ('S0', 'S1'):
            rows = e4_rows(kin, st, kh, src, 'R0')
            for (qdir, vn, vd), r in zip(cases(), rows):
                qd = G.unit(qdir); q = kh * qd; v = VMAG * np.array(vd) / np.linalg.norm(vd); w = (q / H) @ v
                f = source(src, q, v); W = Wmat(q, w)
                ovz = abs(W[:, 3].conj() @ f) / (KAP * np.linalg.norm(q / H) * VMAG)   # = |w kappa| for v || q; defined when w = 0
                ovg = float(np.abs(W[:, :3].conj().T @ f).max() / np.linalg.norm(f))
                res['overlaps'].append(dict(kh=kh, src=src, qdir=list(qdir), v=vn, reslicing_over_omega=float(ovz),
                                            displacement_rel=ovg, forced_resid=r['resid'], err_Ov_R0=r['err_Ov'],
                                            comp_err_R0=[float(abs(c)) for c in r['comp_err']]))
        # readout gauge sensitivity, worst case of 34K: q || (1,0,0), v parallel
        q = kh * np.array([1.0, 0, 0]); v = VMAG * np.array([1.0, 0, 0]); w = (q / H) @ v
        x, Mm, W = solve(kin, st, q, w, source('S1', q, v))
        hg, klg, _ = E.gr_solution(q / H, v); Rg = E.riemann(hg, klg)
        hgm, klgm, _ = E.gr_solution(q / H, -v); R1g = (Rg - E.riemann(hgm, klgm)) / 2
        for rd in ('R0', 'R1', 'R1w_a', 'R1w_b', 'R2'):
            L = readout(rd, q, w); base = riem(L, x, q, w); ch = []
            for j in range(4):
                xg = x + W[:, j] * np.linalg.norm(x) / np.linalg.norm(W[:, j])
                ch.append(float(np.abs(riem(L, xg, q, w) - base).max()))
            res['sensitivity'].append(dict(kh=kh, readout=rd, change_over_Rfull=[c / np.abs(Rg).max() for c in ch],
                                           change_over_ROv=[c / np.abs(R1g).max() for c in ch]))
        print('kh', kh, 'S0 ovz max %.2e S1 ovz max %.2e' % (
            max(o['reslicing_over_omega'] for o in res['overlaps'] if o['kh'] == kh and o['src'] == 'S0'),
            max(o['reslicing_over_omega'] for o in res['overlaps'] if o['kh'] == kh and o['src'] == 'S1')),
            [(s['readout'], '%.1e' % max(s['change_over_Rfull'])) for s in res['sensitivity'] if s['kh'] == kh], flush=True)
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'audit.json').write_text(json.dumps(res, **DUMP))
    print('cpu', res['cpu_s'])


# ---------------- Part B ----------------
def bench():
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); kin = I.make_kin(K.load_J1()); res = dict(table={})
    for src in ('S0', 'S1'):
        for rd in ('R0', 'R1', 'R1w_a', 'R1w_b', 'R2'):
            rows = {kh: e4_rows(kin, st, kh, src, rd) for kh in KHS}
            e = {kh: e4_max(rows[kh]) for kh in KHS}
            i, c = K.worst(rows[0.05])
            res['table']['%s_%s' % (src, rd)] = dict(E4={str(k): v for k, v in e.items()}, ratio=e[0.05] / e[0.025],
                                                     ratio_0125=e[0.025] / e[0.0125],
                                                     worst_case=[rows[0.05][i]['qdir'], rows[0.05][i]['v']], worst_comp=c,
                                                     per_case_005=[r['err_Ov'] for r in rows[0.05]],
                                                     max_resid=max(r['resid'] for kh in KHS for r in rows[kh]))
            t = res['table']['%s_%s' % (src, rd)]
            print(src, rd, ' '.join('%.3e' % e[k] for k in KHS), 'ratio %.3f' % t['ratio'], 'worst', t['worst_case'], c,
                  'resid %.1e' % t['max_resid'], flush=True)
    p = res['table']['S1_R1']
    res['registered'] = dict(E4_005=p['E4']['0.05'], ratio=p['ratio'],
                             verdict='L-PASS' if (p['E4']['0.05'] <= 1e-2 and p['ratio'] >= 3) else 'L-FAIL')
    res['cpu_s'] = time.process_time() - t0
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'bench.json').write_text(json.dumps(res, **DUMP))
    print('registered', res['registered'], 'cpu', res['cpu_s'])


# ---------------- Part C ----------------
def cvec(coef, st, kh=0.0125):
    rows = e4_rows(I.make_kin(coef), st, kh, 'S1', 'R1')
    c = np.concatenate([r['comp_err'] for r in rows]) / kh
    return np.concatenate([c.real, c.imag])


def fit():
    t0 = time.process_time(); mp.mp.dps = 50
    b = json.loads((OUT / 'bench.json').read_text())
    if b['registered']['verdict'] != 'L-FAIL':
        print('L-PASS: Part C not run'); return
    st = E.tet_static(); rng = np.random.default_rng(SEED); res = {}; save = {}
    for R in (1, 2):
        info, x0, Nsp, Wt, xT = K.solve_space(gfwd, R, st)
        d = np.linalg.norm(x0 - xT)
        # J1 itself is the range-1 closest-to-T solution; at range 2 the base is the range-2 closest-to-T solution
        base = x0
        t1 = time.process_time(); c0 = cvec(I.coef_from_vec(base / Wt, R), st); t_one = time.process_time() - t1
        dim = Nsp.shape[1]
        if R == 2 and t_one * dim > 900:
            Q, _ = np.linalg.qr(rng.standard_normal((dim, 300))); B = Nsp @ Q; sub = 300
        else:
            B = Nsp; sub = dim
        step = 1e-3 * d
        Jm = np.zeros((c0.size, B.shape[1]))
        for j in range(B.shape[1]):
            Jm[:, j] = (cvec(I.coef_from_vec((base + step * B[:, j]) / Wt, R), st) - c0) / step
        nb = B.shape[1]; rad = 0.5 * d
        # minimise t s.t. |c0 + J y|_inf <= t, |y|_2 <= rad  (B has orthonormal columns)
        cons = [dict(type='ineq', fun=lambda p: np.concatenate([p[-1] - (c0 + Jm @ p[:-1]), p[-1] + (c0 + Jm @ p[:-1])]),
                     jac=lambda p: np.vstack([np.hstack([-Jm, np.ones((len(c0), 1))]), np.hstack([Jm, np.ones((len(c0), 1))])])),
                dict(type='ineq', fun=lambda p: np.array([rad ** 2 - p[:-1] @ p[:-1]]),
                     jac=lambda p: np.hstack([-2 * p[:-1], [0.0]])[None])]
        p0 = np.r_[np.zeros(nb), np.abs(c0).max()]
        sol = so.minimize(lambda p: p[-1], p0, jac=lambda p: np.r_[np.zeros(nb), 1.0], constraints=cons, method='SLSQP',
                          options=dict(maxiter=500, ftol=1e-12))
        y = sol.x[:-1]; lin_max = float(np.abs(c0 + Jm @ y).max())
        coef = I.coef_from_vec((base + B @ y) / Wt, R); kin = I.make_kin(coef)
        e = {kh: e4_max(e4_rows(kin, st, kh, 'S1', 'R1')) for kh in KHS}
        out = dict(R=R, dim=dim, subspace=sub, c0_max=float(np.abs(c0).max()), lin_min_max=lin_max, opt_success=bool(sol.success),
                   opt_msg=str(sol.message), y_norm_over_d=float(np.linalg.norm(y) / d),
                   E4={str(k): v for k, v in e.items()}, ratio=e[0.05] / e[0.025])
        # G1 reduced set
        reg = json.loads((ROOT / 'results/EX-GC34F/spectrum.json').read_text())['directions']
        names = ['100', '010', '001', 'd01+', 'd01-', 'd02+', 'd02-', 'd12+', 'd12-', 'b1++', 'b1+-', 'b1-+', 'b1--', 'generic']
        KV = G.KernelG(st, kin); rg = np.random.default_rng(SEED + 7); grow = []; nunc = 0
        for nm in names:
            for kh in (0.05, 0.2, 1.0, 3.0):
                q = kh * np.array(reg[nm]); a = F.analyse_posthoc(KV.pencil(q), q, rg)
                nunc += sum(not rr['confirmed'] for rr in a['roots'])
                grow += [(nm, kh, rr['w_over_k']) for rr in a['roots'] if rr['confirmed'] and rr['cls'] == 'growing']
        out.update(G1=len(grow) == 0, n_growing=len(grow), growing=grow[:10], n_unconfirmed=nunc)
        out['verdict'] = ('C-REMOVABLE' if (e[0.05] <= 1e-2 and out['ratio'] >= 3 and out['G1'])
                          else 'C-FIXED' if lin_max >= 0.1 else 'C-PARTIAL')
        res['R%d' % R] = out
        for o, Ko in coef.items():
            save['R%d|%d,%d,%d' % ((R,) + o)] = Ko
        print('R', R, {k: v for k, v in out.items() if k != 'growing'}, 'cpu %.0f' % (time.process_time() - t0), flush=True)
    res['cpu_s'] = time.process_time() - t0
    (OUT / 'fit.json').write_text(json.dumps(res, **DUMP))
    np.savez_compressed(OUT / 'kernel_C.npz', **save)


if __name__ == '__main__':
    dict(audit=audit, bench=bench, fit=fit)[sys.argv[1]]()
