"""EX-GC34I (exploratory): a local lattice kinetic symbol derived from an exact linearised re-slicing identity.
Plan: spec/EX-GC34I_plan.md.

  python src/exgc34i.py derive -> results/EX-GC34I/derive.json + kernels.npz   (solutions, I0)
  python src/exgc34i.py spec   -> results/EX-GC34I/candidates.json             (34G G1-G3, G5 machinery)
  python src/exgc34i.py evolve -> results/EX-GC34I/evolve.json                 (G4, only if earned)
  python src/exgc34i.py check  -> results/EX-GC34I/check.json                  (G7 real-space sum)
Identity (I): K(q) sigma_xi(q) = C_en(q) / l0, sigma_xi = -mu l0 Dhat.(e^{iq.d} - 1) g(q).
Reuses exgc34e/f/g unchanged. Convention e^{-i w t}.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import itertools, json, platform, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import exgc34d as D
import exgc34e as E
import exgc34f as F
import exgc34g as G

OUT = ROOT / 'results/EX-GC34I'
H = E.H; ET = E.ET; NV = 11
IN, IB, IE = F.IN, F.IB, F.IE
L0 = np.array([np.linalg.norm(ET[t] * H) for t in range(7)])
DH = np.array([ET[t] / np.linalg.norm(ET[t]) for t in range(7)])
DUMP = G.DUMP
SEED = 3434 + 9


def S_map(q, w):
    """7 x 11: strain rate of the type-t edge leaving the vertex."""
    S = np.zeros((7, NV), complex)
    for t in range(7):
        S[t, 4 + t] = L0[t] * (-1j * w)
        S[t, 1:4] = -L0[t] * DH[t] * (np.exp(1j * (q @ ET[t])) - 1)
    return S


def g_of(kind, q):
    return (np.exp(1j * q) - 1) / H if kind == 'f' else 1j * np.sin(q) / H


def sigma0(kind, q):
    return -L0 * ((DH * (np.exp(1j * (ET @ q)) - 1)[:, None]) @ g_of(kind, q))


def K_T(q):
    return G.kin_T(q, 1.0)[np.ix_(IE, IE)] / np.outer(L0, L0)


def offsets(R):
    return [o for o in itertools.product(range(-R, R + 1), repeat=3)]


def half(R):
    return [o for o in offsets(R) if o > (0, 0, 0)]


def basis(R):
    """Unknowns: o = 0 symmetric (28), o in positive half (49 each); K_{-o} = K_o^T. Returns list of (o, i, j, weight)."""
    B = []
    for i in range(7):
        for j in range(i, 7):
            B.append(((0, 0, 0), i, j, 1.0 if i == j else np.sqrt(2)))
    for o in half(R):
        for i in range(7):
            for j in range(7):
                B.append((o, i, j, np.sqrt(2)))
    return B


def K_from(coef, R, q):
    """coef: dict o -> 7x7 real for o in offsets(R) (full set)."""
    K = np.zeros((7, 7), complex)
    for o, Ko in coef.items():
        K += Ko * np.exp(1j * (q @ np.array(o, float)))
    return K


def coef_from_vec(c, R):
    B = basis(R); coef = {o: np.zeros((7, 7)) for o in offsets(R)}
    for v, (o, i, j, w) in zip(c, B):
        if o == (0, 0, 0):
            coef[o][i, j] += v; coef[o][j, i] = coef[o][i, j]
        else:
            coef[o][i, j] += v
            mo = tuple(-np.array(o)); coef[mo][j, i] += v
    return coef


def vec_from_coef(coef, R):
    return np.array([coef[o][i, j] for (o, i, j, w) in basis(R)])


def extract_T(R):
    """Real-space coefficients of K_T by FFT on a 4^3 grid (offsets within +-1), embedded in range R."""
    N = 4; ks = 2 * np.pi * np.arange(N) / N
    vals = np.zeros((N, N, N, 7, 7), complex)
    for a in range(N):
        for b in range(N):
            for c in range(N):
                vals[a, b, c] = K_T(np.array([ks[a], ks[b], ks[c]]))
    co = np.fft.fftn(vals, axes=(0, 1, 2)) / N ** 3        # K(q) = sum_o K_o e^{i q.o}  ->  K_o = mean K(q) e^{-iq.o}
    coef = {o: np.zeros((7, 7)) for o in offsets(R)}
    imag_max = 0.0
    for o in offsets(1):
        v = co[tuple(np.array(o) % N)]
        imag_max = max(imag_max, np.abs(v.imag).max())
        coef[o] = v.real.copy()
    return coef, imag_max


def derive_one(kind, R, st, mu_fixed=None):
    B = basis(R); n = len(B)
    N = 2 * R + 6; ks = 2 * np.pi * np.arange(N) / N
    rows = []; rhs = []
    for a, b, c in itertools.product(range(N), repeat=3):
        q = np.array([ks[a], ks[b], ks[c]])
        s0 = sigma0(kind, q)
        Ren = E.static_kernel(st, q[None])[0][IE, IN] / L0
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
    W = np.array([w for (_, _, _, w) in B])
    Aw = A / W                                               # variables c~ = W c (Frobenius metric)
    U, s, Vh = np.linalg.svd(Aw, full_matrices=False)
    tol = s[0] * 1e-11
    r = int((s > tol).sum())
    Ap = (Vh[:r].T / s[:r]) @ U[:, :r].T                     # pseudo-inverse
    cT, imagT = extract_T(R); xT = W * vec_from_coef(cT, R)
    # consistency: is b (scaled by 1/mu) in the range of A?
    bproj = U[:, :r] @ (U[:, :r].T @ b)
    consistent_rel = float(np.linalg.norm(b - bproj) / np.linalg.norm(b))
    # T's own leading-order mu (registered: fixes mu): small-q least squares of s K_T sigma0 = R/l0
    rng = np.random.default_rng(SEED)
    num = 0; den = 0
    for _ in range(50):
        q = 0.02 * rng.standard_normal(3)
        y = K_T(q) @ sigma0(kind, q); z = E.static_kernel(st, q[None])[0][IE, IN] / L0
        num += np.vdot(y, z).real; den += np.vdot(y, y).real
    mu_T = den / num if mu_fixed is None else mu_fixed
    # closest to T at fixed mu = mu_T: x = xT - Ap(A xT - b/mu_T); the joint (K, mu) optimum is reported only
    u1 = Ap @ (Aw @ xT); u2 = Ap @ b
    s_opt = float(u1 @ u2 / (u2 @ u2))
    x = xT - (u1 - u2 / mu_T)
    c = x / W
    coef = coef_from_vec(c, R)
    out = dict(kind=kind, R=R, grid=N, n_unknowns=n, rank=r, solution_space_dim=n - r, consistent_rel=consistent_rel,
               exists=bool(consistent_rel <= 1e-10), mu=float(mu_T), mu_joint_reported=float(1 / s_opt),
               KT_fft_imag_max=float(imagT),
               rel_change_vs_T=float(np.linalg.norm((x - xT)) / np.linalg.norm(xT)))
    supp = {}
    for o in offsets(R):
        d = np.abs(coef[o] - cT[o]).max()
        if d > 1e-12:
            supp[str(o)] = float(d)
    out['correction_support'] = supp
    return out, coef


def kin_from(coef, mu, kind):
    def kin(q, w):
        S = S_map(q, w); K = K_from(coef, None, q)
        return S.conj().T @ K @ S
    return kin


def I0(coef, R, kind, mu, st, rng):
    res = dict()
    iden = []; herm = []
    for _ in range(200):
        q = rng.uniform(-np.pi, np.pi, 3)
        K = K_from(coef, R, q)
        lhs = K @ (mu * sigma0(kind, q)); rhs_ = E.static_kernel(st, q[None])[0][IE, IN] / L0
        iden.append(np.linalg.norm(lhs - rhs_) / np.linalg.norm(rhs_))
        herm.append(np.abs(K - K.conj().T).max() / np.abs(K).max())
    res['identity_rel'] = float(max(iden)); res['hermitian'] = float(max(herm))
    kin = kin_from(coef, mu, kind); nul = []
    for _ in range(50):
        q = rng.uniform(-np.pi, np.pi, 3); w = rng.uniform(-3, 3) / H
        M = kin(q, w) + E.static_kernel(st, q[None])[0]
        v = np.zeros(NV, complex); v[IN] = -1j * w; v[1:4] = mu * g_of(kind, q)
        nul.append(np.linalg.norm(M @ v) / (np.linalg.norm(M, 2) * np.linalg.norm(v)))
    res['null_rel'] = float(max(nul))
    res['I0'] = bool(res['identity_rel'] <= 1e-12 and res['hermitian'] <= 1e-14 and res['null_rel'] <= 1e-12)
    return res


def derive():
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); rng = np.random.default_rng(SEED + 1)
    # internal check: S^dagger K_T S reproduces kin_T
    chk = []
    for _ in range(20):
        q = rng.uniform(-np.pi, np.pi, 3); w = rng.uniform(-3, 3) / H
        a = S_map(q, w).conj().T @ K_T(q) @ S_map(q, w); b = G.kin_T(q, w)
        chk.append(np.abs(a - b).max() / np.abs(b).max())
    res = dict(internal_KT_check=float(max(chk)), solutions={})
    cT1, im = extract_T(1); rec = []
    for _ in range(20):
        q = rng.uniform(-np.pi, np.pi, 3); rec.append(np.abs(K_from(cT1, 1, q) - K_T(q)).max() / np.abs(K_T(q)).max())
    res['KT_realspace_reconstruct'] = float(max(rec))
    save = {}
    for kind in ('f', 'c'):
        for R in (1, 2):
            name = 'I%s%d' % (kind, R)
            if R == 2 and res['solutions']['I%s1' % kind]['exists']:
                continue
            out, coef = derive_one(kind, R, st)
            out.update(I0(coef, R, kind, out['mu'], st, rng))
            res['solutions'][name] = out
            for o, Ko in coef.items():
                save['%s|%d,%d,%d' % ((name,) + o)] = Ko
            print(name, {k: v for k, v in out.items() if k != 'correction_support'}, flush=True)
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'derive.json').write_text(json.dumps(res, **DUMP))
    np.savez_compressed(OUT / 'kernels.npz', **save)
    print('cpu', res['cpu_s'])


def load_cands():
    d = json.loads((OUT / 'derive.json').read_text())['solutions']
    z = np.load(OUT / 'kernels.npz')
    cands = {}
    for name, s in d.items():
        if not s['exists'] or not s['I0']:
            continue
        coef = {}
        for key in z.files:
            nm, o = key.split('|')
            if nm == name:
                coef[tuple(int(v) for v in o.split(','))] = z[key]
        cands[name] = (coef, s)
    return cands


def make_kin(coef):
    def kin(q, w):
        S = S_map(q, w); K = K_from(coef, None, q)
        return S.conj().T @ K @ S
    return kin


def spec():
    cands = load_cands()
    G.CANDS = {n: make_kin(c) for n, (c, s) in cands.items()}
    G.OUT = OUT
    if not G.CANDS:
        print('no candidate with an exact solution; nothing to score'); return
    G.spec()


def evolve():
    cands = load_cands()
    G.CANDS = {n: make_kin(c) for n, (c, s) in cands.items()}
    G.OUT = OUT
    G.evolve()


def check():
    """G7: direct real-space sum 1/2 sum_v sum_o sigma(v)^T K_o sigma(v+o) on a periodic 8^3 lattice vs Fourier form."""
    cands = load_cands()
    N = 8; rng = np.random.default_rng(SEED + 2)
    beta = rng.standard_normal((3, N, N, N)); ldot = rng.standard_normal((7, N, N, N)) * H
    sh = lambda f, off: np.roll(f, shift=tuple(-np.array(off)), axis=(-3, -2, -1))
    sig = np.stack([L0[t] * (ldot[t] - np.einsum('a,a...->...', DH[t], sh(beta, ET[t]) - beta)) for t in range(7)])
    x = np.zeros((11, N, N, N), complex); x[1:4] = beta; x[4:] = ldot / (-1j)
    xh = np.fft.fftn(x, axes=(1, 2, 3)); kk = 2 * np.pi * np.fft.fftfreq(N)
    out = {}
    for name, (coef, s) in cands.items():
        direct = 0.5 * sum(np.einsum('i...,ij,j...->...', sig, Ko, sh(sig, o)).sum() for o, Ko in coef.items())
        kin = make_kin(coef); four = 0.0
        for i in range(N):
            for j in range(N):
                for l in range(N):
                    q = np.array([kk[i], kk[j], kk[l]]); v = xh[:, i, j, l]
                    four += 0.5 * (v.conj() @ kin(q, 1.0) @ v).real
        four /= N ** 3
        out[name] = dict(direct=float(direct), fourier=float(four), rel=float(abs(direct / four - 1)), G7=bool(abs(direct / four - 1) <= 1e-10))
        print(name, out[name], flush=True)
    (OUT / 'check.json').write_text(json.dumps(out, **DUMP))




# ---------------- POST HOC (labelled; after the registered run) ----------------
def posthoc():
    """If1 with mu = -1 (continuum normalisation; registered If1 wave speed^2 = |mu_T| = 0.9776). Same derivation and scoring."""
    global OUT
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); rng = np.random.default_rng(SEED + 5)
    P = OUT / 'posthoc'; P.mkdir(parents=True, exist_ok=True)
    out, coef = derive_one('f', 1, st, mu_fixed=-1.0)
    out.update(I0(coef, 1, 'f', -1.0, st, rng))
    out['label'] = 'POST HOC: mu fixed at -1 instead of the registered T least-squares value'
    print({k: v for k, v in out.items() if k != 'correction_support'}, flush=True)
    (P / 'derive.json').write_text(json.dumps(dict(solutions={'If1m1': out}), **DUMP))
    np.savez_compressed(P / 'kernels.npz', **{'If1m1|%d,%d,%d' % o: Ko for o, Ko in coef.items()})
    OUT = P
    spec(); check(); evolve()
    print('cpu', time.process_time() - t0)


if __name__ == '__main__' and sys.argv[1] not in ('posthoc_evolve', 'posthoc_sync'):
    dict(derive=derive, spec=spec, evolve=evolve, check=check, posthoc=posthoc)[sys.argv[1]]()


def posthoc_evolve():
    """POST HOC: the 34G G4 time-domain recipe for If1m1 although it fails G5 (E4) narrowly."""
    global OUT
    OUT = OUT / 'posthoc'
    cands = load_cands()
    G.CANDS = {n: make_kin(c) for n, (c, s) in cands.items()}
    G.OUT = OUT
    G.evolve(['If1m1'], 'evolve_posthoc.json', 'POST HOC: G4 recipe on If1m1 (mu = -1); not earned under the plan (G5 E4 1.25e-2 > 1e-2)')


if __name__ == '__main__' and sys.argv[1] == 'posthoc_evolve':
    posthoc_evolve()


def posthoc_sync(seed_cut=1e-3, outname='sync.json'):
    """POST HOC: synchronous gauge (n = 0, beta = 0), possible only because re-slicing is exact. Regular ODE
    A_ee e'' = C_ee e. (1) spectrum of A^-1 C over the 34G case set: growth, low-mode count, speeds; (2) zero-sector
    curvature vs displacement control; (3) implicit-midpoint evolution, wave + non-zero seeds, constraints and tidal field."""
    import scipy.linalg as sl
    sys.path.insert(0, str(ROOT / 'src')); import exgc34h as Hh
    t0 = time.process_time(); mp.mp.dps = 50
    st = E.tet_static(); P0 = OUT / 'posthoc'
    z = np.load(P0 / 'kernels.npz')
    coef = {tuple(int(v) for v in k.split('|')[1].split(',')): z[k] for k in z.files}
    KV = G.KernelG(st, make_kin(coef))
    cdef, alld = Hh.cases_def()
    for n in ('100', 'd01+', 'b1++'):
        cdef.append((n, float(np.pi / np.abs(G.unit(alld[n])).max())))
    cases = []
    for name, kh in cdef:
        q = kh * G.unit(alld[name]); k = kh / H; P = KV.pencil(q)
        Aee = P['A'][np.ix_(IE, IE)]; Cee = P['C'][np.ix_(IE, IE)]
        lam, V = np.linalg.eig(np.linalg.solve(Aee, Cee))        # e'' = A^-1 C e ; w^2 = -lam
        w = np.sqrt(-lam + 0j); w = np.where(w.real < 0, -w, w)
        zero = np.abs(lam) <= 1e-6 * k * k
        grow = [complex(v / k) for v, zz in zip(w, zero) if not zz and abs(v.imag) > 1e-8 * abs(v)]
        low = sorted([float(v.real / k) for v, zz in zip(w, zero) if not zz and abs(v.imag) <= 1e-8 * abs(v) and v.real * H < 0.5])
        c = dict(dir=name, kh=kh, n_zero=int(zero.sum()), growing=[[g.real, g.imag] for g in grow], low_speeds=low,
                 light_speed=float(E.light_speed(q)), all_w_over_k=sorted([[float(v.real / k), float(v.imag / k)] for v in w]))
        # zero sector: static displacement-like edge vectors; curvature vs displacement control
        Vz = V[:, zero]; rr = []
        for j in range(Vz.shape[1]):
            x = np.zeros(NV, complex); x[IE] = Vz[:, j]
            h = E.metric_from_x(x, q); zz4 = np.zeros((4, 4), complex)
            rr.append(float(np.linalg.norm(Hh.riemann_t(h, zz4, zz4, q)) / (k * k * np.linalg.norm(h))))
        c['zero_r'] = rr; c['r_disp'] = Hh.r_disp(q)
        cases.append(c)
    res = dict(label='POST HOC: synchronous gauge (n = beta = 0) for If1m1', cases=cases)
    res['n_growing'] = sum(len(c['growing']) for c in cases)
    res['low_mode_bad'] = [(c['dir'], c['kh'], c['low_speeds']) for c in cases if c['kh'] <= 0.05 and
                           (len(c['low_speeds']) != 2 or max(abs(v - c['light_speed']) for v in c['low_speeds']) > 1e-3)]
    res['max_speed_err_kh_le_005'] = max(abs(v - c['light_speed']) for c in cases if c['kh'] <= 0.05 for v in c['low_speeds'])
    res['n_zero_set'] = sorted({c['n_zero'] for c in cases})
    res['zero_r_over_rdisp_max'] = max(max(c['zero_r'] or [0]) / max(c['r_disp'], 1e-300) for c in cases if c['r_disp'] > 1e-12)
    res['zero_r_max_where_rdisp_zero'] = max([max(c['zero_r'] or [0]) for c in cases if c['r_disp'] <= 1e-12] or [0])
    # time domain
    runs = []
    for name in ('generic', 'd01+'):
        q = 0.05 * G.unit(alld[name]); k = 0.05 / H; P = KV.pencil(q)
        Aee = P['A'][np.ix_(IE, IE)]; Cee = P['C'][np.ix_(IE, IE)]; Fm = np.linalg.solve(Aee, Cee)
        L = np.block([[np.zeros((7, 7)), np.eye(7)], [Fm, np.zeros((7, 7))]])
        mu_, Vl = np.linalg.eig(L)
        nz = np.abs(mu_) > seed_cut * k
        iw = [j for j in range(14) if nz[j] and abs(abs(mu_[j]) / k - 1) < 0.01 and mu_[j].imag < 0]
        jw = iw[0]; zw = Vl[:, jw] / np.linalg.norm(Vl[:, jw])
        others = [Vl[:, j] / np.linalg.norm(Vl[:, j]) for j in range(14) if nz[j] and j != jw]
        Z = np.array(others).T
        wmax = np.abs(mu_).max(); T = 40 * 2 * np.pi / abs(mu_[jw])
        CnE = P['C'][IN, IE]; CbE = P['C'][np.ix_(IB, IE)]; BbE = P['B'][np.ix_(IB, IE)]
        out = dict(dir=name, wave_w_over_k=float(abs(mu_[jw]) / k), n_seeds=len(others))
        for tag, fac in (('dt', 0.02), ('dt2', 0.01)):
            ns = int(np.ceil(T / (fac / wmax))); dt = T / ns
            I14 = np.eye(14); Ms = np.linalg.solve(I14 - 0.5 * dt * L, I14 + 0.5 * dt * L)
            every = max(1, ns // 2000); Mev = np.linalg.matrix_power(Ms, every)
            z1 = zw.copy(); z2 = zw + 1e-8 * Z.sum(1); c0 = np.linalg.lstsq(Z, z2 - z1, rcond=None)[0]
            lap = mom = amp = grow = 0.0
            for i in range(0, ns + 1, every):
                for zz in (z1, z2):
                    e = zz[:7]; ed = zz[7:]
                    lap = max(lap, abs(CnE @ e) / (np.linalg.norm(CnE) * np.linalg.norm(e)))
                    mom = max(mom, np.linalg.norm(CbE @ e + 1j * BbE @ ed) / (np.linalg.norm(CbE, 2) * np.linalg.norm(e) + np.linalg.norm(BbE, 2) * np.linalg.norm(ed)))
                amp = max(amp, abs(np.linalg.norm(z1) - 1))
                cc = np.linalg.lstsq(Z, z2 - z1, rcond=None)[0]; grow = max(grow, float(np.max(np.abs(cc) / np.abs(c0))))
                z1 = Mev @ z1; z2 = Mev @ z2
            out[tag] = dict(dt=dt, steps=ns, sample_every=every, max_lapse=float(lap), max_mom=float(mom), wave_amp_change=float(amp), max_seed_growth=grow)
        out['G4_thresholds_met'] = all(out[t]['max_lapse'] <= 1e-10 and out[t]['max_mom'] <= 1e-10 and out[t]['wave_amp_change'] <= 1e-8
                                       and out[t]['max_seed_growth'] <= 1 + 1e-6 for t in ('dt', 'dt2'))
        runs.append(out)
    res['time_domain'] = runs
    res['cpu_s'] = time.process_time() - t0
    res['seed_cut_over_k'] = seed_cut
    (P0 / outname).write_text(json.dumps(res, **DUMP))
    print({k: v for k, v in res.items() if k not in ('cases',)})


if __name__ == '__main__' and sys.argv[1] == 'posthoc_sync':
    posthoc_sync()
