"""EX-GC34P (exploratory): Y1 linear gravity sourced by the 34O scalar packet (no feedback on the scalar).
Plan: spec/EX-GC34P_plan.md.

  python src/exgc34p.py norm       -> results/EX-GC34P/norm.json    (P0)
  python src/exgc34p.py compat     -> results/EX-GC34P/compat.json  (P1)
  python src/exgc34p.py curv DIR   -> results/EX-GC34P/curv_DIR.json (P2; DIR = axis | generic)
  python src/exgc34p.py recompute  -> results/EX-GC34P/recompute.json
Reuses exgc34n (G_up, PM, OFFS, PAIRS), exgc34o.msources2 (34O sources), exgc34n.msources (34N null control).
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import json, platform, sys, time
from pathlib import Path
import numpy as np
import scipy.fft as sf

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import exgc34n as N
import exgc34o as O
import exgc34m as Mm

OUT = ROOT / 'results/EX-GC34P'
LBOX = 32.0; K0 = 0.6; SIG = 3.0; MASS = 1.0; KAP = N.KAP
NCUT = 25                      # |n_a| <= 25  <=>  |k_a| <= 4.91
NREF = 96                      # continuum reference grid
IH = list(range(4, 10)); IC = [0, 1, 2, 3]
WT = np.array([1, 1, 1, 2, 2, 2.0])        # tidal-norm weights (off-diagonal entries counted twice)
DUMP = N.DUMP


def kgrid(Nn):
    return 2 * np.pi * np.fft.fftfreq(Nn, d=LBOX / Nn)


def mask_modes():
    n = np.arange(-NCUT, NCUT + 1)
    I, J, K = np.meshgrid(n, n, n, indexing='ij')
    keep = ~((I == 0) & (J == 0) & (K == 0))
    return np.stack([I[keep], J[keep], K[keep]], 1)          # (M, 3) integer mode numbers


# ---------------- packet with fixed physical normalisation ----------------
def amp(KX, KY, KZ, dname):
    d = np.array(Mm.DIRS[dname]); d = d / np.linalg.norm(d); k0 = K0 * d
    x0 = np.full(3, LBOX / 2) - 6 * d
    return ((2 * np.pi / LBOX) ** 3 * (SIG ** 2 / (2 * np.pi)) ** 1.5
            * np.exp(-((KX - k0[0]) ** 2 + (KY - k0[1]) ** 2 + (KZ - k0[2]) ** 2) * SIG ** 2 / 2)
            * np.exp(-1j * (KX * x0[0] + KY * x0[1] + KZ * x0[2])))


class Packet:
    def __init__(self, Nn, dname, lattice=True):
        self.N = Nn; self.h = LBOX / Nn; kk = kgrid(Nn)
        self.KK = np.meshgrid(kk, kk, kk, indexing='ij')
        self.a = amp(*self.KK, dname)
        if lattice:
            self.w = np.sqrt(MASS ** 2 + sum(4 * np.sin(K * self.h / 2) ** 2 for K in self.KK) / self.h ** 2)
        else:
            self.w = np.sqrt(MASS ** 2 + sum(K ** 2 for K in self.KK))
        self.scale = 1.0

    def field(self, t, s=1):
        Z = self.a * np.exp(-1j * s * self.w * t) * self.N ** 3 * self.scale
        phi = sf.ifftn(Z).real; phid = sf.ifftn(-1j * s * self.w * Z).real
        return phi, phid

    def cgrads(self, t, s=1):
        Z = self.a * np.exp(-1j * s * self.w * t) * self.N ** 3 * self.scale
        return [sf.ifftn(1j * K * Z).real for K in self.KK]


def cont_T(P, t, s=1):
    """Continuum T^{mu nu} on the reference grid (spectral derivatives)."""
    phi, phid = P.field(t, s); g = P.cgrads(t, s)
    g2 = sum(x ** 2 for x in g)
    rho = 0.5 * (phid ** 2 + g2 + MASS ** 2 * phi ** 2)
    Lp = 0.5 * (-phid ** 2 + g2 + MASS ** 2 * phi ** 2)
    p = [-phid * x for x in g]
    T = [[g[a] * g[b] - (Lp if a == b else 0) for b in range(3)] for a in range(3)]
    return rho, p, T


# ---------------- source vectors on the mode set ----------------
def gather(arr_hat, Nn, modes):
    return arr_hat[modes[:, 0] % Nn, modes[:, 1] % Nn, modes[:, 2] % Nn]


def lat_f(P, t, s, modes, kvec, src=O.msources2):
    phi, phid = P.field(t, s); h = P.h; V = h ** 3
    s_n, s_b, s_d, s_o = src(phi, phid, h)
    comps = [s_n] + list(s_b) + list(s_d) + list(s_o)
    f = np.empty((len(modes), 10), complex)
    for j, (c, off) in enumerate(zip(comps, N.OFFS)):
        f[:, j] = gather(sf.fftn(c), P.N, modes) / P.N ** 3 * np.exp(-1j * kvec @ (off * h)) * 2 * KAP / V
    return f


def cont_f(P, t, s, modes):
    rho, p, T = cont_T(P, t, s)
    Tup = np.zeros((len(modes), 4, 4), complex)
    Tup[:, 0, 0] = gather(sf.fftn(rho), P.N, modes)
    for a in range(3):
        Tup[:, 0, 1 + a] = Tup[:, 1 + a, 0] = gather(sf.fftn(p[a]), P.N, modes)
        for b in range(a, 3):
            Tup[:, 1 + a, 1 + b] = Tup[:, 1 + b, 1 + a] = gather(sf.fftn(T[a][b]), P.N, modes)
    Tup /= P.N ** 3
    return KAP * Tup.reshape(-1, 16) @ N.PM


# ---------------- batched Y1 kernel ----------------
def kernel_batch(Kv, w):
    """M(K, w) for many modes, (M, 10, 10), from exgc34n.G_up's formula (vectorised)."""
    Mn = len(Kv); kl = np.concatenate([np.full((Mn, 1), -w), Kv], 1)
    ETA = N.ETA; ku = kl @ ETA; k2 = np.einsum('mi,mi->m', kl, ku)
    out = np.zeros((Mn, 10, 10))
    for j in range(10):
        hj = N.PM[:, j].reshape(4, 4); tr = np.trace(ETA @ hj)
        kh = ku @ hj.T; khk = np.einsum('mi,mi->m', kh, ku)
        G = 0.5 * (k2[:, None, None] * hj - kl[:, :, None] * kh[:, None, :] - kh[:, :, None] * kl[:, None, :]
                   + kl[:, :, None] * kl[:, None, :] * tr + ETA[None] * khk[:, None, None] - ETA[None] * (k2 * tr)[:, None, None])
        G = np.einsum('ij,mjk,kl->mil', ETA, G, ETA)
        out[:, :, j] = G.reshape(Mn, 16) @ N.PM
    return out


def mv(A, x):
    """Batched real-or-complex matrix times complex vector."""
    if np.iscomplexobj(A):
        return np.matmul(A, x[..., None])[..., 0]
    y = np.matmul(A, np.stack([x.real, x.imag], -1)); return y[..., 0] + 1j * y[..., 1]


class System:
    def __init__(self, Kv):
        M0 = kernel_batch(Kv, 0.0); Mp = kernel_batch(Kv, 1.0); Mm_ = kernel_batch(Kv, -1.0)
        C = M0; B = (Mp - Mm_) / 2; A = (Mp + Mm_) / 2 - M0
        ix = np.ix_(range(len(Kv)), IH, IH); ic = np.ix_(range(len(Kv)), IC, IH)
        self.Ahh = A[:, IH][:, :, IH]; self.Bhh = B[:, IH][:, :, IH]; self.Chh = C[:, IH][:, :, IH]
        self.Ach = A[:, IC][:, :, IH]; self.Bch = B[:, IC][:, :, IH]; self.Cch = C[:, IC][:, :, IH]
        self.Ainv = np.linalg.inv(self.Ahh)
        self.Bmax = float(np.abs(self.Bhh).max()); self.prep()

    def prep(self):
        mm = np.matmul
        self.AiC = mm(self.Ainv, self.Chh); self.AiB = 1j * mm(self.Ainv, self.Bhh)
        self.RC = self.Cch - mm(self.Ach, self.AiC); self.RB = 1j * self.Bch - mm(self.Ach, self.AiB)
        self.RA = mm(self.Ach, self.Ainv)

    def acc(self, h, hd, fh):
        out = mv(self.AiC, h) - mv(self.Ainv, fh)
        return out + mv(self.AiB, hd) if self.Bmax > 0 else out

    def cres(self, h, hd, f):
        hdd = self.acc(h, hd, f[:, IH])
        r = mv(self.RC, h) + mv(self.RA, f[:, IH]) - f[:, IC]
        r = r + mv(self.RB, hd)                   # B_ch is non-zero even though B_hh = 0
        return r, hdd

    def init(self, f):
        AiC = np.einsum('mij,mjk->mik', self.Ainv, self.Chh); AiB = np.einsum('mij,mjk->mik', self.Ainv, self.Bhh)
        C0 = np.concatenate([self.Cch - np.einsum('mij,mjk->mik', self.Ach, AiC),
                             1j * (self.Bch - np.einsum('mij,mjk->mik', self.Ach, AiB))], 2)
        rhs = f[:, IC] - np.einsum('mij,mj->mi', self.Ach, np.einsum('mij,mj->mi', self.Ainv, f[:, IH]))
        z = np.einsum('mij,mj->mi', np.linalg.pinv(C0), rhs)
        sv = np.linalg.svd(C0, compute_uv=False)
        return z[:, :6], z[:, 6:], float((sv[:, -1] / sv[:, 0]).min())


def evolve(fsrc, S, T, dt, tout):
    """RK4; fsrc(t) -> (M,10) source. Returns tidal (-1/2 hdd) at tout and max relative constraint residual."""
    f0 = fsrc(0.0); h, hd, cmin = S.init(f0)
    r0, _ = S.cres(h, hd, f0); cmax = float(np.linalg.norm(r0) / np.linalg.norm(f0[:, IC]))
    ns = int(round(T / dt)); outs = {}; t = 0.0; fa = f0; init_res = cmax
    for i in range(ns + 1):
        if any(abs(t - to) < 1e-9 for to in tout):
            _, hdd = S.cres(h, hd, fa); outs[round(t, 6)] = -0.5 * hdd
        if i == ns:
            break
        fm = fsrc(t + dt / 2); fb = fsrc(t + dt)
        k1h, k1v = hd, S.acc(h, hd, fa[:, IH])
        k2h, k2v = hd + dt / 2 * k1v, S.acc(h + dt / 2 * k1h, hd + dt / 2 * k1v, fm[:, IH])
        k3h, k3v = hd + dt / 2 * k2v, S.acc(h + dt / 2 * k2h, hd + dt / 2 * k2v, fm[:, IH])
        k4h, k4v = hd + dt * k3v, S.acc(h + dt * k3h, hd + dt * k3v, fb[:, IH])
        h = h + dt / 6 * (k1h + 2 * k2h + 2 * k3h + k4h); hd = hd + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
        t = (i + 1) * dt; fa = fb
        r, _ = S.cres(h, hd, fa); cmax = max(cmax, float(np.linalg.norm(r) / np.linalg.norm(fa[:, IC])))
    return outs, cmax, init_res, cmin


def lattice_K(modes, h):
    k = modes * 2 * np.pi / LBOX
    return k, 2 * np.sin(k * h / 2) / h


def setup_lattice(h, dname, modes):
    Nn = int(round(LBOX / h)); P = Packet(Nn, dname, True)
    k, Kv = lattice_K(modes, h)
    s_n = O.msources2(*P.field(0.0), h)[0]; E0 = -float(s_n.sum())
    return P, k, Kv, E0


def tidal_err(Rl, Rc):
    d = np.sum(WT * np.abs(Rl - Rc) ** 2); n = np.sum(WT * np.abs(Rc) ** 2)
    return float(np.sqrt(d / n))


def econt(P, t=0.0, s=1):
    rho, p, T = cont_T(P, t, s); return float(rho.sum() * P.h ** 3), rho


# ---------------- P0 ----------------
def norm():
    t0 = time.process_time(); res = dict(rows=[])
    for dname in Mm.DIRS:
        Pc = Packet(NREF, dname, False); Ec, rhoc = econt(Pc)
        X = np.meshgrid(*[np.arange(NREF) * Pc.h] * 3, indexing='ij')

        def width(rho, X):
            m = rho.sum(); c = [np.sum(rho * x) / m for x in X]
            return float(np.sqrt(sum(np.sum(rho * (x - ci) ** 2) for x, ci in zip(X, c)) / m)), c

        def kc(P):
            wts = np.abs(P.a) ** 2 * P.w ** 2; kmag = np.sqrt(sum(K ** 2 for K in P.KK))
            return float(np.sum(wts * kmag) / np.sum(wts))
        wc, cc = width(rhoc, X); kcc = kc(Pc)
        row = dict(dir=dname, E_cont=Ec, width_cont=wc, kc_cont=kcc, lat=[])
        for h in (1.0, 0.5, 0.25):
            P, k, Kv, El = setup_lattice(h, dname, mask_modes())
            Nn = P.N; Xl = np.meshgrid(*[np.arange(Nn) * h] * 3, indexing='ij')
            P.scale = np.sqrt(Ec / El)
            rho = -O.msources2(*P.field(0.0), h)[0] / h ** 3
            wl, cl = width(rho, Xl)
            old = Mm.packet(Nn, h, dname, 0.0)
            Eold = -float(O.msources2(old[0], old[1], h)[0].sum())
            row['lat'].append(dict(h=h, E_fixed_unscaled=El, E_old_normalisation=Eold, width=wl,
                                   d_width=abs(wl / wc - 1), kc=kc(P), d_kc=abs(kc(P) / kcc - 1)))
        for a, b in ((0, 1), (1, 2)):
            row['lat'][b]['old_E_ratio_vs_coarser'] = row['lat'][a]['E_old_normalisation'] / row['lat'][b]['E_old_normalisation']
        row['dE_fixed_05_025'] = abs(row['lat'][1]['E_fixed_unscaled'] / row['lat'][2]['E_fixed_unscaled'] - 1)
        res['rows'].append(row); print(json.dumps(row, **DUMP), flush=True)
    R = res['rows']
    res['P0'] = bool(all(r['dE_fixed_05_025'] <= 1e-2 for r in R) and
                     all(l['d_width'] <= 1e-2 and l['d_kc'] <= 1e-2 for r in R for l in r['lat'] if l['h'] < 1))
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True); (OUT / 'norm.json').write_text(json.dumps(res, **DUMP))
    print('P0', res['P0'], 'cpu', res['cpu_s'])


def lattice_run(h, dname, s, dt, T, tout, src=O.msources2, scale=None):
    modes = mask_modes(); P, k, Kv, El = setup_lattice(h, dname, modes)
    P.scale = scale if scale is not None else 1.0
    S = System(Kv)
    out = evolve(lambda t: lat_f(P, t, s, modes, k, src), S, T, dt, tout)
    return out, S


def escale(dname, h):
    Pc = Packet(NREF, dname, False); Ec, _ = econt(Pc)
    _, _, _, El = setup_lattice(h, dname, mask_modes())
    return float(np.sqrt(Ec / El))


# ---------------- P1 ----------------
def compat():
    t0 = time.process_time(); res = dict(rows=[])
    for dname in Mm.DIRS:
        for s in (1, -1):
            for tag, src in (('34O', O.msources2), ('34N', N.msources)):
                for dt in (0.05, 0.025):
                    (outs, cmax, cinit, cmin), S = lattice_run(0.5, dname, s, dt, 8.0, [], src)
                    row = dict(dir=dname, prep='A' if s == 1 else 'B', src=tag, h=0.5, dt=dt, cres=cmax, cinit=cinit,
                               C0_min_rel_sv=cmin, Bhh_max=S.Bmax)
                    res['rows'].append(row); print(row, 'cpu %.0f' % (time.process_time() - t0), flush=True)
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True); (OUT / 'compat.json').write_text(json.dumps(res, **DUMP))


def score():
    """P1 verdict: h = 0.5 rows from compat.json; h = 0.25 (and reported h = 1) residuals from the curv runs."""
    res = json.loads((OUT / 'compat.json').read_text()); t0 = time.process_time()
    cv = {d: {(r['h'], r['prep']): r['cres'] for r in json.loads((OUT / f'curv_{d}.json').read_text())['lat_runs']} for d in Mm.DIRS}
    res['from_curv'] = {d: [dict(h=k[0], prep=k[1], cres=v) for k, v in cv[d].items()] for d in cv}
    R = res['rows']
    g = lambda **kw: [r for r in R if all(r[k_] == v for k_, v in kw.items())]
    ok = True
    for dname in Mm.DIRS:
        for p in ('A', 'B'):
            o1 = g(dir=dname, prep=p, src='34O', h=0.5, dt=0.05)[0]['cres']; o2 = g(dir=dname, prep=p, src='34O', h=0.5, dt=0.025)[0]['cres']
            n1 = g(dir=dname, prep=p, src='34N', h=0.5, dt=0.05)[0]['cres']; n2 = g(dir=dname, prep=p, src='34N', h=0.5, dt=0.025)[0]['cres']
            q = cv[dname][(0.25, p)]
            ok &= (o1 <= 1e-6 and o1 / o2 >= 8 and n1 >= 1e-4 and abs(n2 / n1 - 1) <= 0.2 and q <= 1e-6)
    res['P1'] = bool(ok)
    OUT.mkdir(parents=True, exist_ok=True); (OUT / 'compat.json').write_text(json.dumps(res, **DUMP))
    print('P1', res['P1'], 'cpu', res['cpu_s'])


# ---------------- P2 ----------------
def curv(dname):
    t0 = time.process_time(); res = dict(dir=dname, rows=[]); tout = [4.0, 8.0]; modes = mask_modes()
    kc = modes * 2 * np.pi / LBOX
    Pc = Packet(NREF, dname, False); Sc = System(kc)
    ref = {}
    for s in (1, -1):
        ref[s] = evolve(lambda t: cont_f(Pc, t, s, modes), Sc, 8.0, 0.05, tout)
        print('cont', s, 'cres %.2e' % ref[s][1], 'cpu %.0f' % (time.process_time() - t0), flush=True)
    res['cont_cres'] = [ref[1][1], ref[-1][1]]
    # dropped source fraction (lattice h = 0.25 and continuum, t = 0)
    for h in (1.0, 0.5, 0.25):
        sc = escale(dname, h); lat = {}
        for s in (1, -1):
            (outs, cmax, cinit, cmin), S = lattice_run(h, dname, s, 0.05, 8.0, tout, O.msources2, sc)
            lat[s] = outs; res.setdefault('lat_runs', []).append(dict(h=h, prep='A' if s == 1 else 'B', cres=cmax, cinit=cinit, C0_min_rel_sv=cmin))
            print('lat', h, s, 'cres %.2e' % cmax, 'cpu %.0f' % (time.process_time() - t0), flush=True)
        for t in tout:
            A, B = lat[1][t], lat[-1][t]; Ac, Bc = ref[1][0][t], ref[-1][0][t]
            res['rows'].append(dict(h=h, t=t, err_A=tidal_err(A, Ac), err_B=tidal_err(B, Bc),
                                    err_AmB=tidal_err(A - B, Ac - Bc),
                                    AmB_over_A=float(np.sqrt(np.sum(WT * np.abs(Ac - Bc) ** 2) / np.sum(WT * np.abs(Ac) ** 2)))))
            print(res['rows'][-1], flush=True)
            if h == 0.5 and t == 8.0:
                res['probe'] = dict(modes=modes[[0, 1000, 60000]], A=A[[0, 1000, 60000]], Ac=Ac[[0, 1000, 60000]])
    # dropped fraction
    for h in (0.25,):
        Nn = int(round(LBOX / h)); P = Packet(Nn, dname, True); phi, phid = P.field(0.0)
        s_n = O.msources2(phi, phid, h)[0]; F = sf.fftn(s_n); tot = np.sum(np.abs(F) ** 2) - np.abs(F[0, 0, 0]) ** 2
        kept = np.sum(np.abs(gather(F, Nn, modes)) ** 2)
        res['dropped_fraction_energy_source_h025'] = float(1 - kept / tot)
    R = res['rows']; g = lambda h, t: next(r for r in R if r['h'] == h and r['t'] == t)
    ok = True
    for t in tout:
        for key in ('err_A', 'err_B', 'err_AmB'):
            e2, e1 = g(0.25, t)[key], g(0.5, t)[key]
            ok &= (e2 <= 2e-2 and e1 / e2 >= 3)
    res['P2'] = bool(ok); res['cpu_s'] = time.process_time() - t0
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    OUT.mkdir(parents=True, exist_ok=True); (OUT / f'curv_{dname}.json').write_text(json.dumps(res, **DUMP))
    print('P2', dname, res['P2'], 'cpu', res['cpu_s'])


def recompute():
    """Separate path for 3 probe modes (h = 0.5, axis, preparation A, t = 8): per-mode kernel from exgc34n.G_up (loop,
    not the batched kernel), sources by direct DFT sums (no FFT), initial data by lstsq, DOP853 integration."""
    from scipy.integrate import solve_ivp
    t0 = time.process_time(); h = 0.5; Nn = 64; sc = escale('axis', h)
    (outs, _, _, _), _ = lattice_run(h, 'axis', 1, 0.05, 8.0, [8.0], O.msources2, sc)
    allm = mask_modes(); Rm = outs[8.0]; mag = np.sqrt(np.sum(WT * np.abs(Rm) ** 2, 1))
    pick = list(np.argsort(mag)[-3:]) + [int(np.argsort(mag)[len(mag) // 2])]     # 3 largest + a median-size mode
    modes = allm[pick]; main = Rm[pick]
    P = Packet(Nn, 'axis', True); P.scale = sc
    k = modes * 2 * np.pi / LBOX; Kv = 2 * np.sin(k * h / 2) / h
    X = np.arange(Nn) * h
    def M_of(Kq, w):
        kl = np.array([-w, *Kq]); out = np.zeros((10, 10))
        for j in range(10):
            out[:, j] = N.PM.T @ N.G_up(N.PM[:, j].reshape(4, 4), kl).reshape(16)
        return out
    mats = []
    for Kq in Kv:
        C = M_of(Kq, 0.0); B = (M_of(Kq, 1.0) - M_of(Kq, -1.0)) / 2; A = (M_of(Kq, 1.0) + M_of(Kq, -1.0)) / 2 - C
        mats.append((A, B, C))
    def fvec(t):
        phi, phid = P.field(t, 1); s_ = O.msources2(phi, phid, h)
        comps = [s_[0]] + list(s_[1]) + list(s_[2]) + list(s_[3]); out = np.zeros((len(modes), 10), complex)
        for j, (c, off) in enumerate(zip(comps, N.OFFS)):
            for m, kq in enumerate(k):
                ph = [np.exp(-1j * kq[a] * (X + off[a] * h)) for a in range(3)]
                out[m, j] = np.einsum('ijk,i,j,k->', c, ph[0], ph[1], ph[2]) / Nn ** 3 * 2 * KAP / h ** 3
        return out
    def rhs(t, y):
        f = fvec(t); y = y.view(complex).reshape(len(modes), 12); dy = np.zeros_like(y)
        for m, (A, B, C) in enumerate(mats):
            hh, hd = y[m, :6], y[m, 6:]
            dy[m, :6] = hd
            dy[m, 6:] = np.linalg.solve(A[np.ix_(IH, IH)], C[np.ix_(IH, IH)] @ hh + 1j * B[np.ix_(IH, IH)] @ hd - f[m, IH])
        return dy.reshape(-1).view(float)
    f0 = fvec(0.0); y0 = np.zeros((len(modes), 12), complex)
    for m, (A, B, C) in enumerate(mats):
        Ai = np.linalg.inv(A[np.ix_(IH, IH)])
        C0 = np.hstack([C[np.ix_(IC, IH)] - A[np.ix_(IC, IH)] @ Ai @ C[np.ix_(IH, IH)],
                        1j * (B[np.ix_(IC, IH)] - A[np.ix_(IC, IH)] @ Ai @ B[np.ix_(IH, IH)])])
        rh = f0[m, IC] - A[np.ix_(IC, IH)] @ Ai @ f0[m, IH]
        y0[m] = np.linalg.lstsq(C0, rh, rcond=None)[0]
    sol = solve_ivp(rhs, (0, 8.0), y0.reshape(-1).view(float), method='DOP853', rtol=1e-10, atol=1e-20)
    y = sol.y[:, -1].view(complex).reshape(len(modes), 12); f8 = fvec(8.0); tid = []
    for m, (A, B, C) in enumerate(mats):
        hdd = np.linalg.solve(A[np.ix_(IH, IH)], C[np.ix_(IH, IH)] @ y[m, :6] + 1j * B[np.ix_(IH, IH)] @ y[m, 6:] - f8[m, IH])
        tid.append(-0.5 * hdd)
    tid = np.array(tid)
    rel = [float(np.linalg.norm(tid[m] - main[m]) / np.linalg.norm(tid[m])) for m in range(len(modes))]
    out = dict(modes=modes, mode_size_rel_to_largest=mag[pick] / mag.max(), rel_diff_vs_main=rel, nfev=sol.nfev,
               note='first attempt used 3 arbitrary probe modes stored by curv (edge modes of negligible size); replaced by '
                    'the 3 largest modes and a median one from a fresh main-path run', cpu_s=time.process_time() - t0)
    (OUT / 'recompute.json').write_text(json.dumps(out, **DUMP)); print(out)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'norm':
        norm()
    elif cmd == 'compat':
        compat()
    elif cmd == 'recompute':
        recompute()
    elif cmd == 'score':
        score()
    elif cmd == 'curv':
        curv(sys.argv[2])


def posthoc_spectrum():
    """POST HOC (labelled): how the tidal norm is spread over modes (h = 0.5, axis, A, t = 4 and 8)."""
    t0 = time.process_time(); sc = escale('axis', 0.5)
    (outs, _, _, _), _ = lattice_run(0.5, 'axis', 1, 0.05, 8.0, [4.0, 8.0], O.msources2, sc)
    m = mask_modes(); kmag = np.linalg.norm(m, axis=1) * 2 * np.pi / LBOX; res = dict(note='POST HOC', rows=[])
    for t, R in outs.items():
        w = np.sum(WT * np.abs(R) ** 2, 1); tot = w.sum()
        res['rows'].append(dict(t=t, frac_lowest_shell=float(w[np.abs(m).max(1) <= 1].sum() / tot),
                                frac_k_below_0p5=float(w[kmag < 0.5].sum() / tot), frac_k_0p5_to_2=float(w[(kmag >= 0.5) & (kmag < 2)].sum() / tot),
                                frac_k_above_2=float(w[kmag >= 2].sum() / tot)))
    res['cpu_s'] = time.process_time() - t0
    (OUT / 'posthoc_spectrum.json').write_text(json.dumps(res, **DUMP)); print(res)


if __name__ == '__main__' and sys.argv[1] == 'posthoc':
    posthoc_spectrum()
