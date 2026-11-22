"""EX-GC34A (exploratory): one shared link structure for every field -> one light cone?

Every field moves only through the same lattice structure: a site lapse/volume part (Ms = sqrt(gamma)/N,
Vs = N sqrt(gamma)) and nine link weights per site (3 axis + 6 face-diagonal directions), whose second moment is the
metric tensor A = N sqrt(gamma) gamma^{-1}:  sum_d w_d d d^T = A.  No field has its own speed coefficient.

Lagrangian of a scalar f (h^3 per site):
    L = sum_s h^3 [ Ms f_t^2/2 - Vs mu^2 f^2/2 ] - sum_s sum_d h^3 w_d(s+d/2) (f(s+d) - f(s))^2 / (2 h^2)
Fields: F1 massless scalar, F2 massive scalar, F3 two-component field with field-space metric G and an internal
mass term, F4 O(4) sigma model (nonlinear), F5 control with its own (scaled) link weights.

  python src/exgc34a.py smoke                 # flat-background method check (allowed before registration)
  python src/exgc34a.py m1  [--only TAG]       # uniform backgrounds, packet velocities
  python src/exgc34a.py m2  [--only TAG]       # bump backgrounds, deflection
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import argparse, json, platform, sys, time, hashlib
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC34A'
RAW = ROOT / 'local_raw/EX-GC34A'

# link directions in units of h: 3 axes, 6 face diagonals
DIRS = np.array([(1, 0, 0), (0, 1, 0), (0, 0, 1),
                 (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)])
PAIR = {(0, 1): (3, 4), (0, 2): (5, 6), (1, 2): (7, 8)}   # (i,j) -> (d for e_i+e_j, d for e_i-e_j)


def split_weights(A, beta_div=12.0):
    """Nine weights with second moment A (A[...,i,j] arrays). Null-space choice beta_ij = (A_ii + A_jj)/beta_div.
    beta_div = 12 gives the fourth-moment-isotropic stencil (axis 1/3, diagonal 1/6) for A = I."""
    w = [None] * 9
    beta = {}
    for (i, j), (dp, dm) in PAIR.items():
        b = (A[..., i, i] + A[..., j, j]) / beta_div
        beta[(i, j)] = b
        w[dp] = b + A[..., i, j] / 2
        w[dm] = b - A[..., i, j] / 2
    for i in range(3):
        s = sum(2 * beta[p] for p in PAIR if i in p)
        w[i] = A[..., i, i] - s
    w = np.stack(w)
    if np.min(w) <= 0:
        raise ValueError('non-positive link weight; stencil unstable')
    return w


class Grid:
    def __init__(self, n, lo, h):
        self.n = tuple(n); self.h = h; self.lo = np.array(lo, float)
        self.ax = [lo[a] + h * np.arange(n[a]) for a in range(3)]
        self.X = np.meshgrid(*self.ax, indexing='ij')

    def pos(self, shift=(0, 0, 0)):
        return [self.X[a] + self.h * shift[a] for a in range(3)]


class Background:
    """Builds Ms, Vs (site) and w[d] (link midpoints) from metric functions.
    kind: 'flat' | 'aniso' (uniform gamma_ij) | 'bump' (N=1+Phi or exp(Phi); gamma_ij = exp(2B) delta)."""
    def __init__(self, G, kind, beta_div=12.0, gam=None, A=0.05, a=2.0, gppn=1.0, nordstrom=False, scale=1.0):
        self.G = G; self.kind = kind; self.beta_div = beta_div
        self.gam = None if gam is None else np.array(gam, float)
        self.Aamp = A; self.a = a; self.gppn = gppn; self.nord = nordstrom; self.scale = scale
        self.Ms = self.site(G.pos())[0]; self.Vs = self.site(G.pos())[1]
        self.w = np.stack([self.link(d) for d in range(9)])

    # metric pieces at arbitrary positions (arrays)
    def lapse_B(self, P):
        if self.kind != 'bump':
            z = np.zeros_like(P[0]); return 1.0 + z, z
        r = np.sqrt(P[0] ** 2 + P[1] ** 2 + P[2] ** 2)
        Phi = -self.Aamp / np.sqrt(r ** 2 + self.a ** 2)
        if self.nord:
            return np.exp(Phi), Phi              # conformally flat: N = e^Phi, gamma = e^{2 Phi} delta
        return 1.0 + Phi, -self.gppn * Phi       # weak-field PPN form, gamma_ij = exp(-2 gppn Phi) delta

    def site(self, P):
        N, B = self.lapse_B(P)
        if self.kind == 'aniso':
            sg = np.sqrt(np.linalg.det(self.gam))
            return sg / N, N * sg
        return np.exp(3 * B) / N, N * np.exp(3 * B)

    def Amat(self, P):
        N, B = self.lapse_B(P)
        shp = P[0].shape
        if self.kind == 'aniso':
            sg = np.sqrt(np.linalg.det(self.gam)); gi = np.linalg.inv(self.gam)
            A = np.broadcast_to(gi * sg, shp + (3, 3)) * N[..., None, None]
        else:
            A = np.zeros(shp + (3, 3))
            f = N * np.exp(B)                    # N sqrt(gamma) gamma^{-1} = N e^{3B} e^{-2B}
            for i in range(3):
                A[..., i, i] = f
        return A * self.scale

    def link(self, d):
        dv = DIRS[d]
        P = self.G.pos(0.5 * dv)
        return split_weights(self.Amat(P), self.beta_div)[d]


def shift(f, d, sgn):
    """f(s + sgn*d), periodic."""
    dv = DIRS[d] * sgn
    return np.roll(f, shift=tuple(-dv), axis=(-3, -2, -1))


def lap_w(f, w, h):
    """sum_d [w_d(s) (f(s+d)-f(s)) + w_d(s-d) (f(s-d)-f(s))] / h^2 ; w_d(s) is the link s -> s+d."""
    out = np.zeros_like(f)
    for d in range(9):
        fp = shift(f, d, +1); fm = shift(f, d, -1); wm = shift(w[d], d, -1)
        out += w[d] * (fp - f) + wm * (fm - f)
    return out / h ** 2


def link_energy_density(f, w, h):
    """Per-site energy density of the gradient term, each link split half to each end. f may have leading comps."""
    e = np.zeros(f.shape[-3:])
    for d in range(9):
        df2 = (shift(f, d, +1) - f) ** 2
        if df2.ndim == 4:
            df2 = df2.sum(0)
        L = 0.5 * w[d] * df2 / h ** 2
        e += 0.5 * L + 0.5 * shift(L, d, -1)
    return e


# ---------------- dispersion and packets ----------------
def omega_flat(kvec, wloc, Ms, Vs, mu, h):
    """Lattice dispersion for constant coefficients; kvec shape (...,3); wloc 9 weights."""
    q = np.zeros(kvec.shape[:-1])
    for d in range(9):
        q += wloc[d] * 4 * np.sin(h * (kvec @ DIRS[d]) / 2) ** 2 / h ** 2
    return np.sqrt((q + Vs * mu ** 2) / Ms)


def forward_packet(G, x0, k0, sigma, amp, wloc, Ms, Vs, mu):
    """Complex packet z = amp*env*exp(i k.x) and its forward time derivative using the local constant dispersion."""
    X = G.X
    r2 = sum((X[a] - x0[a]) ** 2 for a in range(3))
    z = amp * np.exp(-r2 / (2 * sigma ** 2)) * np.exp(1j * sum(k0[a] * (X[a] - x0[a]) for a in range(3)))
    kk = np.meshgrid(*[2 * np.pi * np.fft.fftfreq(G.n[a], d=G.h) for a in range(3)], indexing='ij')
    K = np.stack(kk, -1)
    om = omega_flat(K, wloc, Ms, Vs, mu, G.h)
    zt = np.fft.ifftn(-1j * om * np.fft.fftn(z))
    return z, zt


# ---------------- fields ----------------
class Scalar:
    def __init__(self, bg, mu, w=None):
        self.bg = bg; self.mu = mu; self.w = bg.w if w is None else w

    def accel(self, f):
        bg = self.bg
        return (lap_w(f, self.w, bg.G.h) - bg.Vs * self.mu ** 2 * f) / bg.Ms

    def edens(self, f, v):
        bg = self.bg
        return 0.5 * bg.Ms * v ** 2 + 0.5 * bg.Vs * self.mu ** 2 * f ** 2 + link_energy_density(f, self.w, bg.G.h)


class TwoComp:
    """u = (u1,u2): L = 1/2 u_t.G u_t Ms - 1/2 grad terms with the same G - Vs eps (u1-u2)^2/2."""
    def __init__(self, bg, g, eps):
        self.bg = bg; self.g = g; self.eps = eps
        self.Gm = np.array([[1, g], [g, 1]]); self.Gi = np.linalg.inv(self.Gm)

    def accel(self, u):
        bg = self.bg
        L = np.stack([lap_w(u[0], bg.w, bg.G.h), lap_w(u[1], bg.w, bg.G.h)])
        d = u[0] - u[1]
        pot = bg.Vs * self.eps * np.stack([d, -d])
        F = np.einsum('ab,b...->a...', self.Gm, L) - pot
        return np.einsum('ab,b...->a...', self.Gi, F) / bg.Ms

    def modes(self, u, v):
        return (u[0] + u[1]) / 2, (u[0] - u[1]) / 2, (v[0] + v[1]) / 2, (v[0] - v[1]) / 2


class Sigma:
    """O(4) sigma model, n.n = 1."""
    def __init__(self, bg):
        self.bg = bg

    def accel(self, n, v):
        bg = self.bg
        F = np.stack([lap_w(n[c], bg.w, bg.G.h) for c in range(4)])
        lam = (n * F).sum(0) + bg.Ms * (v * v).sum(0)
        return (F - lam * n) / bg.Ms

    def edens(self, n, v):
        bg = self.bg
        return 0.5 * bg.Ms * (v * v).sum(0) + link_energy_density(n, bg.w, bg.G.h)


def centroid(e, G):
    s = e.sum()
    return [float((e * G.X[a]).sum() / s) for a in range(3)], float(s * G.h ** 3)


def run(kind_field, bg, x0, k0, sigma, T, dt, every, amp=1e-3, mu=0.0, two=(0.4, 0.6), w_own=None,
        save_final=None):
    G = bg.G
    i0 = tuple(int(round((x0[a] - G.lo[a]) / G.h)) for a in range(3))
    wloc = bg.w[(slice(None),) + i0]; Ms0 = bg.Ms[i0]; Vs0 = bg.Vs[i0]
    hist = {'t': []}
    t0 = time.time()
    if kind_field in ('scalar', 'control'):
        F = Scalar(bg, mu, w_own)
        wl = wloc if w_own is None else w_own[(slice(None),) + i0]
        z, zt = forward_packet(G, x0, k0, sigma, amp, wl, Ms0, Vs0, mu)
        f = z.real.copy(); v = zt.real.copy()
        acc = F.accel
        def meas():
            return {'c': centroid(F.edens(f, v), G)}
    elif kind_field == 'two':
        g, eps = two
        F = TwoComp(bg, g, eps)
        mua = np.sqrt(2 * eps / (1 - g))
        zs, zts = forward_packet(G, x0, k0, sigma, amp, wloc, Ms0, Vs0, 0.0)
        za, zta = forward_packet(G, x0, k0, sigma, amp, wloc, Ms0, Vs0, mua)
        f = np.stack([zs.real + za.real, zs.real - za.real]); v = np.stack([zts.real + zta.real, zts.real - zta.real])
        acc = F.accel
        Ss = Scalar(bg, 0.0); Sa = Scalar(bg, mua)
        def meas():
            s, a, vs, va = F.modes(f, v)
            return {'c_sym': centroid(Ss.edens(s, vs), G), 'c_anti': centroid(Sa.edens(a, va), G)}
    elif kind_field == 'sigma':
        F = Sigma(bg)
        z, zt = forward_packet(G, x0, k0, sigma, 1.0, wloc, Ms0, Vs0, 0.0)
        e = np.array([1, 0.6j, 0.4]); e = e / np.sqrt((abs(e) ** 2).sum())
        pi = np.stack([(amp * z * e[c]).real for c in range(3)])
        pv = np.stack([(amp * zt * e[c]).real for c in range(3)])
        p0 = np.sqrt(1 - (pi ** 2).sum(0))
        f = np.concatenate([p0[None], pi]); v = np.concatenate([(-(pi * pv).sum(0) / p0)[None], pv])
        acc = lambda n: F.accel(n, v)
        def meas():
            return {'c': centroid(F.edens(f, v), G)}
    else:
        raise ValueError(kind_field)

    nstep = int(round(T / dt))
    a = acc(f)
    for it in range(nstep + 1):
        if it % every == 0:
            m = meas(); hist['t'].append(it * dt)
            for kname, val in m.items():
                hist.setdefault(kname, []).append(val[0]); hist.setdefault(kname + '_E', []).append(val[1])
        if it == nstep:
            break
        v += 0.5 * dt * a
        f += dt * v
        if kind_field == 'sigma':
            f /= np.sqrt((f * f).sum(0))
            v -= (v * f).sum(0) * f
        a = acc(f)
        v += 0.5 * dt * a
        if kind_field == 'sigma':
            v -= (v * f).sum(0) * f
    if save_final is not None:
        RAW.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(save_final, f=f, v=v)
    hist['cpu_s'] = time.process_time(); hist['wall_s'] = time.time() - t0
    return hist


def fit_velocity(t, c, t_from):
    t = np.asarray(t); c = np.asarray(c); m = t >= t_from
    return [float(np.polyfit(t[m], c[m, a], 1)[0]) for a in range(3)]


def platform_info():
    import scipy
    return {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
            'machine': platform.platform()}


def save(tag, r):
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / (tag + '.json'); q = p.with_suffix('.tmp')
    q.write_text(json.dumps(r, indent=1) + '\n'); q.replace(p)


# ---------------- smoke test (method check on a flat background) ----------------
def smoke():
    h = 0.2; n = (96, 64, 64); G = Grid(n, (-9.6, -6.4, -6.4), h)
    bg = Background(G, 'flat')
    k0 = np.array([3.0, 0, 0]); x0 = np.array([-5.0, 0, 0])
    hist = run('scalar', bg, x0, k0, 2.0, T=8.0, dt=0.25 * h, every=4)
    vm = fit_velocity(hist['t'], hist['c'], 2.0)
    print('measured v', vm, 'E drift', hist['c_E'][-1] / hist['c_E'][0] - 1, 'wall', hist['wall_s'])
    return vm


def gam_B1():
    def R(ax, t):
        c, s_ = np.cos(t), np.sin(t); M = np.eye(3); i, j = [k for k in range(3) if k != ax]
        M[i, i] = c; M[j, j] = c; M[i, j] = -s_; M[j, i] = s_; return M
    Q = R(0, np.radians(20)) @ R(2, np.radians(30))
    return Q @ np.diag([1.25, 1.0, 0.8]) @ Q.T


H = 0.2; DT = 0.25 * H; SIG = 2.0; KMAG = 3.0
M1_BG = {'B0': dict(kind='flat'), 'B1': dict(kind='aniso', beta_div=12.0), 'B1p': dict(kind='aniso', beta_div=16.0)}
M1_DIR = {'x': (1, 0, 0), 'diag': (1, 1, 1), 'gen': (0.8, 0.5, 0.33)}
M1_FIELDS = ['F1', 'F2', 'F3', 'F4', 'F5']
M2_BG = {'B2': dict(kind='bump', gppn=0.0), 'B3': dict(kind='bump', gppn=1.0), 'B4': dict(kind='bump', nordstrom=True)}
M2_FIELDS = ['F1', 'F2', 'F3', 'F4']


def make_bg(G, spec, scale=1.0):
    spec = dict(spec)
    if spec['kind'] == 'aniso':
        spec['gam'] = gam_B1()
    return Background(G, scale=scale, **spec)


def field_run(fid, bg, x0, k0, T, every, final):
    if fid == 'F1':
        return run('scalar', bg, x0, k0, SIG, T, DT, every, save_final=final)
    if fid == 'F2':
        return run('scalar', bg, x0, k0, SIG, T, DT, every, mu=1.0, save_final=final)
    if fid == 'F3':
        return run('two', bg, x0, k0, SIG, T, DT, every, two=(0.4, 0.6), save_final=final)
    if fid == 'F4':
        return run('sigma', bg, x0, k0, SIG, T, DT, every, amp=0.5, save_final=final)
    if fid == 'F5':
        bg5 = make_bg(bg.G, bg._spec, scale=1.10)
        return run('control', bg, x0, k0, SIG, T, DT, every, w_own=bg5.w, save_final=final)


def do_m1(only=None):
    n = (80, 80, 80); G = Grid(n, (-8.0, -8.0, -8.0), H)
    for bname, spec in M1_BG.items():
        bg = make_bg(G, spec); bg._spec = spec
        for dname, dv in M1_DIR.items():
            k0 = KMAG * np.array(dv, float) / np.linalg.norm(dv)
            for fid in M1_FIELDS:
                tag = f'M1_{bname}_{dname}_{fid}'
                if only and only not in tag:
                    continue
                if (OUT / (tag + '.json')).exists():
                    continue
                hist = field_run(fid, bg, np.zeros(3), k0, 5.0, 4, RAW / (tag + '.npz'))
                save(tag, dict(tag=tag, bg=bname, dir=dname, field=fid, k0=k0.tolist(), h=H, dt=DT, sigma=SIG,
                               T=5.0, hist=hist, platform=platform_info()))
                print(tag, 'wall %.1f' % hist['wall_s'], flush=True)


def do_m2(only=None):
    lo = (-24.0, -6.0, -6.0); n = (240, 100, 60); G = Grid(n, lo, H)
    for bname, spec in M2_BG.items():
        bg = make_bg(G, spec); bg._spec = spec
        for fid in M2_FIELDS:
            tag = f'M2_{bname}_{fid}'
            if only and only not in tag:
                continue
            if (OUT / (tag + '.json')).exists():
                continue
            hist = field_run(fid, bg, np.array([-18.0, 6.0, 0.0]), np.array([KMAG, 0, 0]), 36.0, 10,
                             RAW / (tag + '.npz'))
            save(tag, dict(tag=tag, bg=bname, field=fid, h=H, dt=DT, sigma=SIG, T=36.0, hist=hist,
                           platform=platform_info()))
            print(tag, 'wall %.1f' % hist['wall_s'], flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('cmd'); ap.add_argument('--only', default=None)
    A = ap.parse_args()
    if A.cmd == 'smoke':
        smoke()
    elif A.cmd == 'm1':
        do_m1(A.only)
    elif A.cmd == 'm2':
        do_m2(A.only)
