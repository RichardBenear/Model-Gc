"""EX-GC26C (exploratory): Skyrmion sourcing a dynamical Nordström scalar field, common action.

L = L_m[phi; W = exp(2 Phi)] + sum_i h^3 (Phi_t^2 - |grad_h Phi|^2) / (8 pi G_N)
=> Phi_tt = lap_h Phi - 4 pi G_N rho,  rho_i = 2 W_i (e2_i - k2_i) / (12 pi^2)
Matter: the EX-GC26E1 kernel with W and Wd = 2 Phi_t W (unchanged). Phi lives on a larger box (L_PHI) with Phi = 0 on
its wall; the matter grid (L = 7) sits at its centre with the same h. Massless Skyrme reference model, not Model Gc.

  python src/exgc26c.py validate
  python src/exgc26c.py run --B 1|2 [--GN 0.01] [--T 12]
"""
import os, sys
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import argparse, hashlib, json, platform, time
from pathlib import Path
import numpy as np
import scipy
import exgc26e1 as e1

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC26C'
RAW = ROOT / 'local_raw/EX-GC26C'
NORM = e1.NORM
L_PHI = 17.0


def save(tag, r):
    OUT.mkdir(parents=True, exist_ok=True); p = OUT / (tag + '.json'); q = p.with_suffix('.tmp')
    q.write_text(json.dumps(r, indent=1, allow_nan=False) + '\n'); q.replace(p)


class Field:
    def __init__(self, G, Lphi):
        self.h = G.h; self.n = int(round(2 * Lphi / G.h)) + 1; self.off = (self.n - G.n) // 2
        assert (self.n - G.n) % 2 == 0
        ax = (np.arange(self.n) - (self.n - 1) / 2) * G.h
        assert np.allclose(ax[self.off:self.off + G.n], G.ax)
        self.ax = ax; X = np.meshgrid(ax, ax, ax, indexing='ij'); self.r = np.sqrt(X[0] ** 2 + X[1] ** 2 + X[2] ** 2)
        self.sub = (slice(self.off, self.off + G.n),) * 3

    def lap(self, F):
        out = np.zeros_like(F); c = F[1:-1, 1:-1, 1:-1]
        out[1:-1, 1:-1, 1:-1] = (F[2:, 1:-1, 1:-1] + F[:-2, 1:-1, 1:-1] + F[1:-1, 2:, 1:-1] + F[1:-1, :-2, 1:-1]
                                 + F[1:-1, 1:-1, 2:] + F[1:-1, 1:-1, :-2] - 6 * c) / self.h ** 2
        return out

    def energy(self, F, Pi, GN):
        g2 = sum(np.sum(np.diff(F, axis=a) ** 2) for a in range(3)) / self.h ** 2
        return float(self.h ** 3 * (np.sum(Pi ** 2) + g2) / (8 * np.pi * GN))


def e2_density(G, p):
    """Exact site gradient density e2_i = sum_o w m_o(i) tr_o(i) of the eight-orientation stencil (as in the kernel)."""
    n = G.n; q = {}; val = {}
    for a in range(3):
        for s in (1, -1):
            d = (np.roll(p, -s, axis=a + 1) - p) / G.h; m = np.ones((n, n, n), bool)
            idx = [slice(None)] * 3; idx[a] = slice(n - 1, n) if s == 1 else slice(0, 1); m[tuple(idx)] = False
            q[a, s] = np.sum(d * d, axis=0) * m; val[a, s] = m.astype(float)
    e2 = np.zeros((n, n, n))
    for a in range(3):
        other = np.ones((n, n, n))
        for b in range(3):
            if b != a: other = other * (val[b, 1] + val[b, -1])
        e2 += (q[a, 1] * val[a, 1] + q[a, -1] * val[a, -1]) * other
    return e2 / 8


def rho(G, p, v, W):
    return 2 * W * (e2_density(G, p) - np.sum(v * v, axis=0)) / NORM


def deriv(G, Fd, s, GN):
    p, v, F, Pi = s
    Fs = F[Fd.sub]; W = np.exp(2 * Fs); Wd = 2 * Pi[Fd.sub] * W
    a = e1.accel(G, p, v, W, Wd, 0., 2)
    PiD = Fd.lap(F); PiD[Fd.sub] -= 4 * np.pi * GN * rho(G, p, v, W)
    PiD[0, :, :] = PiD[-1, :, :] = PiD[:, 0, :] = PiD[:, -1, :] = PiD[:, :, 0] = PiD[:, :, -1] = 0.
    return v, a, Pi, PiD


def step(G, Fd, s, dt, GN):
    p, v, F, Pi = s
    def add(s0, k, c):
        q = G.fix(s0[0] + c * k[0]); return (q, e1.ref.tangent(q, s0[1] + c * k[1]) * G.inner, s0[2] + c * k[2], s0[3] + c * k[3])
    k1 = deriv(G, Fd, s, GN); k2 = deriv(G, Fd, add(s, k1, dt / 2), GN); k3 = deriv(G, Fd, add(s, k2, dt / 2), GN); k4 = deriv(G, Fd, add(s, k3, dt), GN)
    k = tuple((a + 2 * b + 2 * c + d) / 6 for a, b, c, d in zip(k1, k2, k3, k4))
    return add(s, k, dt)


def hm(G, p, v, F, Fd):
    W = np.exp(2 * F[Fd.sub]); Hd = e1.fields(G, p, v, W)[0]; fac = G.h ** 3 / NORM
    return float(Hd.sum() * fac), Hd


def shells(Fd, F, rmax=16.0, dr=0.2):
    rb = np.arange(0, rmax + dr, dr); idx = np.digitize(Fd.r.ravel(), rb) - 1; f = F.ravel(); out = []
    for k in range(len(rb) - 1):
        sel = idx == k
        if sel.sum() < 6: continue
        x = f[sel]; rr = Fd.r.ravel()[sel]
        out.append(dict(r=float(rr.mean()), mean=float(x.mean()), std=float(x.std()), rphi=float((x * rr).mean()), n=int(sel.sum())))
    return out


def validate():
    t0 = time.process_time(); e1.configure(0.2); e1.load()
    G = e1.ref.Grid(.25, 1.5); Fd = Field(G, 3.0); rng = np.random.default_rng(263); GN = 0.05
    p = G.fix(G.hedgehog(e1.ref.PROF) + .05 * rng.normal(size=(4,) + G.r.shape)); v = e1.ref.tangent(p, rng.normal(scale=.05, size=p.shape)) * G.inner
    F = .03 * rng.normal(size=Fd.r.shape); Pi = .03 * rng.normal(size=Fd.r.shape)
    for A in (F, Pi):
        A[0, :, :] = A[-1, :, :] = A[:, 0, :] = A[:, -1, :] = A[:, :, 0] = A[:, :, -1] = 0.
    s = (p, v, F, Pi)
    def total(s):
        return hm(G, s[0], s[1], s[2], Fd)[0] + Fd.energy(s[2], s[3], GN)
    k = deriv(G, Fd, s, GN); ep = 1e-6
    def push(c):
        q = G.fix(s[0] + c * k[0]); return (q, e1.ref.tangent(q, s[1] + c * k[1]) * G.inner, s[2] + c * k[2], s[3] + c * k[3])
    rate = (total(push(ep)) - total(push(-ep))) / (2 * ep)
    # e2_density vs the orientation sum used in fields()
    _, _, _, e2f = e1.fields(G, p, v, np.ones(G.r.shape)); e2err = float(abs(e2_density(G, p) - e2f).max())
    # discrete Green's function: interior point source, compare r*Phi with the continuum far from walls
    Fd2 = Field(G, 6.0); src = np.zeros(Fd2.r.shape); c = Fd2.n // 2; src[c, c, c] = 1 / Fd2.h ** 3
    from scipy.sparse.linalg import cg
    from scipy.sparse import diags, kron, identity
    m = Fd2.n - 2; D = diags([1, -2, 1], [-1, 0, 1], shape=(m, m)) / Fd2.h ** 2; I = identity(m)
    Lap = kron(kron(D, I), I) + kron(kron(I, D), I) + kron(kron(I, I), D)
    sol, info = cg(-Lap.tocsr(), -(-4 * np.pi * src[1:-1, 1:-1, 1:-1]).ravel(), rtol=1e-10, maxiter=4000)
    Phi = np.zeros(Fd2.r.shape); Phi[1:-1, 1:-1, 1:-1] = -sol.reshape(m, m, m)   # lap Phi = 4 pi delta -> Phi = -1/r + wall terms
    sh = shells(Fd2, Phi, rmax=3.0)
    out = dict(energy_rate=float(rate), energy_scale=float(total(s)), e2_density_error=e2err, green_cg_info=int(info),
               green_rphi=[(round(x['r'], 3), round(x['rphi'], 5)) for x in sh if x['r'] > .5],
               platform=f'Python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__} {platform.platform()}',
               cpu_s=time.process_time() - t0)
    save('validation', out); print(json.dumps(out, indent=1), flush=True)
    assert abs(rate) < 1e-7 * abs(out['energy_scale']) + 1e-9 and e2err < 1e-12


def run(B, GN, T):
    t0 = time.process_time(); e1.configure(0.2); e1.load(); G = e1.ref.Grid(e1.H, e1.L); Fd = Field(G, L_PHI)
    st = json.loads((e1.OUT / f'static_B{B}.json').read_text()); rp = e1.restpath(B)
    assert hashlib.sha256(rp.read_bytes()).hexdigest() == st['sha256']
    p = np.load(rp); v = np.zeros_like(p); F = np.zeros(Fd.r.shape); Pi = np.zeros_like(F); s = (p, v, F, Pi)
    dt = .05 * G.h; steps = round(T / dt); every = round(.5 / dt); rows = []; profiles = {}; fac = G.h ** 3 / NORM
    d0 = e1.a2.geometric(p); tag = f'B{B}_GN{GN:g}'
    for i in range(steps + 1):
        t = i * dt
        if i % every == 0 or i == steps:
            p, v, F, Pi = s; H, Hd = hm(G, p, v, F, Fd); W = np.exp(2 * F[Fd.sub])
            Mg = float(np.sum(rho(G, p, v, W)) * G.h ** 3)
            rows.append(dict(t=t, H_m=H, H_phi=Fd.energy(F, Pi, GN), M_g=Mg, phi_centre=float(F[(Fd.n // 2,) * 3]),
                             centroid=[float(np.sum(Hd * x) / Hd.sum()) for x in G.X], E2W=float(np.sum(W * e2_density(G, p)) * fac)))
            if abs(t - round(t)) < 1e-9 and round(t) % 4 == 0:
                profiles[str(round(t))] = shells(Fd, F)
            print(tag, round(t, 2), round(time.process_time() - t0), rows[-1], flush=True)
        if i < steps:
            s = step(G, Fd, s, dt, GN)
    p, v, F, Pi = s; RAW.mkdir(parents=True, exist_ok=True); end = RAW / f'{tag}_end.npz'; np.savez(end, p=p, v=v, F=F, Pi=Pi)
    r = dict(B=B, GN=GN, h=G.h, L=e1.L, L_phi=L_PHI, T=steps * dt, dt=dt, E=st['E'], E2=st['E2'], rows=rows, profiles=profiles,
             degree_initial=d0, degree_final=e1.a2.geometric(p), rest_sha256=st['sha256'], end_sha256=hashlib.sha256(end.read_bytes()).hexdigest(),
             cpu_s=time.process_time() - t0, platform=f'Python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__} {platform.platform()}')
    save(tag, r); print('DONE', tag, r['cpu_s'], flush=True)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('phase', choices=['validate', 'run']); ap.add_argument('--B', type=int, choices=[1, 2], default=1)
    ap.add_argument('--GN', type=float, default=0.01); ap.add_argument('--T', type=float, default=12.0); a = ap.parse_args()
    if a.phase == 'validate': validate()
    else: run(a.B, a.GN, a.T)


if __name__ == '__main__':
    main()
