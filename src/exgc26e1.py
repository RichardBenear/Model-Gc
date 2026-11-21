"""EX-GC26E1 (exploratory): B = 1 vs B = 2 Skyrmion fall in a ramped uniform Nordstrom field.

Massless Skyrme reference model (c = 1), not Model Gc; lattice construction and coupling are ours.
Stencil, integrator and optimizer inherited from EX-GC26B2/B3 (ChatGPT). New here: the site weight
W_i = exp(2 Phi(x_i, t)) on the quadratic (sigma-model) kinetic and gradient terms (the quartic term is
conformally invariant), Phi = lambda(t) g x, lambda = 26B3 half-cosine ramp over t in [0, 1]; a B = 2 rational-map
start; baryon-density centroid. Control: 26B2 squared-deficit coupling with one fixed f for both species.

  python src/exgc26e1.py validate
  python src/exgc26e1.py static --B 1|2 [--maxiter N]
  python src/exgc26e1.py run --B 1|2 --mode energy|control|null
  add --h 0.15 for EX-GC26E2 (records in results/EX-GC26E2, raw in local_raw/EX-GC26E2)
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '1'
import argparse, ctypes, hashlib, itertools, json, platform, subprocess, sys, time, types
try:
    import resource
except ImportError:      # Windows (PI's PC): no CPU rlimit; stub so the inherited exgc26b2 import works
    resource = None
    sys.modules['resource'] = types.ModuleType('resource')
from pathlib import Path
import numpy as np
import scipy
import exgc26a_skyrmion as ref
import exgc26a2_qualify as a2
import exgc26b_mass as bm
import exgc26b2 as b2

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC26E1'
RAW = ROOT / 'local_raw/EX-GC26E1'
NORM = ref.NORM
H, L = 0.2, 7.0
GRAV = 0.0125          # field gradient g: nominal F/E = 2g = 0.025, as in 26B3
T_END = 5.0
LIB = None
PTR = a2.PTR


def configure(h):
    """h = 0.2 -> EX-GC26E1 (default); h = 0.15 -> EX-GC26E2. Same box, field, ramp, T and dt/h."""
    global H, OUT, RAW
    H = h; stage = 'EX-GC26E1' if h == 0.2 else 'EX-GC26E2'
    OUT = ROOT / 'results' / stage; RAW = ROOT / 'local_raw' / stage


def save(tag, r):
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / (tag + '.json'); q = p.with_suffix('.tmp')
    q.write_text(json.dumps(r, indent=1, allow_nan=False) + '\n'); q.replace(p)


def load():
    global LIB
    RAW.mkdir(parents=True, exist_ok=True)
    so = RAW / 'stencilW.so'; src = Path(__file__).with_name('exgc26e1_stencil.cpp')
    if not so.exists() or so.stat().st_mtime < src.stat().st_mtime:
        tmp = so.with_name(f'stencilW-{os.getpid()}.so')
        subprocess.run(['c++', '-O3', '-march=native', '-shared', '-fPIC', '-std=c++17', str(src), '-o', str(tmp)], check=True)
        tmp.replace(so)
    LIB = ctypes.CDLL(str(so))
    LIB.evaluateW.argtypes = [ctypes.c_int, ctypes.c_double, PTR, PTR, PTR, PTR]
    LIB.accelerationW.argtypes = [ctypes.c_int, ctypes.c_double, PTR, PTR, PTR, PTR, ctypes.c_double, ctypes.c_int, PTR]
    LIB.evaluateW.restype = LIB.accelerationW.restype = None
    b2.load()   # 26B2 kernel, for the W = 1 equality check and the inherited optimizer path


def cvec(a):
    return np.ascontiguousarray(np.moveaxis(a, 0, -1))


def evaluate(G, p, W=None):
    W = np.ones(p.shape[1:]) if W is None else W
    x = cvec(p); y = np.empty_like(x); e = np.zeros(2)
    LIB.evaluateW(G.n, G.h, x, np.ascontiguousarray(W), y, e)
    return e, np.moveaxis(y, -1, 0)


def accel(G, p, v, W, Wd, f=0., power=2):
    x = cvec(p); y = cvec(v); z = np.empty_like(x)
    LIB.accelerationW(G.n, G.h, x, y, np.ascontiguousarray(W), np.ascontiguousarray(Wd), f, power, z)
    return np.moveaxis(z, -1, 0)


def ramp(t):
    return .5 * (1 - np.cos(np.pi * np.clip(t, 0, 1)))


def dramp(t):
    return .5 * np.pi * np.sin(np.pi * t) if 0 < t < 1 else 0.


def weights(G, t, g):
    if g == 0:
        one = np.ones(G.r.shape); return one, np.zeros_like(one)
    W = np.exp(2 * g * ramp(t) * G.X[0])
    return W, 2 * g * dramp(t) * G.X[0] * W


# ---------- independent NumPy reference (validation only) ----------
def numpy_ref(G, p, v, W, Wd, f, power=2):
    E = np.zeros(2); grad = np.zeros_like(p); R = np.zeros_like(p); M = np.zeros((4, 4) + p.shape[1:]); vv = bm.dot(v, v)
    for s in itertools.product((-1, 1), repeat=3):
        m, d, D = b2.orient(G, p, s); _, dv, _ = b2.orient(G, v, s)
        tr = sum(D[i, i] for i in range(3)); e4 = .5 * (tr * tr - np.sum(D * D, axis=(0, 1)))
        E += np.array([(W * tr * m).sum(), (e4 * m).sum()]) / 8
        vd = [bm.dot(v, x) for x in d]; sd = sum(bm.dot(x, y) for x, y in zip(d, dv))
        M += ((W + tr) * np.eye(4)[:, :, None, None, None] - sum(x[:, None] * x[None, :] for x in d)) / 8   # d is masked, so W I survives everywhere as in the kernel
        R -= (2 * sd * v - sum(dv[i] * vd[i] + d[i] * bm.dot(dv[i], v) for i in range(3))) * m / 8
        for i in range(3):
            q = (W + tr) * d[i] - sum(D[i, j] * d[j] for j in range(3))
            grad += b2.adj(G, 2 * q, i, s[i], m) / 8
            R += b2.adj(G, vv * d[i] - vd[i] * v - q, i, s[i], m) / 8
    R -= Wd * v
    R[0] -= .5 * f * G.X[0] * power * (1 - p[0]) ** (power - 1)
    A = np.moveaxis(M, (0, 1), (-2, -1)); rhs = np.stack([np.moveaxis(R, 0, -1), np.moveaxis(p, 0, -1)], -1)
    sol = np.linalg.solve(A, rhs); a = np.moveaxis(sol[..., 0], -1, 0); c = np.moveaxis(sol[..., 1], -1, 0)
    return E * G.h ** 3 / NORM, grad * G.h ** 3 / NORM, (a + c * (-vv - bm.dot(p, a)) / bm.dot(p, c)) * G.inner


def fields(G, p, v, W):
    """Coupled energy density, canonical momentum density, baryon density, gradient-energy density e2."""
    Hd = np.zeros(p.shape[1:]); P = np.zeros((3,) + Hd.shape); Bd = np.zeros_like(Hd); e2 = np.zeros_like(Hd); vv = bm.dot(v, v)
    for s in itertools.product((-1, 1), repeat=3):
        m, d, D = b2.orient(G, p, s); tr = sum(D[i, i] for i in range(3))
        vd = np.array([bm.dot(v, x) for x in d]); K = vv * (W + tr) - np.sum(vd ** 2, axis=0)
        e = W * tr + .5 * (tr * tr - np.sum(D ** 2, axis=(0, 1)))
        Hd += (K + e) * m / 8; e2 += tr * m / 8
        for i in range(3):
            P[i] -= 2 * ((W + tr) * vd[i] - sum(D[i, j] * vd[j] for j in range(3))) * m / 8
        Bd += np.linalg.det(np.moveaxis(np.stack([p, d[0], d[1], d[2]], 0), (0, 1), (-1, -2))) * m / 8
    return Hd, P, Bd, e2


def measure(G, p, v, t, g, f):
    W, _ = weights(G, t, g); Hd, P, Bd, e2 = fields(G, p, v, W); fac = G.h ** 3 / NORM
    lam = ramp(t); dq = (1 - p[0]) ** 2
    V = -f * lam * float(np.sum(G.X[0] * dq) * fac)
    return dict(t=t, H=float(Hd.sum() * fac), V=V, P=(P.sum((1, 2, 3)) * fac).tolist(),
                centroid=[float(np.sum(Hd * x) / Hd.sum()) for x in G.X],
                bcentroid=[float(np.sum(Bd * x) / Bd.sum()) for x in G.X],
                B=float(Bd.sum() * G.h ** 3 / (2 * np.pi ** 2)),
                F_energy=float(-2 * g * lam * np.sum(W * e2) * fac),   # -dU/dX for U = sum W e2 (first order in shift)
                F_control=float(f * lam * np.sum(dq) * fac), Q=float(np.sum(dq) * fac),
                E2=float(np.sum(e2) * fac))


def step(G, p, v, t, dt, g, f):
    def a(q, w, tt):
        W, Wd = weights(G, tt, g); return accel(G, q, w, W, Wd, f * ramp(tt), 2)
    k1p = v; k1v = a(p, v, t)
    k2p = v + dt * k1v / 2; k2v = a(G.fix(p + dt * k1p / 2), k2p, t + dt / 2)
    k3p = v + dt * k2v / 2; k3v = a(G.fix(p + dt * k2p / 2), k3p, t + dt / 2)
    k4p = v + dt * k3v; k4v = a(G.fix(p + dt * k3p), k4p, t + dt)
    q = G.fix(p + dt * (k1p + 2 * k2p + 2 * k3p + k4p) / 6)
    return q, ref.tangent(q, v + dt * (k1v + 2 * k2v + 2 * k3v + k4v) / 6) * G.inner


# ---------- initial states ----------
def rational_map(G, B, scale):
    X = np.array(G.X); r = G.r; rs = np.maximum(r, 1e-12)
    ct = X[2] / rs; st = np.sqrt(np.maximum(0, 1 - ct ** 2)); ph = np.arctan2(X[1], X[0])
    z = np.where(ct > -1 + 1e-14, st / np.maximum(1 + ct, 1e-14), 1e14) * np.exp(1j * ph)   # tan(theta/2) e^{i phi}
    R = z ** B; a2_ = np.abs(R) ** 2
    n = np.array([2 * R.real, 2 * R.imag, 1 - a2_]) / (1 + a2_)
    f = np.pi * np.exp(-r ** 2 / (2 * scale ** 2))
    p = np.empty((4,) + r.shape); p[0] = np.cos(f); p[1:] = np.sin(f) * n
    return G.fix(p)


def restpath(B):
    return RAW / f'rest_B{B}_h{H:g}_L{L:g}.npy'


def static(B, maxiter):
    t0 = time.process_time(); G = ref.Grid(H, L)
    p = G.hedgehog(ref.PROF) if B == 1 else rational_map(G, 2, 1.3)
    a2.evaluate = lambda G_, q, mode=0: evaluate(G_, q)
    p, info = a2.optimize(G, p, 0, maxiter)
    np.save(restpath(B), p); fac = G.h ** 3 / NORM
    Hd, _, Bd, e2 = fields(G, p, np.zeros_like(p), np.ones(G.r.shape))
    wall = (np.abs(np.array(G.X)).max(0) > L - 1.5)
    info.update(B=B, h=H, L=L, actual_L=float(G.ax[-1]), geometric=a2.geometric(p), B_density=float(Bd.sum() * G.h ** 3 / (2 * np.pi ** 2)),
                Q=float(np.sum((1 - p[0]) ** 2) * fac), wall_energy_fraction=float(Hd[wall].sum() / Hd.sum()),
                bcentroid=[float(np.sum(Bd * x) / Bd.sum()) for x in G.X], centroid=[float(np.sum(Hd * x) / Hd.sum()) for x in G.X],
                pred_energy_a_over_g=float(-2 * info['E2'] / info['E']), sha256=hashlib.sha256(restpath(B).read_bytes()).hexdigest(),
                cpu_total=time.process_time() - t0)
    save(f'static_B{B}', info); print('static', B, {k: info[k] for k in ('E', 'E2', 'E4', 'ratio', 'projected_gradient_max', 'iterations', 'Q', 'wall_energy_fraction', 'pred_energy_a_over_g')}, info['geometric']['degrees'], flush=True)


def run(B, mode):
    t0 = time.process_time(); G = ref.Grid(H, L)
    st = json.loads((OUT / f'static_B{B}.json').read_text()); assert hashlib.sha256(restpath(B).read_bytes()).hexdigest() == st['sha256']
    p = np.load(restpath(B)); v = np.zeros_like(p); fac = G.h ** 3 / NORM
    s1 = json.loads((OUT / 'static_B1.json').read_text())
    fctl = 2 * GRAV * s1['E2'] / s1['Q']   # one constant for both species: B = 1 control force = its energy-coupling force
    g = GRAV if mode == 'energy' else 0.; f = -fctl if mode == 'control' else 0.   # sign: control pushes toward -x like gravity
    dt = .05 * H; steps = round(T_END / dt); every = max(1, round(.1 / dt)); rows = []; d0 = a2.geometric(p)
    tag = f'B{B}_{mode}'
    for i in range(steps + 1):
        t = i * dt
        if i % every == 0 or i == steps:
            rows.append(measure(G, p, v, t, g, f))
            if (i // every) % 10 == 0 and i % every == 0: print(tag, round(t, 2), round(time.process_time() - t0), rows[-1]['bcentroid'][0], flush=True)
        if i < steps:
            p, v = step(G, p, v, t, dt, g, f)
    end = RAW / (tag + '_end.npy'); np.save(end, p)
    r = dict(B=B, mode=mode, h=H, L=L, T=steps * dt, dt=dt, g=g, f=f, f_control_constant=fctl, E=st['E'], E2=st['E2'], Q=st['Q'],
             rows=rows, degree_initial=d0, degree_final=a2.geometric(p), rest_sha256=st['sha256'],
             end_sha256=hashlib.sha256(end.read_bytes()).hexdigest(), cpu_s=time.process_time() - t0,
             platform=f'Python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__} {platform.platform()}')
    save(tag, r); print('DONE', tag, r['cpu_s'], flush=True)


def validate():
    t0 = time.process_time(); rng = np.random.default_rng(2651); G = ref.Grid(.25, 1)
    p = G.fix(G.vacuum() + .17 * rng.normal(size=(4, G.n, G.n, G.n))); v = ref.tangent(p, rng.normal(scale=.05, size=p.shape)) * G.inner
    one = np.ones(G.r.shape); zero = np.zeros_like(one)
    # 1. W = 1 equals the 26B2 kernel
    b2.POWER = 2
    eq_E = float(abs(evaluate(G, p)[0] - b2.evaluate(G, p)[0]).max()); eq_g = float(abs(evaluate(G, p)[1] - b2.evaluate(G, p)[1]).max())
    eq_a = float(abs(accel(G, p, v, one, zero, .07, 2) - b2.accel(G, p, v, .07)).max())
    # 2. independent NumPy reference with nontrivial W, Wd, f
    W = np.exp(2 * .3 * G.X[0] + .2 * G.X[1]); Wd = .4 * G.X[0] * W - .1 * G.X[2]
    er, gr, ar = numpy_ref(G, p, v, W, Wd, .07); e, gk = evaluate(G, p, W); ak = accel(G, p, v, W, Wd, .07, 2)
    # 3. directional derivative of the static energy with W
    q = ref.tangent(p, rng.normal(size=p.shape)) * G.inner; q /= np.linalg.norm(q); eps = 1e-5
    fd = (evaluate(G, G.fix(p + eps * q), W)[0].sum() - evaluate(G, G.fix(p - eps * q), W)[0].sum()) / (2 * eps)
    dg = abs(fd - np.sum(gk * q)) / abs(fd)
    # 4. energy conservation rate with static W and control force
    f = .07; fac = G.h ** 3 / NORM
    def total(pp, vv_):
        Hd = fields(G, pp, vv_, W)[0]; return Hd.sum() * fac - f * np.sum(G.X[0] * (1 - pp[0]) ** 2) * fac
    a = accel(G, p, v, W, zero, f, 2); ep = 1e-6
    rate = (total(G.fix(p + ep * v), ref.tangent(G.fix(p + ep * v), v + ep * a)) - total(G.fix(p - ep * v), ref.tangent(G.fix(p - ep * v), v - ep * a))) / (2 * ep)
    # 5. mirror y -> -y with W depending on x only
    Wx = np.exp(.5 * G.X[0]); pr = p[:, :, ::-1, :].copy(); vr = v[:, :, ::-1, :].copy()
    mir = float(abs(accel(G, pr, vr, Wx, zero) - accel(G, p, v, Wx, zero)[:, :, ::-1, :]).max())
    # 6. time-dependent W: dH/dt must equal -dL/dt|explicit = -sum Wd (vv - e2)  (checks the Wd*v term in the kernel)
    G2 = ref.Grid(.25, 1.5); p2 = G2.fix(G2.hedgehog(ref.PROF)); v2 = np.zeros_like(p2); fac2 = G2.h ** 3 / NORM
    Wt = lambda tt: (np.exp(.6 * tt * G2.X[0]), .6 * G2.X[0] * np.exp(.6 * tt * G2.X[0]))
    def ledger_rate(pp, ww, tt):
        Wn, Wdn = Wt(tt); _, _, _, e2 = fields(G2, pp, ww, Wn); return float(-np.sum(Wdn * (bm.dot(ww, ww) - e2)) * fac2)
    dt2 = .005; Hs = []; rs = []; tt = 0.
    for k in range(41):
        Wn, _ = Wt(tt); Hs.append(fields(G2, p2, v2, Wn)[0].sum() * fac2); rs.append(ledger_rate(p2, v2, tt))
        if k == 40: break
        def acc(q, w, t_):
            Wn, Wdn = Wt(t_); return accel(G2, q, w, Wn, Wdn, 0., 2)
        k1v = acc(p2, v2, tt); k2p = v2 + dt2 * k1v / 2; k2v = acc(G2.fix(p2 + dt2 * v2 / 2), k2p, tt + dt2 / 2)
        k3p = v2 + dt2 * k2v / 2; k3v = acc(G2.fix(p2 + dt2 * k2p / 2), k3p, tt + dt2 / 2); k4p = v2 + dt2 * k3v; k4v = acc(G2.fix(p2 + dt2 * k3p), k4p, tt + dt2)
        p2 = G2.fix(p2 + dt2 * (v2 + 2 * k2p + 2 * k3p + k4p) / 6); v2 = ref.tangent(p2, v2 + dt2 * (k1v + 2 * k2v + 2 * k3v + k4v) / 6) * G2.inner; tt += dt2
    rs = np.array(rs); integ = float(np.sum((rs[1:] + rs[:-1]) / 2) * dt2)
    ledger = dict(H_change=float(Hs[-1] - Hs[0]), integrated_rate=integ, rel_mismatch=float(abs(Hs[-1] - Hs[0] - integ) / abs(Hs[-1] - Hs[0])))
    vac = G.vacuum()
    out = dict(ramp_ledger=ledger, W1_vs_26B2_E=eq_E, W1_vs_26B2_grad=eq_g, W1_vs_26B2_accel=eq_a, numpy_E=float(abs(e - er).max()), numpy_grad=float(abs(gk - gr).max()),
               numpy_accel=float(abs(ak - ar).max()), directional_grad_rel=float(dg), energy_rate=float(rate), mirror_accel=mir,
               vacuum_E=float(evaluate(G, vac, W)[0].sum()), vacuum_accel=float(abs(accel(G, vac, zero * vac, W, Wd, f)).max()),
               platform=f'Python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__} {platform.platform()}', cpu_s=time.process_time() - t0)
    save('validation', out); print(out, flush=True)
    assert max(eq_E, eq_g, eq_a) < 1e-12 and max(out['numpy_E'], out['numpy_grad'], out['numpy_accel'], mir, out['vacuum_E'], out['vacuum_accel']) < 1e-10
    assert dg < 1e-5 and abs(rate) < 1e-7 and ledger['rel_mismatch'] < 1e-3


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('phase', choices=['validate', 'static', 'run'])
    ap.add_argument('--B', type=int, choices=[1, 2], default=1); ap.add_argument('--maxiter', type=int, default=300)
    ap.add_argument('--mode', choices=['energy', 'control', 'null'], default='null'); ap.add_argument('--cpu-limit', type=int, default=1500)
    ap.add_argument('--h', type=float, choices=[0.2, 0.15], default=0.2)
    a = ap.parse_args(); configure(a.h)
    if resource is not None and a.cpu_limit > 0:
        resource.setrlimit(resource.RLIMIT_CPU, (a.cpu_limit, a.cpu_limit + 1))
    load()
    if a.phase == 'validate': validate()
    elif a.phase == 'static': static(a.B, a.maxiter)
    else: run(a.B, a.mode)


if __name__ == '__main__':
    main()
