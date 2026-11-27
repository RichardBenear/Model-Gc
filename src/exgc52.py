"""EX-GC52: a moving source in the 51 model. 52A charge of a moving Dirac packet; 52B steady comoving field of a moving
source (full rule with d/dt -> -V d/dx); 52C semiclassical forces (co-moving probe R(V), probe at rest); 52D dynamic
co-moving pair. Plan: spec/EX-GC52_plan.md. EXPLORATORY. Our construction; not LaViolette's Model G.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np
from scipy.fft import dstn, idstn

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G
import exgc49 as F49

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC52"
M, SIG = 0.6, 5.0
VS = (0.0, 0.1, 0.2, 0.3)
fN, Nf = G.fN, G.Nf


def k_for_v(V, m, r, N0=1.0):
    """k0 along x with lattice group velocity V (bisection on the rising branch)."""
    if V == 0:
        return 0.0
    v = lambda k: (G.disp(np.array([k + 1e-6, 0, 0]), N0, m, r) - G.disp(np.array([k - 1e-6, 0, 0]), N0, m, r)) / 2e-6
    lo, hi = 0.0, 0.7
    for _ in range(80):
        mid = (lo + hi) / 2; lo, hi = (mid, hi) if v(mid) < V else (lo, mid)
    return (lo + hi) / 2


def packet_k(n, m, r, x0, sigma, k0, N0=1.0, chi=(1, 0, 0, 0)):
    """positive-energy Gaussian packet centred at x0 with mean momentum k0 (vector), zero on the boundary."""
    k = 2 * np.pi * np.fft.fftfreq(n); K = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
    S = np.sin(K); W = np.sum(1 - np.cos(K), -1)
    h = N0 ** 2 * np.einsum("...j,jab->...ab", S, np.array(G.AL)) + (N0 * m + N0 ** 2 * r * W)[..., None, None] * G.BE
    E = G.disp(K, N0, m, r); P = 0.5 * (np.eye(4) + h / E[..., None, None])
    d = K - np.asarray(k0, float); d = (d + np.pi) % (2 * np.pi) - np.pi
    g = np.exp(-np.sum(d ** 2, -1) * sigma ** 2 / 2) * np.exp(-1j * (K @ np.asarray(x0, float)))
    a = P @ np.asarray(chi, complex); a = a / np.linalg.norm(a, axis=-1, keepdims=True) * g[..., None]
    psi = np.fft.ifftn(a, axes=(0, 1, 2))
    for ax in range(3):
        sl = [slice(None)] * 4; sl[ax] = 0; psi[tuple(sl)] = 0; sl[ax] = -1; psi[tuple(sl)] = 0
    w = np.abs(g) ** 2; w = w / w.sum()
    return psi / np.linalg.norm(psi), K, w


def EN_avg(K, w, N0, m, r):
    h = 1e-5; return float(np.sum(w * (G.disp(K, N0 + h, m, r) - G.disp(K, N0 - h, m, r)) / (2 * h)))


# ------------------------------------------------------------------ 52A
def part_A(n=96):
    rows = {}
    for r in (0.0, 1.0):
        for V in VS:
            k0 = k_for_v(V, M, r); c = n // 2; psi, K, w = packet_k(n, M, r, (c, c, c), SIG, (k0, 0, 0))
            one = np.ones((n, n, n)); Q = float(G.charge(psi, one, M, r).sum()); E = G.energy(psi, one, M, r)
            pred = EN_avg(K, w, 1.0, M, r); gam = 1 / np.sqrt(1 - V * V)
            rows[f"{r}_{V}"] = dict(r=r, V=V, k0=k0, Q=Q, E=E, EN_pred=pred, rel=Q / pred - 1, Q_over_cont=Q / (gam * M * (1 + V * V)), Q_over_E=Q / E)
            print("52A r %.1f V %.1f k0 %.4f Q %.6f E %.6f pred %.6f rel %.1e Q/(gam m(1+V2)) %.4f" % (r, V, k0, Q, E, pred, Q / pred - 1, Q / (gam * M * (1 + V * V))), flush=True)
    return rows


# ------------------------------------------------------------------ steady comoving solver
def Tx(N, f):
    """(1/N) d_x[(1/N) d_x f] on the interior, compact stencil with link-mean N."""
    I = (slice(1, -1),) * 3; Np = (N[2:, 1:-1, 1:-1] + N[1:-1, 1:-1, 1:-1]) / 2; Nm = (N[:-2, 1:-1, 1:-1] + N[1:-1, 1:-1, 1:-1]) / 2
    fx = f[I]; return ((f[2:, 1:-1, 1:-1] - fx) / Np - (fx - f[:-2, 1:-1, 1:-1]) / Nm) / N[I]


def lap_aniso(f, V):
    I = (slice(1, -1),) * 3; fx = f[I]
    lx = f[2:, 1:-1, 1:-1] + f[:-2, 1:-1, 1:-1] - 2 * fx; ly = f[1:-1, 2:, 1:-1] + f[1:-1, :-2, 1:-1] - 2 * fx; lz = f[1:-1, 1:-1, 2:] + f[1:-1, 1:-1, :-2] - 2 * fx
    return (1 - V * V) * lx + ly + lz


def moving(psis, ms, r, kappa, n, centre, V, tol=1e-11, maxit=300, qfix=None, gauss=True):
    """steady comoving solution of op(N, f) - V^2 Tx(N, f) = kappa q; Dirichlet faces from the linear moving point solution."""
    c = np.array(centre, float); x = np.arange(n); X = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1) - c
    gam = 1 / np.sqrt(1 - V * V); Rg = np.sqrt((gam * X[..., 0]) ** 2 + X[..., 1] ** 2 + X[..., 2] ** 2); Rg = np.where(Rg == 0, 1e-9, Rg)
    N = np.ones((n, n, n)); mI = n - 2; k = np.arange(1, mI + 1); l1 = 2 * np.cos(np.pi * k / (mI + 1)) - 2
    lam = (1 - V * V) * l1[:, None, None] + l1[None, :, None] + l1[None, None, :]; I = (slice(1, -1),) * 3; Fh, Gh = [], []
    for it in range(maxit):
        q = qfix if qfix is not None else sum(G.charge(p, N, mm, r) for p, mm in zip(psis, ms)); Q = q.sum()
        Nb = np.exp(-kappa * Q * gam / (4 * np.pi * Rg)) if gauss else np.ones_like(Rg)
        f = fN(N); fb = fN(Nb); f[0], f[-1], f[:, 0], f[:, -1], f[:, :, 0], f[:, :, -1] = fb[0], fb[-1], fb[:, 0], fb[:, -1], fb[:, :, 0], fb[:, :, -1]
        N = Nf(f); opf = F49.op(N, f) - V * V * Tx(N, f)
        res = opf - kappa * q[I]; rn = float(np.abs(res).max())
        if rn < tol:
            break
        fbz = f.copy(); fbz[I] = 0; bc = lap_aniso(fbz, V)
        R = opf - lap_aniso(f, V); g = idstn(dstn(kappa * q[I] - R - bc, type=1) / lam, type=1)
        x0 = f[I].ravel(); gx = g.ravel(); Fk = gx - x0; Fh.append(Fk); Gh.append(gx)
        if len(Fh) > 6:
            Fh.pop(0); Gh.pop(0)
        if len(Fh) > 1:
            dF = np.array([Fh[i + 1] - Fh[i] for i in range(len(Fh) - 1)]).T; dG = np.array([Gh[i + 1] - Gh[i] for i in range(len(Gh) - 1)]).T
            xn = gx - dG @ np.linalg.lstsq(dF, Fk, rcond=None)[0]
        else:
            xn = gx
        f[I] = xn.reshape(mI, mI, mI); N = Nf(f)
    return N, dict(iters=it + 1, resid=rn, Q=float(Q))


# ------------------------------------------------------------------ 52B + 52C
def part_BC(n=129, u=0.01, rho_p=16, Xpass=32):
    c = n // 2; res = {}; Q0 = None; fields = {}
    for V in VS:
        k0 = k_for_v(V, M, 0.0); psi, K, w = packet_k(n, M, 0.0, (c, c, c), SIG, (k0, 0, 0))
        if Q0 is None:
            Q0 = float(G.charge(psi, np.ones((n, n, n)), M, 0.0).sum()); kap = 4 * np.pi * 12 * u / Q0
        t0 = time.time(); N, info = moving([psi], [M], 0.0, kap, n, (c, c, c), V); fields[V] = N
        gam = 1 / np.sqrt(1 - V * V); Q = info["Q"]; psi_ = np.log(N); rr = np.arange(12, 37)
        amp_rho = psi_[c, c + rr, c] * rr / (-kap * Q * gam / (4 * np.pi)); amp_rho2 = psi_[c, c, c + rr] * rr / (-kap * Q * gam / (4 * np.pi))
        amp_x = psi_[c + rr, c, c] * rr / (-kap * Q / (4 * np.pi)); amp_xm = psi_[c - rr, c, c] * rr / (-kap * Q / (4 * np.pi))
        contraction = float(np.mean(psi_[c, c + rr, c] * rr) / np.mean(0.5 * (psi_[c + rr, c, c] + psi_[c - rr, c, c]) * rr))
        # 52C co-moving probe: semiclassical transverse force at (c, c + rho_p, c), probe = same packet shape, velocity V
        pp, Kp, wp = packet_k(n, M, 0.0, (c, c + rho_p, c), SIG, (k0, 0, 0)); rho_ = np.sum(np.abs(pp) ** 2, -1)
        dyN = float(np.sum(rho_ * (np.roll(N, -1, 1) - np.roll(N, 1, 1)) / 2)); Nloc = float(np.sum(rho_ * N))
        EN = EN_avg(Kp, wp, Nloc, M, 0.0); pdot = -EN * dyN
        # probe at rest beside the passing source: forces along the pass at transverse distance rho_p (source frame X = probe - source)
        Xs = np.arange(-Xpass, Xpass + 1); EN0 = EN_avg(*packet_k(n, M, 0.0, (c, c, c), SIG, (0, 0, 0))[1:], 1.0, M, 0.0)
        Fy = np.array([-EN0 * (N[c + X, c + rho_p + 1, c] - N[c + X, c + rho_p - 1, c]) / 2 for X in Xs])
        Fx = np.array([-EN0 * (N[c + X + 1, c + rho_p, c] - N[c + X - 1, c + rho_p, c]) / 2 for X in Xs])
        res[str(V)] = dict(V=V, k0=k0, gamma=gam, kappa=kap, info=info, Q=Q, N_centre=float(N[c, c, c]),
                           amp_rho_mean=float(np.mean(amp_rho)), amp_rho_z_mean=float(np.mean(amp_rho2)), amp_x_mean=float(np.mean(amp_x)), amp_xminus_mean=float(np.mean(amp_xm)),
                           amp_rho_spread=float(np.ptp(amp_rho)), contraction=contraction, contraction_over_gamma=contraction / gam,
                           probe=dict(EN=EN, dyN=dyN, Nloc=Nloc, pdot=pdot), rest_probe=dict(X=Xs.tolist(), Fy=Fy.tolist(), Fx=Fx.tolist(), EN0=EN0), wall=time.time() - t0)
        print("52B V %.1f it %d res %.1e Q %.6f | amp rho %.4f (z %.4f) amp x %.4f/%.4f | contraction/gam %.4f | pdot %.4e  %.0fs" % (
            V, info["iters"], info["resid"], Q, np.mean(amp_rho), np.mean(amp_rho2), np.mean(amp_x), np.mean(amp_xm), contraction / gam, pdot, time.time() - t0), flush=True)
    p0 = res["0.0"]["probe"]["pdot"]; Q0s = res["0.0"]["Q"]
    for V in VS[1:]:
        R = res[str(V)]; gam = R["gamma"]
        R["R"] = gam * R["probe"]["pdot"] / p0
        R["R_pred_lattice"] = gam ** 2 * (R["Q"] / Q0s) * (R["probe"]["EN"] / res["0.0"]["probe"]["EN"])
        R["R_cont"] = gam ** 4 * (1 + V * V) ** 2
        # probe at rest: transverse impulse over the pass vs GR linearised point formula with m = Q(0) (h00 = 2 G m gam^2 (1+V^2)/sqrt(gam^2 X^2 + rho^2))
        X = np.array(R["rest_probe"]["X"], float); kap = R["kappa"]; Gm = kap * Q0s / (4 * np.pi)
        dPhi_dy = -Gm * gam ** 2 * (1 + V * V) * rho_p / (gam ** 2 * X ** 2 + rho_p ** 2) ** 1.5
        Fy_GR = R["rest_probe"]["EN0"] * dPhi_dy                       # force = -E_N dN/dy, N = 1 - Phi
        R["rest_impulse_ratio"] = float(np.trapezoid(R["rest_probe"]["Fy"], X) / np.trapezoid(Fy_GR, X))
        dPhi_dX = -Gm * gam ** 2 * X * gam ** 2 / (gam ** 2 * X ** 2 + rho_p ** 2) ** 1.5
        R["rest_Fx_GR_harmonic"] = (R["rest_probe"]["EN0"] * dPhi_dX * (1 - 3 * V * V) / (1 + V * V) * (1 + V * V)).tolist()   # GR: gam^2 (1 - 3V^2) dPhi0/dX (unscored)
        print("52C V %.1f R %.5f pred(lattice) %.5f cont %.5f | rest transverse impulse / GR %.4f" % (V, R["R"], R["R_pred_lattice"], R["R_cont"], R["rest_impulse_ratio"]), flush=True)
    np.save(OUT / "partial_fields_B.npy", np.stack([fields[V][:, :, n // 2] for V in VS]))
    return res


# ------------------------------------------------------------------ 52D
def evolve_g(psis, ms, r, f, g, kappa, T, dt, record):
    I = (slice(1, -1),) * 3
    def rhs(ps, f, g):
        N = Nf(f); dps = [-1j * G.Hm(p, N, mm, r) for p, mm in zip(ps, ms)]
        q = sum(G.charge(p, N, mm, r) for p, mm in zip(ps, ms)); df = N * g; dg = np.zeros_like(g); dg[I] = N[I] * (F49.op(N, f) - kappa * q[I])
        for a in range(3):
            sl = [slice(None)] * 3; sl[a] = 0; df[tuple(sl)] = 0; sl[a] = -1; df[tuple(sl)] = 0
        return dps, df, dg
    t = 0.0; hist = [record(psis, f, 0.0)]
    while t < T - 1e-9:
        def add(ps, f, g, k, h):
            return [p + h * kp for p, kp in zip(ps, k[0])], f + h * k[1], g + h * k[2]
        k1 = rhs(psis, f, g); k2 = rhs(*add(psis, f, g, k1, dt / 2)); k3 = rhs(*add(psis, f, g, k2, dt / 2)); k4 = rhs(*add(psis, f, g, k3, dt))
        psis = [p + dt / 6 * (a + 2 * b + 2 * c + d) for p, a, b, c, d in zip(psis, k1[0], k2[0], k3[0], k4[0])]
        f = f + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]); g = g + dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]); t += dt; hist.append(record(psis, f, t))
    return psis, f, hist


def momentum(p, ax):
    n = p.shape[0]; pk = np.fft.fftn(p, axes=(0, 1, 2)); w = np.sum(np.abs(pk) ** 2, -1); k = 2 * np.pi * np.fft.fftfreq(n)
    sh = [1, 1, 1]; sh[ax] = -1; return float(np.sum(w * k.reshape(sh)) / w.sum())


def initial_field(ps, ms, r, kap, n, c, V):
    """steady comoving field of the rigid frozen packets: q from N = 1, then two refreshes of q in the solved N."""
    q = sum(G.charge(p, np.ones((n, n, n)), mm, r) for p, mm in zip(ps, ms))
    for _ in range(3):
        N, info = moving(ps, ms, r, kap, n, (c, c, c), V, qfix=q, tol=1e-10)
        q = sum(G.charge(p, N, mm, r) for p, mm in zip(ps, ms))
    f = fN(N); g = np.zeros_like(f)
    g[1:-1] = -V * (f[2:] - f[:-2]) / 2 / N[1:-1]          # g = (1/N) d_t f with d_t = -V d_x
    return f, g, info


def part_D(n=65, sep=14, u=0.01, T=16.0, dt=0.1, Vs=(0.0, 0.2)):
    c = n // 2; res = {}
    for V in Vs:
        k0 = k_for_v(V, M, 0.0); p1, _, _ = packet_k(n, M, 0.0, (c, c - sep // 2, c), SIG, (k0, 0, 0)); p2, _, _ = packet_k(n, M, 0.0, (c, c + sep // 2, c), SIG, (k0, 0, 0))
        if V == Vs[0]:
            Q1 = float(G.charge(p1, np.ones((n, n, n)), M, 0.0).sum()); kap = 4 * np.pi * sep * u / Q1
        runs = {}
        for lab, ps in (("pair", [p1, p2]), ("one", [p1])):
            t0 = time.time(); f, g, info = initial_field(ps, [M] * len(ps), 0.0, kap, n, c, V)
            _, _, hist = evolve_g(ps, [M] * len(ps), 0.0, f, g, kap, T, dt, lambda ps, f, t: (t, [momentum(p, 1) for p in ps], [G.energy(p, Nf(f), M, 0.0) for p in ps]))
            runs[lab] = hist; print("52D V %.1f %s %.0fs (init it %d res %.1e)" % (V, lab, time.time() - t0, info["iters"], info["resid"]), flush=True)
        t = np.array([h[0] for h in runs["pair"]]); P1 = np.array([h[1][0] for h in runs["pair"]]) - np.array([h[1][0] for h in runs["one"]])
        F1 = float(np.polyfit(t, P1, 2)[1]); res[str(V)] = dict(V=V, k0=k0, kappa=kap, F1=F1, t=t.tolist(), P1=P1.tolist())
        print("52D V %.1f F1 %.5e" % (V, F1), flush=True)
    if "0.0" in res:
        for V in Vs[1:]:
            gam = 1 / np.sqrt(1 - V * V); res[str(V)]["R_dyn"] = gam * res[str(V)]["F1"] / res["0.0"]["F1"]; print("52D R_dyn(%.1f) %.5f" % (V, res[str(V)]["R_dyn"]), flush=True)
    return res


def save(name, obj, t0):
    OUT.mkdir(parents=True, exist_ok=True)
    rec = dict(study="EX-GC52", part=name, label="EXPLORATORY", platform=platform.platform(), python=platform.python_version(), numpy=np.__version__,
               wall_s=time.time() - t0, result=obj)
    (OUT / f"{name}.json").write_text(json.dumps(rec, indent=1, default=float))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for p in sys.argv[1:]:
        t0 = time.time(); save(p, {"A": part_A, "BC": part_BC, "D": part_D}[p](), t0); print(p, "done %.0fs" % (time.time() - t0), flush=True)
