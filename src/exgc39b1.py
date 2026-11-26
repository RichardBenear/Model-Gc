"""EX-GC39B B1: do emergent particles keep their proper size, clock and fall under 37B's event weights?
Plan: spec/EX-GC39B_plan.md. Our construction, not LaViolette's Model G.

37B rule: one-site events at rate N, two-site coincidences at rate N^2. Weak-field metric -N^2 dt^2 + N^-2 dx^2,
proper length dl = dx / N. Equivalence principle: coordinate size ~ N, proper size and proper rates fixed.
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import ex008_longrange as E8

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC39B"


# ---------------------------------------------------------------- B1a: Gc1 reaction-diffusion particle
class Scaled(E8.Model):
    """Steady Gc1 with reactions at rate N and hops at rate N^2: divide by N, so every diffusion term carries N."""
    def __init__(self, R, dr, k, cx, cy, N):
        super().__init__(R, dr, k, cx, cy)
        self.LapN = (N * self.LapN).tocsr()
        self.LapD = (N * self.LapD).tocsr()
        self.bD = N * self.bD


def b1a(R=60.0, dr=0.05, k=0.1):
    M1, u1, res, ok = E8.particle(R, dr, k, 1, 0)
    base = E8.analyse(M1, u1)
    rows = []
    for N in (0.8, 0.9, 1.0, 1.1):
        M = Scaled(R, dr, k, 1, 0, N)
        s = np.sqrt(N)
        u0 = np.array([np.interp(M.r / s, M1.r, u1[i]) for i in range(4)])
        u0[3] = np.interp(M.r / s, M1.r, u1[3] - E8.GINF, right=0.0) + E8.GINF
        u, res, ok = M.newton(u0)
        a = E8.analyse(M, u)
        lam = E8.leading_eigs(M, u, 0)[0][0]           # rate in units of the reaction rate (= proper rate)
        rows.append(dict(N=N, ok=bool(ok), residual=res, W_rms_coord=a["W_radius_rms"],
                         W_rms_proper=a["W_radius_rms"] / N, core_W=a["core_W"], lam_proper=lam, lam_coord=N * lam))
        print(rows[-1], flush=True)
    Ns = np.array([r["N"] for r in rows]); Rp = np.array([r["W_rms_proper"] for r in rows])
    expo = float(np.polyfit(np.log(Ns), np.log(Rp), 1)[0])
    lam = np.array([r["lam_proper"] for r in rows])
    return dict(rows=rows, proper_radius_exponent=expo, lam_proper_max_rel_spread=float(np.ptp(lam) / abs(lam.mean())),
                check_N1=dict(core_W_ex008=base["core_W"], core_W_scaled=rows[2]["core_W"],
                              cons_ex008=base["consumption_M"]))


# ---------------------------------------------------------------- B1b: phi^4 kink on a 1+1D lattice
def kink_run(lam, v, Nfun, rule, h=0.1, dt=0.02, L=30.0, T=40.0, perturb=0.0):
    x = np.arange(-L, L + h / 2, h)
    N = Nfun(x)
    w = np.sqrt(2) / (v * np.sqrt(lam))
    phi = v * np.tanh(x / (w * N[np.argmin(np.abs(x))]) if rule == "handshake" else x / w)
    phi = phi + perturb * np.exp(-(x / w) ** 2) / np.cosh(x / w) ** 2    # shape-mode kick (odd shape mode ~ sech^2 tanh;
    old = phi.copy()                                                      # even bump excites the breathing/shape sector)
    face = N[1:] * N[:-1]
    pot_w = N ** 2 if rule == "handshake" else N

    def acc(p):
        a = np.zeros_like(p)
        flux = face * (p[1:] - p[:-1]) / h ** 2
        a[1:-1] = N[1:-1] ** 2 * (flux[1:] - flux[:-1]) - pot_w[1:-1] * lam * p[1:-1] * (p[1:-1] ** 2 - v * v)
        return a
    def energy(p, pold):
        pt = (p - pold) / dt
        return float(np.sum(0.5 * pt[1:-1] ** 2 / N[1:-1] ** 2 * h) + np.sum(0.5 * face * ((p[1:] - p[:-1]) / h) ** 2 * h)
                     + np.sum((pot_w / N ** 2)[1:-1] * lam / 4 * (p[1:-1] ** 2 - v * v) ** 2 * h))
    phi = old + 0.5 * dt * dt * acc(old)
    ts, xc, mid = [], [], []
    E0 = energy(phi, old)
    i0 = np.argmin(np.abs(x))
    for n in range(int(T / dt)):
        new = 2 * phi - old + dt * dt * acc(phi)
        old, phi = phi, new
        j = np.where(np.sign(phi[:-1]) != np.sign(phi[1:]))[0]
        if len(j):
            j = j[np.argmin(np.abs(x[j]))]
            xc.append(x[j] - phi[j] * h / (phi[j + 1] - phi[j]))
        else:
            xc.append(np.nan)
        ts.append((n + 1) * dt)
        mid.append(phi[i0 + int(round(w / h))])
    return np.array(ts), np.array(xc), np.array(mid), energy(phi, old), E0, x, phi


def width_of(x, phi, v):
    """Coordinate width: distance between phi = +-v tanh(1/2) crossings, divided by 2 artanh(tanh(.5))=1 -> w."""
    t = v * np.tanh(0.5)
    def cross(level):
        j = np.where((phi[:-1] - level) * (phi[1:] - level) <= 0)[0][0]
        return x[j] + (level - phi[j]) * (x[j + 1] - x[j]) / (phi[j + 1] - phi[j])
    return cross(t) - cross(-t)


def relax_static(lam, v, N, rule, h=0.1, L=30.0):
    """Static kink in uniform N by gradient flow of the static equation; returns coordinate width."""
    x = np.arange(-L, L + h / 2, h)
    w = np.sqrt(2) / (v * np.sqrt(lam))
    p = v * np.tanh(x / w)
    pot_w = N ** 2 if rule == "handshake" else N
    for _ in range(200000):
        r = np.zeros_like(p)
        r[1:-1] = N ** 2 * N * N * (p[2:] - 2 * p[1:-1] + p[:-2]) / h ** 2 - pot_w * lam * p[1:-1] * (p[1:-1] ** 2 - v * v)
        p += 0.2 * h * h / (N ** 4) * r
        if np.max(np.abs(r)) < 1e-11:
            break
    return width_of(x, p, v), x, p


def shape_freq(lam, v, N, rule, h=0.1):
    """Small-oscillation frequencies about the static kink in uniform N (coordinate), from the linearised operator."""
    wd, x, p = relax_static(lam, v, N, rule, h)
    pot_w = N ** 2 if rule == "handshake" else N
    n = len(x) - 2
    main = 2 * N ** 4 / h ** 2 + pot_w * lam * (3 * p[1:-1] ** 2 - v * v)
    off = -N ** 4 / h ** 2 * np.ones(n - 1)
    from scipy.linalg import eigh_tridiagonal
    ev = eigh_tridiagonal(main, off, select="i", select_range=(0, 1))[0]
    return wd, np.sqrt(np.clip(ev, 0, None))


def b1b():
    out = dict(uniform={}, fall={})
    for rule in ("handshake", "lapse_only"):
        for lam, v in ((1.0, 1.0), (8.0, 0.5)):
            rows = []
            for N in (0.8, 0.9, 1.0, 1.1):
                wd, om = shape_freq(lam, v, N, rule)
                rows.append(dict(N=N, width_coord=wd, width_proper=wd / N, shape_omega_coord=float(om[1]),
                                 shape_omega_proper=float(om[1] / N), zero_mode=float(om[0])))
            Ns = np.array([r["N"] for r in rows])
            ew = float(np.polyfit(np.log(Ns), np.log([r["width_proper"] for r in rows]), 1)[0])
            eo = float(np.polyfit(np.log(Ns), np.log([r["shape_omega_proper"] for r in rows]), 1)[0])
            out["uniform"][f"{rule}_lam{lam:g}_v{v:g}"] = dict(rows=rows, proper_width_exponent=ew, proper_omega_exponent=eo)
            print(rule, lam, v, ew, eo, flush=True)
    g = 1e-3
    for rule in ("handshake", "lapse_only"):
        for lam, v in ((1.0, 1.0), (8.0, 0.5)):
            ts, xc, _, E1, E0, _, _ = kink_run(lam, v, lambda x: 1 + g * x, rule)
            tsf, xcf, _, E1f, E0f, _, _ = kink_run(lam, v, lambda x: np.ones_like(x), rule)
            d = xc - xcf          # subtract any flat-lattice drift
            c = np.polyfit(ts, d, 2)
            out["fall"][f"{rule}_lam{lam:g}_v{v:g}"] = dict(accel=float(2 * c[0]), accel_over_minus_g=float(-2 * c[0] / g),
                                                            flat_max_drift=float(np.nanmax(np.abs(xcf))),
                                                            flat_energy_rel=float(abs(E1f - E0f) / E0f))
            print(rule, lam, v, out["fall"][f"{rule}_lam{lam:g}_v{v:g}"], flush=True)
    return out


# ---------------------------------------------------------------- B1c: 3D Skyrmion hedgehog, radial
def skyrme(N, p, mpi=0.5, R=12.0, n=600):
    r = (np.arange(n) + 0.5) * R / n
    dr = R / n
    w2, w4, w0 = 1.0, N ** (p - 4), N ** -2.0

    def full(fin):
        return np.concatenate([[np.pi], fin, [0.0]])
    rf = np.concatenate([[0.0], r, [R]])

    def energy(fin):
        f = full(fin)
        d = np.diff(rf); fm = 0.5 * (f[1:] + f[:-1]); fp = np.diff(f) / d; rm = 0.5 * (rf[1:] + rf[:-1])
        s2 = np.sin(fm) ** 2; s2p = np.sin(2 * fm)
        e = (w2 * (fp ** 2 * rm ** 2 + 2 * s2) + w4 * s2 * (2 * fp ** 2 + s2 / rm ** 2)
             + w0 * 2 * mpi ** 2 * rm ** 2 * (1 - np.cos(fm)))
        dfm = w2 * 2 * s2p + w4 * (s2p * (2 * fp ** 2 + s2 / rm ** 2) + s2 * s2p / rm ** 2) + w0 * 2 * mpi ** 2 * rm ** 2 * np.sin(fm)
        dfp = w2 * 2 * fp * rm ** 2 + w4 * 4 * s2 * fp
        g = np.zeros_like(f)
        g[1:] += d * (0.5 * dfm) + dfp
        g[:-1] += d * (0.5 * dfm) - dfp
        return 4 * np.pi * float(np.sum(e * d)), 4 * np.pi * g[1:-1]
    f0 = np.pi * np.exp(-r / (1.0 * N))
    res = minimize(energy, f0, jac=True, method="L-BFGS-B", options=dict(maxiter=50000, maxfun=10 ** 6, ftol=1e-15, gtol=1e-10))
    f = full(res.x)
    fm = 0.5 * (f[1:] + f[:-1]); fp = np.diff(f) / np.diff(rf); rm = 0.5 * (rf[1:] + rf[:-1]); d = np.diff(rf)
    B = -(2 / np.pi) * np.sin(fm) ** 2 * fp * d               # baryon number per shell
    r2 = float(np.sum(B * rm ** 2) / np.sum(B))
    return dict(N=N, p=p, E=res.fun, B=float(np.sum(B)), r_rms_coord=np.sqrt(r2), r_rms_proper=np.sqrt(r2) / N,
                E_proper=res.fun / N, ok=bool(res.success))


def b1c():
    out = {}
    for p in (6, 4):
        rows = [skyrme(N, p) for N in (0.9, 1.0, 1.1)]
        Ns = np.array([r["N"] for r in rows])
        out[f"p{p}"] = dict(rows=rows,
                            proper_size_exponent=float(np.polyfit(np.log(Ns), np.log([r["r_rms_proper"] for r in rows]), 1)[0]),
                            proper_energy_exponent=float(np.polyfit(np.log(Ns), np.log([r["E_proper"] for r in rows]), 1)[0]))
        print(p, out[f"p{p}"]["proper_size_exponent"], out[f"p{p}"]["proper_energy_exponent"], flush=True)
    out["grid_check_N1_p6"] = dict(n600=skyrme(1.0, 6)["E"], n1200=skyrme(1.0, 6, n=1200)["E"])
    return out


if __name__ == "__main__":
    part = sys.argv[1]
    t0 = time.process_time()
    res = dict(b1a=b1a, b1b=b1b, b1c=b1c)[part]()
    res["cpu_seconds"] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
