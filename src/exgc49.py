"""EX-GC49: firing-rate changes carried by the same two-site coincidences as light.
49A: static mass, sum_y N_x N_y [f(N_y) - f(N_x)] = s_x on a 129^3 Dirichlet cube, four transported variables f; PPN beta
from N^2 = 1 - 2U + 2 beta U^2. 49B: light vs firing-rate perturbation packets along a chain through the f = N well.
Plan: spec/EX-GC49_plan.md (amendment A1). EXPLORATORY. Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np
from scipy.fft import dstn, idstn
from scipy.special import erf
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC49"
SIG = 2.0; Q10S = (0.015, 0.03, 0.045, 0.06)

# transported variable f, its inverse, the harmonic F(N), its inverse, predicted beta
VARS = {
    "N":      dict(f=lambda N: N,               fi=lambda f: f,                      Fi=lambda F: np.cbrt(1 + 3 * F), beta=-0.5),
    "lnN":    dict(f=lambda N: np.log(N),       fi=lambda f: np.exp(f),              Fi=lambda F: np.sqrt(1 + 2 * F),  beta=0.0),
    "invN":   dict(f=lambda N: -1 / N,          fi=lambda f: -1 / f,                 Fi=lambda F: 1 + F,               beta=0.5),
    "slow":   dict(f=lambda N: -0.5 / N ** 2,   fi=lambda f: np.sqrt(-0.5 / f),      Fi=lambda F: np.exp(F),           beta=1.0),
}


def geometry(n=129):
    c = n // 2; x = np.arange(n) - c; X, Y, Zc = np.meshgrid(x, x, x, indexing="ij"); r = np.sqrt(X ** 2 + Y ** 2 + Zc ** 2)
    return c, r


def F_cont(r, Q):
    r = np.where(r == 0, 1e-12, r); return -Q * erf(r / (np.sqrt(2) * SIG)) / r


def nb_sum(A, B):
    """sum over the 6 neighbours y of interior x: A_y * B_x-ish helpers done by caller; returns list of neighbour slices."""
    pass


def op(N, f):
    """sum_y N_x N_y (f_y - f_x) on the interior (shape n-2 cubed)."""
    I = (slice(1, -1),) * 3; Nx = N[I]; fx = f[I]; out = np.zeros_like(Nx)
    for ax in range(3):
        for d in (-1, 1):
            sl = [slice(1, -1)] * 3; sl[ax] = slice(1 + d, N.shape[ax] - 1 + d); sl = tuple(sl)
            out += Nx * N[sl] * (f[sl] - fx)
    return out


def lap(f):
    I = (slice(1, -1),) * 3; fx = f[I]; out = np.zeros_like(fx)
    for ax in range(3):
        for d in (-1, 1):
            sl = [slice(1, -1)] * 3; sl[ax] = slice(1 + d, f.shape[ax] - 1 + d); out += f[tuple(sl)] - fx
    return out


def solve(var, Q, n=129, tol=1e-11, maxit=400, verbose=False):
    V = VARS[var]; c, r = geometry(n); m = n - 2
    w = np.exp(-r ** 2 / (2 * SIG ** 2)); s = 4 * np.pi * Q * w / w.sum(); sI = s[1:-1, 1:-1, 1:-1]
    N = V["Fi"](F_cont(r, Q))                   # initial guess = continuum solution; boundary fixed to it
    f = V["f"](N)
    k = np.arange(1, m + 1); lam1 = 2 * np.cos(np.pi * k / (m + 1)) - 2
    lam = lam1[:, None, None] + lam1[None, :, None] + lam1[None, None, :]
    fb = f.copy(); fb[1:-1, 1:-1, 1:-1] = 0; bcontrib = lap(fb)    # boundary contribution to the discrete Laplacian
    hist = []; Fh = []; Gh = []
    for it in range(maxit):
        res = op(N, f) - sI; rn = float(np.abs(res).max()); hist.append(rn)
        if verbose:
            print(var, Q, it, rn, flush=True)
        if rn < tol:
            break
        R = op(N, f) - lap(f)                                     # nonlinear part
        rhs = sI - R - bcontrib
        g = idstn(dstn(rhs, type=1) / lam, type=1)
        # Anderson mixing (depth 5) on the interior f
        x0 = f[1:-1, 1:-1, 1:-1].ravel(); gx = g.ravel(); Fk = gx - x0
        Fh.append(Fk); Gh.append(gx)
        if len(Fh) > 6:
            Fh.pop(0); Gh.pop(0)
        if len(Fh) > 1:
            dF = np.array([Fh[i + 1] - Fh[i] for i in range(len(Fh) - 1)]).T; dG = np.array([Gh[i + 1] - Gh[i] for i in range(len(Gh) - 1)]).T
            gam = np.linalg.lstsq(dF, Fk, rcond=None)[0]; xn = gx - dG @ gam
        else:
            xn = gx
        f[1:-1, 1:-1, 1:-1] = xn.reshape(m, m, m); N = V["fi"](f)
    return N, r, c, hist


def samples(N, r, c):
    xs = np.arange(10, 41); ax = [N[c + i, c, c] for i in xs] + [N[c, c - i, c] for i in xs] + [N[c, c, c + i] for i in xs]
    nd = np.arange(6, 24); dg = [N[c + i, c + i, c + i] for i in nd]
    return dict(axis_r=np.concatenate([xs] * 3).astype(float), axis_N=np.array(ax), diag_r=np.sqrt(3) * nd, diag_N=np.array(dg))


def fit_one(rr, NN):
    x = 1 / rr; A = np.stack([x, x ** 2, x ** 3, x ** 4], 1); a = np.linalg.lstsq(A, NN ** 2 - 1, rcond=None)[0]
    return dict(A=float(-a[0] / 2), beta=float(2 * a[1] / a[0] ** 2), coef=a.tolist())


def fit_joint(sets):
    nS = len(sets)
    def resid(p):
        beta = p[0]; out = []
        for i, (rr, NN) in enumerate(sets):
            A, c3, c4 = p[1 + 3 * i: 4 + 3 * i]; x = 1 / rr
            out.append(1 - 2 * A * x + 2 * beta * A ** 2 * x ** 2 + c3 * x ** 3 + c4 * x ** 4 - NN ** 2)
        return np.concatenate(out)
    p0 = [0.0] + sum([[fit_one(*s)["A"], 0.0, 0.0] for s in sets], [])
    sol = least_squares(resid, p0, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    return dict(beta=float(sol.x[0]), A=[float(sol.x[1 + 3 * i]) for i in range(nS)], rms=float(np.sqrt(np.mean(sol.fun ** 2))))


def part_A():
    res = {}
    for var in VARS:
        rows = {}; sets_l = []; sets_c = []
        for q10 in Q10S:
            Q = 10 * q10; t0 = time.time(); N, r, c, hist = solve(var, Q)
            S = samples(N, r, c); rr = np.concatenate([S["axis_r"], S["diag_r"]]); NN = np.concatenate([S["axis_N"], S["diag_N"]])
            Nc = VARS[var]["Fi"](F_cont(rr, Q))
            fl = fit_one(rr, NN); fc = fit_one(rr, Nc); fax = fit_one(S["axis_r"], S["axis_N"]); fdg = fit_one(S["diag_r"], S["diag_N"])
            rows[str(q10)] = dict(Q=Q, iters=len(hist), resid=hist[-1], wall=time.time() - t0, fit_lattice=fl, fit_axis=fax, fit_diag=fdg, fit_cont=fc,
                                  A_over_Q=fl["A"] / Q, maxdiff_lat_cont=float(np.abs(NN - Nc).max()), N_center=float(N[c, c, c]))
            sets_l.append((rr, NN)); sets_c.append((rr, Nc))
            print("49A", var, q10, "it %d res %.1e | beta lat %.4f axis %.4f diag %.4f cont %.4f | A/Q %.5f | |dN| %.2e | %.0fs" % (
                len(hist), hist[-1], fl["beta"], fax["beta"], fdg["beta"], fc["beta"], fl["A"] / Q, rows[str(q10)]["maxdiff_lat_cont"], time.time() - t0), flush=True)
        jl = fit_joint(sets_l); jc = fit_joint(sets_c)
        res[var] = dict(beta_pred=VARS[var]["beta"], rows=rows, joint_lattice=jl, joint_cont=jc)
        print("49A", var, "JOINT beta lattice %.4f cont %.4f pred %.2f" % (jl["beta"], jc["beta"], VARS[var]["beta"]), flush=True)
    return res


# ------------------------------------------------------------------ 49B
def radial_profile(var="N", q10=0.045):
    N, r, c, _ = solve(var, 10 * q10); rr = np.arange(0, c + 1, dtype=float)
    # radial profile from the axis (plus continuum beyond the box)
    return np.array([N[c + i, c, c] for i in range(c + 1)]), 10 * q10


def part_B(b=5.0, q10=0.045, L=1201, x0=-120.0, x1=120.0, dt=0.05):
    Nax, Q = radial_profile("N", q10); x = np.arange(L) - L // 2
    rr = np.sqrt(x.astype(float) ** 2 + b ** 2)
    Nw = np.where(rr <= 64, np.interp(rr, np.arange(65.0), Nax), VARS["N"]["Fi"](F_cont(rr, Q)))
    def rhs_fn(Nb, kind):
        Nl, Nr = np.roll(Nb, 1), np.roll(Nb, -1)
        if kind == "light":
            def acc(u):
                return Nb ** 2 * (Nb * Nl * (np.roll(u, 1) - u) + Nb * Nr * (np.roll(u, -1) - u))
        else:
            def acc(d):
                lin = Nb * Nl * (np.roll(d, 1) - d) + Nb * Nr * (np.roll(d, -1) - d)
                lin += d * (Nl * (Nl - Nb) + Nr * (Nr - Nb)) + Nb * (np.roll(d, 1) * (Nl - Nb) + np.roll(d, -1) * (Nr - Nb))
                return Nb ** 2 * lin
        return acc
    res = {}
    for k0 in (0.25, 0.5, 1.0):
        wdt = 8 / k0; row = {}
        for bg in ("flat", "well"):
            Nb = np.ones(L) if bg == "flat" else Nw
            N0 = Nb[np.argmin(np.abs(x - x0))]; om = N0 ** 2 * 2 * np.sin(k0 / 2)
            env = np.exp(-(x - x0) ** 2 / (2 * wdt ** 2)); u0 = env * np.cos(k0 * (x - x0)); v0 = om * env * np.sin(k0 * (x - x0))
            for kind in ("light", "dN"):
                acc = rhs_fn(Nb, kind); u, v = u0.copy(), v0.copy(); t = 0.0; prev = None; tcross = None
                while t < 600:
                    k1u, k1v = v, acc(u); k2u, k2v = v + dt / 2 * k1v, acc(u + dt / 2 * k1u)
                    k3u, k3v = v + dt / 2 * k2v, acc(u + dt / 2 * k2u); k4u, k4v = v + dt * k3v, acc(u + dt * k3u)
                    u = u + dt / 6 * (k1u + 2 * k2u + 2 * k3u + k4u); v = v + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v); t += dt
                    e = u ** 2 + (v / np.maximum(Nb ** 2 * 2 * np.sin(k0 / 2), 1e-9)) ** 2; cen = float(np.sum(e * x) / e.sum())
                    if prev is not None and prev < x1 <= cen:
                        tcross = t - dt * (cen - x1) / (cen - prev); break
                    prev = cen
                row[f"{bg}_{kind}"] = tcross
        row["delay_light"] = row["well_light"] - row["flat_light"]; row["delay_dN"] = row["well_dN"] - row["flat_dN"]
        row["delay_rel_diff"] = row["delay_dN"] / row["delay_light"] - 1; row["arrival_diff"] = row["well_dN"] - row["well_light"]
        res[str(k0)] = row
        print("49B k0", k0, {kk: (round(vv, 5) if vv is not None else None) for kk, vv in row.items()}, flush=True)
    ks = np.array([0.25, 0.5, 1.0]); ad = np.abs([res[str(q)]["arrival_diff"] for q in ks])
    slope = float(np.polyfit(np.log(ks), np.log(ad), 1)[0]) if np.all(ad > 0) else None
    return dict(b=b, q10=q10, rows=res, arrival_diff_logslope=slope)


def save(name, obj, t0):
    OUT.mkdir(parents=True, exist_ok=True)
    rec = dict(study="EX-GC49", part=name, label="EXPLORATORY", platform=platform.platform(), python=platform.python_version(), numpy=np.__version__,
               wall_s=time.time() - t0, result=obj)
    (OUT / f"{name}.json").write_text(json.dumps(rec, indent=1, default=float))


if __name__ == "__main__":
    for p in sys.argv[1:] or ["A", "B"]:
        t0 = time.time(); save(p, {"A": part_A, "B": part_B}[p](), t0); print(p, "done %.0fs" % (time.time() - t0), flush=True)
