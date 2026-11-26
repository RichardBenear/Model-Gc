"""EX-GC39C: scalar radiation from an orbiting pair if the lapse n = N - 1 is a free wave field,
n_tt - lap n = -4 pi G rho (static N = 1 - GM/r). Compared with the exact radiated power for the prescribed source
and with GR's quadrupole formula (predicted ratio 1/6). Plan: spec/EX-GC39C_plan.md.
Our construction, not LaViolette's Model G.
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.fft import dstn, idstn

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC39C"
G = 1.0


def sphere_points(n=1500):
    i = np.arange(n) + 0.5
    th = np.arccos(1 - 2 * i / n); ph = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)


def positions(t, masses, a_list, Om):
    return [np.array([s * a * np.cos(Om * t), s * a * np.sin(Om * t), 0.0]) for (a, s) in zip(a_list, (1, -1))]


def deposit(rho, c, t, masses, a_list, Om, sig, half=9, ramp=1.0, mean_sub=None):
    rho[:] = 0.0
    sl = slice(c - half, c + half + 1)
    ax = np.arange(-half, half + 1, dtype=float)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    sub = np.zeros_like(X)
    for M, p in zip(masses, positions(t, masses, a_list, Om)):
        sub += M * np.exp(-((X - p[0]) ** 2 + (Y - p[1]) ** 2 + (Z - p[2]) ** 2) / (2 * sig * sig)) / (2 * np.pi * sig * sig) ** 1.5
    if mean_sub is not None:          # remove the orbit-averaged (static, non-radiating) density
        sub = sub - mean_sub
    rho[sl, sl, sl] = ramp * sub


def lap(f, out):
    out[:] = 0.0
    out[1:-1, 1:-1, 1:-1] = (f[2:, 1:-1, 1:-1] + f[:-2, 1:-1, 1:-1] + f[1:-1, 2:, 1:-1] + f[1:-1, :-2, 1:-1]
                             + f[1:-1, 1:-1, 2:] + f[1:-1, 1:-1, :-2] - 6 * f[1:-1, 1:-1, 1:-1])
    return out


def run(masses, a_list, Om, sig=1.0, n=160, dt=0.25, sponge=20, t_on=150.0, periods=2, radii=(30, 40, 50), static=False):
    c = n // 2
    idx = np.arange(n)
    dist = np.minimum(idx, n - 1 - idx).astype(float)
    s1 = np.clip((sponge - dist) / sponge, 0, None) ** 2
    damp = 0.3 * np.maximum(np.maximum(s1[:, None, None], s1[None, :, None]), s1[None, None, :])
    n0 = np.zeros((n, n, n)); n1 = np.zeros_like(n0); L = np.zeros_like(n0); rho = np.zeros_like(n0)
    T_avg = periods * 2 * np.pi / (2 * Om) if not static else 0.0
    steps = int((t_on + T_avg) / dt)
    pts = sphere_points()
    samples = {(R, dR): [] for R in radii for dR in (-1.0, 0.0, 1.0)}; tsamp = []
    t = 0.0
    a = 1 / (1 + 0.5 * damp * dt); b = (1 - 0.5 * damp * dt)
    mean_sub = None
    if not static:
        # orbit-averaged density over 64 phases; only rho - <rho> is fed in (the static part does not radiate),
        # and it is switched on smoothly over t_ramp so that start-up transients are small
        acc = None
        for ph in range(64):
            deposit(rho, c, ph * 2 * np.pi / Om / 64, masses, a_list, Om, sig)
            s_ = rho[c - 9:c + 10, c - 9:c + 10, c - 9:c + 10].copy()
            acc = s_ if acc is None else acc + s_
        mean_sub = acc / 64
    t_ramp = 50.0
    for k in range(steps):
        ramp = 1.0 if static else float(np.sin(0.5 * np.pi * min(t / t_ramp, 1.0)) ** 2)
        deposit(rho, c, t, masses, a_list, Om, sig, ramp=ramp, mean_sub=mean_sub)
        lap(n1, L)
        n2 = a * (2 * n1 - b * n0 + dt * dt * (L - 4 * np.pi * G * rho))
        n2[0], n2[-1], n2[:, 0], n2[:, -1], n2[:, :, 0], n2[:, :, -1] = 0, 0, 0, 0, 0, 0
        if not static and t >= t_on and k % 2 == 0:
            for R in radii:
                for dR in (-1.0, 0.0, 1.0):
                    samples[(R, dR)].append(map_coordinates(n1, (c + (R + dR) * pts).T, order=1))
            tsamp.append(t)
        n0, n1 = n1, n2
        t += dt
    if static:
        return n1, rho
    # outgoing power from the radiated harmonics (s = 2, 4) of n on each sphere: P = (1/4piG) <ndot^2> 4 pi R^2.
    # A constant + linear trend absorbs the static field and slow transients. Exact for an outgoing spherical wave.
    ts = np.array(tsamp)
    cols = [np.ones_like(ts), ts - ts.mean()]
    for sh in (2, 4):
        cols += [np.cos(sh * Om * ts), np.sin(sh * Om * ts)]
    X = np.column_stack(cols)
    power, resid, power_ndot2 = {}, {}, {}
    for R in radii:
        co = {}
        for dR in (-1.0, 0.0, 1.0):
            Y = np.array(samples[(R, dR)])
            co[dR], *_ = np.linalg.lstsq(X, Y, rcond=None)
            if dR == 0.0:
                resid[R] = float(np.sqrt(np.mean((Y - X @ co[dR]) ** 2)) / np.sqrt(np.mean((X[:, 2:] @ co[dR][2:]) ** 2)))
        P = 0.0; P2 = 0.0
        for j, sh in enumerate((2, 4)):
            w = sh * Om
            a, b = co[0.0][2 + 2 * j], co[0.0][3 + 2 * j]
            da = (co[1.0][2 + 2 * j] - co[-1.0][2 + 2 * j]) / 2
            db = (co[1.0][3 + 2 * j] - co[-1.0][3 + 2 * j]) / 2
            # time-averaged flux -(1/4piG) <ndot d_r n> = -(1/4piG) (w/2)(b a' - a b')
            P += float(np.mean(-(w / 2) * (b * da - a * db))) / (4 * np.pi * G) * 4 * np.pi * R * R
            P2 += float(np.mean(w * w * (a * a + b * b) / 2)) / (4 * np.pi * G) * 4 * np.pi * R * R
        power[R] = P; power_ndot2[R] = P2
    resid = dict(fit_rel_residual=resid, power_ndot2_estimator=power_ndot2)
    return power, resid


def exact_power(masses, a_list, Om, sig=1.0, smax=8, nt=256):
    """P = sum_s (G w_s^2 / 8 pi) oint |rho_s(w_s n)|^2 dOmega, rho(t) = rho0 + sum Re[rho_s e^{-i s Om t}]."""
    pts = sphere_points(4000)
    ts = 2 * np.pi * np.arange(nt) / nt / Om
    P = 0.0; parts = {}
    for s in range(1, smax + 1):
        w = s * Om
        k = w * pts
        acc = np.zeros(len(pts), complex)
        for t in ts:
            for M, p in zip(masses, positions(t, masses, a_list, Om)):
                acc += M * np.exp(-1j * (k @ p)) * np.exp(1j * w * t)
        rs = 2 * acc / nt * np.exp(-0.5 * (w * sig) ** 2)
        Ps = G * w * w / (8 * np.pi) * float(np.mean(np.abs(rs) ** 2)) * 4 * np.pi
        parts[s] = Ps; P += Ps
    return P, parts


def quad_formulas(masses, a_list, Om):
    mu = masses[0] * masses[1] / (masses[0] + masses[1]); d = a_list[0] + a_list[1]
    return dict(P_scalar_quadrupole=(16 / 15) * G * mu ** 2 * d ** 4 * Om ** 6,
                P_GR=(32 / 5) * G * mu ** 2 * d ** 4 * Om ** 6,
                finite_size_factor=float(np.exp(-(2 * Om * 1.0) ** 2)))


def static_check():
    nfield, rho = run([1.0, 0.0], [0.0, 0.0], 0.0, static=True, t_on=200.0)
    n = nfield.shape[0]
    # exact Dirichlet lattice Poisson via DST-I on the interior
    inner = rho[1:-1, 1:-1, 1:-1]
    m = inner.shape[0]
    k = np.pi * np.arange(1, m + 1) / (m + 1)
    lam = 2 * np.cos(k) - 2
    Lsum = lam[:, None, None] + lam[None, :, None] + lam[None, None, :]
    sol = idstn(dstn(4 * np.pi * G * inner, type=1) / Lsum, type=1)
    c = n // 2
    rs = np.arange(10, 41, 5)
    out = []
    for r in rs:
        out.append(dict(r=int(r), wave=float(nfield[c + r, c, c]), poisson=float(sol[c + r - 1, c - 1, c - 1]),
                        minus_GM_over_r=-G / r))
    rel = max(abs(o["wave"] - o["poisson"]) / abs(o["poisson"]) for o in out)
    grad = [(sol[c + r, c - 1, c - 1] - sol[c + r - 2, c - 1, c - 1]) / 2 * r * r / G for r in (10, 20)]
    return dict(rows=out, max_rel_wave_vs_poisson=float(rel), r2_dn_dr_over_GM_at_10_20=[float(g) for g in grad])


if __name__ == "__main__":
    part = sys.argv[1]
    t0 = time.process_time()
    OUT.mkdir(parents=True, exist_ok=True)
    Om = 0.08
    if part == "static":
        res = static_check()
    else:
        if part == "equal":
            masses, a_list = [1.0, 1.0], [2.5, 2.5]
        else:   # unequal: M1 = 2 M2, a1 M1 = a2 M2, same separation 5
            masses, a_list = [2.0, 1.0], [5 / 3, 10 / 3]
        fl, sd = run(masses, a_list, Om)
        Pex, parts = exact_power(masses, a_list, Om)
        q = quad_formulas(masses, a_list, Om)
        res = dict(masses=masses, radii_orbit=a_list, Omega=Om, power=fl, fit_rel_residual=sd, P_exact=Pex, P_exact_harmonics=parts, **q,
                   lattice_over_exact={R: v / Pex for R, v in fl.items()},
                   exact_over_scalar_quadrupole=Pex / q["P_scalar_quadrupole"],
                   lattice_over_GR={R: v / q["P_GR"] for R, v in fl.items()})
    res["cpu_seconds"] = time.process_time() - t0
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps(res, indent=1, default=float)[:2500])
