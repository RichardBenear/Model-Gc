"""EX-GC51: one declared model end to end. First-quantised Wilson-Dirac matter with coincidence weights (bonds N_vN_w,
mass on sites N_v); N carried by coincidences as the light slowness f = -1/(2N^2), sourced by q = dH_m/dN; light by
the 37B rule. 51A field signal speed vs light; 51B static G_eff, beta; 51C probe fall; 51D two-body reciprocity;
51E Shapiro ratio. Plan: spec/EX-GC51_plan.md (amendment A1). EXPLORATORY. Our construction; not LaViolette's Model G.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np
from scipy.fft import dstn, idstn

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc48 as D48
import exgc49 as F49

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC51"
GAM = D48.GAM; AL = GAM[1:4]; BE = GAM[4]
AL2I = [a / 2j for a in AL]
fN = lambda N: -0.5 / N ** 2          # slowness transport variable
Nf = lambda f: np.sqrt(-0.5 / f)


def shp(a, ax, d):
    """value at x + d e_ax, zero beyond the box."""
    out = np.roll(a, -d, axis=ax); sl = [slice(None)] * a.ndim
    sl[ax] = slice(-1, None) if d == 1 else slice(0, 1); out[tuple(sl)] = 0; return out


def mv(A, psi):
    return psi @ A.T


def Hm(psi, N, m, r):
    out = m * N[..., None] * mv(BE, psi)
    for a in range(3):
        pw = shp(psi, a, 1); Wf = N * shp(N, a, 1)
        out += Wf[..., None] * (mv(AL2I[a], pw) + (r / 2) * mv(BE, psi - pw))
        pv = shp(psi, a, -1); Wb = N * shp(N, a, -1)
        out += Wb[..., None] * (-mv(AL2I[a], pv) + (r / 2) * mv(BE, psi - pv))
    return out


def charge(psi, N, m, r):
    """q_x = d<H>/dN_x (bond energies weighted by the neighbour's N)."""
    rb = np.real(np.sum(psi.conj() * mv(BE, psi), -1)); q = m * rb
    for a in range(3):
        pw = shp(psi, a, 1)
        e = 2 * np.real(np.sum(psi.conj() * mv(AL2I[a], pw), -1)) - r * np.real(np.sum(psi.conj() * mv(BE, pw), -1)) + (r / 2) * (rb + shp(rb, a, 1))
        q += shp(N, a, 1) * e                 # bond (x, x+a) contributes N_{x+a} e to q_x
        q += np.roll(N, 1, axis=a) * np.roll(e, 1, axis=a) * (np.arange(N.shape[a]) > 0).reshape([-1 if i == a else 1 for i in range(3)])  # bond (x-a, x): N_{x-a} e(x-a)
    return q


def energy(psi, N, m, r):
    return float(np.real(np.vdot(psi, Hm(psi, N, m, r))))


def disp(k, N, m, r):
    S2 = np.sum(np.sin(k) ** 2, -1); W = np.sum(1 - np.cos(k), -1)
    return np.sqrt(N ** 4 * S2 + (N * m + N ** 2 * r * W) ** 2)


def packet(n, m, r, N0, x0, sigma, chi=(1, 0, 0, 0)):
    """positive-energy (w.r.t. uniform N0) Gaussian packet at rest centred at x0 (index coords), zero on the boundary."""
    k = 2 * np.pi * np.fft.fftfreq(n); K = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
    S = np.sin(K); W = np.sum(1 - np.cos(K), -1)
    h = N0 ** 2 * np.einsum("...j,jab->...ab", S, np.array(AL)) + (N0 * m + N0 ** 2 * r * W)[..., None, None] * BE
    E = disp(K, N0, m, r); P = 0.5 * (np.eye(4) + h / E[..., None, None])
    g = np.exp(-np.sum(K ** 2, -1) * sigma ** 2 / 2) * np.exp(-1j * (K @ np.asarray(x0, float)))
    a = P @ np.asarray(chi, complex); a = a / np.linalg.norm(a, axis=-1, keepdims=True) * g[..., None]
    psi = np.fft.ifftn(a, axes=(0, 1, 2))
    for ax in range(3):
        sl = [slice(None)] * 4; sl[ax] = 0; psi[tuple(sl)] = 0; sl[ax] = -1; psi[tuple(sl)] = 0
    w = np.abs(g) ** 2; w = w / w.sum()
    return psi / np.linalg.norm(psi), K, w


# ------------------------------------------------------------------ static field for frozen matter
def static(psis, ms, r, kappa, n, centre, tol=1e-11, maxit=300, gauss=True, qfix=None):
    c = np.array(centre, float); x = np.arange(n); X = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1)
    rr = np.sqrt(np.sum((X - c) ** 2, -1)); rr = np.where(rr == 0, 1e-9, rr)
    N = np.ones((n, n, n)); mI = n - 2
    k = np.arange(1, mI + 1); l1 = 2 * np.cos(np.pi * k / (mI + 1)) - 2; lam = l1[:, None, None] + l1[None, :, None] + l1[None, None, :]
    I = (slice(1, -1),) * 3; Fh, Gh = [], []
    for it in range(maxit):
        q = qfix if qfix is not None else sum(charge(p, N, mm, r) for p, mm in zip(psis, ms)); Q = q.sum()
        Nb = np.exp(-kappa * Q / (4 * np.pi * rr)) if gauss else np.ones_like(rr)
        f = fN(N); fb = fN(Nb); f[0], f[-1], f[:, 0], f[:, -1], f[:, :, 0], f[:, :, -1] = fb[0], fb[-1], fb[:, 0], fb[:, -1], fb[:, :, 0], fb[:, :, -1]
        N = Nf(f)
        res = F49.op(N, f) - kappa * q[I]; rn = float(np.abs(res).max())
        if rn < tol:
            break
        fbz = f.copy(); fbz[I] = 0; bc = F49.lap(fbz)
        R = F49.op(N, f) - F49.lap(f); g = idstn(dstn(kappa * q[I] - R - bc, type=1) / lam, type=1)
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


def dirs_samples(N, c, rmin, rmax):
    out = {}
    for name, d in (("axis", (1, 0, 0)), ("face", (1, 1, 0)), ("body", (1, 1, 1))):
        d = np.array(d); L = np.linalg.norm(d); ss = [s for s in range(1, 200) if rmin <= s * L <= rmax and c + s < N.shape[0] - 1]
        out[name] = (np.array([s * L for s in ss]), np.array([N[c + s * d[0], c + s * d[1], c + s * d[2]] for s in ss]))
    return out


# ------------------------------------------------------------------ 51B
SRC = dict(m=0.5, sigma=3.0, r=1.0)
U12 = (0.01, 0.02, 0.035, 0.05)


def part_B(n=129):
    c = n // 2; psi, _, _ = packet(n, SRC["m"], SRC["r"], 1.0, (c, c, c), SRC["sigma"])
    Q0 = float(charge(psi, np.ones((n, n, n)), SRC["m"], SRC["r"]).sum()); rows = {}; fields = {}
    for u in U12:
        kap = 4 * np.pi * 12 * u / Q0; t0 = time.time()
        N, info = static([psi], [SRC["m"]], SRC["r"], kap, n, (c, c, c)); fields[u] = N
        S = dirs_samples(N, c, 12, 48); fits = {k: F49.fit_one(*v) for k, v in S.items()}
        rows[str(u)] = dict(kappa=kap, Q0=Q0, info=info, fits=fits, Geff={k: v["A"] * 4 * np.pi / (kap * info["Q"]) for k, v in fits.items()},
                            N_centre=float(N[c, c, c]), wall=time.time() - t0)
        print("51B U12", u, "kappa %.4f Q %.6f (Q0 %.6f) it %d res %.1e | Geff %s | beta raw %s" % (kap, info["Q"], Q0, info["iters"], info["resid"],
              {k: round(v, 5) for k, v in rows[str(u)]["Geff"].items()}, {k: round(v["beta"], 4) for k, v in fits.items()}), flush=True)
    sep = {}
    for d in ("axis", "face", "body"):
        ks = np.array([rows[str(u)]["kappa"] for u in U12]); a1 = np.array([rows[str(u)]["fits"][d]["coef"][0] for u in U12]); a2 = np.array([rows[str(u)]["fits"][d]["coef"][1] for u in U12])
        cc = np.linalg.lstsq(np.stack([ks, ks ** 2, ks ** 3], 1), a2, rcond=None)[0]; k1 = np.mean(a1 / ks); sep[d] = dict(beta=float(2 * cc[1] / k1 ** 2), lin=float(cc[0]))
    print("51B beta separated", {k: round(v["beta"], 4) for k, v in sep.items()}, flush=True)
    np.save(OUT / "partial_field_B.npy", np.stack([fields[u] for u in U12]))
    return dict(rows=rows, beta_sep=sep, n=n)


# ------------------------------------------------------------------ 51A
def part_A(n=97, u=0.02, T=60.0, dt=0.1):
    c = n // 2; psi, _, _ = packet(n, SRC["m"], SRC["r"], 1.0, (c, c, c), SRC["sigma"])
    Q0 = float(charge(psi, np.ones((n, n, n)), SRC["m"], SRC["r"]).sum()); kap = 4 * np.pi * 12 * u / Q0
    Ns, info = static([psi], [SRC["m"]], SRC["r"], kap, n, (c, c, c), gauss=False)   # same flat Dirichlet boundary as the dynamic run
    I = (slice(1, -1),) * 3
    pts = {f"axis{s}": (c + s, c, c) for s in range(10, 41, 5)}; pts.update({f"body{s}": (c + int(round(s / np.sqrt(3))),) * 3 for s in range(10, 41, 5)})
    def run(kind):
        f = fN(np.ones((n, n, n))); g = np.zeros_like(f); u_ = np.zeros_like(f); v_ = np.zeros_like(f); rec = {k: [] for k in pts}; ts = []
        q1 = charge(psi, np.ones((n, n, n)), SRC["m"], SRC["r"])
        def rhs(f, g):
            N = Nf(f); q = charge(psi, N, SRC["m"], SRC["r"]) if kind == "field" else q1
            if kind == "field":
                df = N * g; dg = np.zeros_like(g); dg[I] = N[I] * (F49.op(N, f) - kap * q[I])
            else:
                one = np.ones_like(f); df = g.copy(); dg = np.zeros_like(g); dg[I] = F49.lap(f) - kap * q[I]
            df[0] = df[-1] = 0; df[:, 0] = df[:, -1] = 0; df[:, :, 0] = df[:, :, -1] = 0
            return df, dg
        t = 0.0
        if kind == "light":
            f = np.zeros_like(f)
        while t < T - 1e-9:
            k1 = rhs(f, g); k2 = rhs(f + dt / 2 * k1[0], g + dt / 2 * k1[1]); k3 = rhs(f + dt / 2 * k2[0], g + dt / 2 * k2[1]); k4 = rhs(f + dt * k3[0], g + dt * k3[1])
            f = f + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); g = g + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]); t += dt; ts.append(t)
            for k_, p in pts.items():
                rec[k_].append(float(Nf(f[p]) - 1) if kind == "field" else float(f[p]))
        return np.array(ts), rec
    out = {}
    # light (linear, N = 1): its static target is the linear Poisson solution; the field's target is Ns
    Nl, _ = None, None
    for kind in ("field", "light"):
        ts, rec = run(kind); rows = {}
        if kind == "light":
            # static linear target: solve lap f = kappa q with zero boundary by DST
            mI = n - 2; kk = np.arange(1, mI + 1); l1 = 2 * np.cos(np.pi * kk / (mI + 1)) - 2; lam = l1[:, None, None] + l1[None, :, None] + l1[None, None, :]
            q1 = charge(psi, np.ones((n, n, n)), SRC["m"], SRC["r"]); fl = np.zeros((n, n, n)); fl[I] = idstn(dstn(kap * q1[I], type=1) / lam, type=1)
        for k_, p in pts.items():
            target = (Ns[p] - 1) if kind == "field" else fl[p]; y = np.array(rec[k_]); i = np.argmax(y / target >= 0.5)
            th = ts[i - 1] + (0.5 * target - y[i - 1]) / (y[i] - y[i - 1]) * (ts[i] - ts[i - 1]) if i > 0 else None
            rr = float(np.sqrt(sum((pp - c) ** 2 for pp in p))); rows[k_] = dict(r=rr, t_half=th, target=float(target))
        rs = np.array([v["r"] for v in rows.values()]); th = np.array([v["t_half"] for v in rows.values()])
        slope = float(np.polyfit(rs, th, 1)[0]); out[kind] = dict(rows=rows, slope=slope, c=1 / slope)
        print("51A", kind, "c = %.4f" % (1 / slope), flush=True)
    out["c_ratio"] = out["field"]["c"] / out["light"]["c"]; out["kappa"] = kap
    return out


# ------------------------------------------------------------------ evolution helpers (51C, 51D)
def evolve(psis, ms, r, f, kappa, T, dt, dynamic, record):
    n = f.shape[0]; I = (slice(1, -1),) * 3; g = np.zeros_like(f)
    def rhs(ps, f, g):
        N = Nf(f); dps = [-1j * Hm(p, N, mm, r) for p, mm in zip(ps, ms)]
        if not dynamic:
            return dps, None, None
        q = sum(charge(p, N, mm, r) for p, mm in zip(ps, ms)); df = N * g; dg = np.zeros_like(g); dg[I] = N[I] * (F49.op(N, f) - kappa * q[I])
        for a in range(3):
            sl = [slice(None)] * 3; sl[a] = 0; df[tuple(sl)] = 0; sl[a] = -1; df[tuple(sl)] = 0
        return dps, df, dg
    t = 0.0; hist = [record(psis, f, 0.0)]
    while t < T - 1e-9:
        def add(ps, f, g, k, h):
            return [p + h * kp for p, kp in zip(ps, k[0])], (f + h * k[1] if dynamic else f), (g + h * k[2] if dynamic else g)
        k1 = rhs(psis, f, g); k2 = rhs(*add(psis, f, g, k1, dt / 2)); k3 = rhs(*add(psis, f, g, k2, dt / 2)); k4 = rhs(*add(psis, f, g, k3, dt))
        psis = [p + dt / 6 * (a + 2 * b + 2 * c + d) for p, a, b, c, d in zip(psis, k1[0], k2[0], k3[0], k4[0])]
        if dynamic:
            f = f + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]); g = g + dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        t += dt; hist.append(record(psis, f, t))
    return psis, f, hist


def centroid(p):
    n = p.shape[0]; rho = np.sum(np.abs(p) ** 2, -1); x = np.arange(n)
    return float(np.sum(rho * x[:, None, None]) / rho.sum())


def momentum_x(p):
    n = p.shape[0]; pk = np.fft.fftn(p, axes=(0, 1, 2)); w = np.sum(np.abs(pk) ** 2, -1); k = 2 * np.pi * np.fft.fftfreq(n)
    return float(np.sum(w * k[:, None, None]) / w.sum())


def fit_acc(ts, dx):
    ts = np.asarray(ts); A = np.stack([ts ** 2 / 2, ts ** 3, ts ** 4], 1); return float(np.linalg.lstsq(A, dx, rcond=None)[0][0])


def semi(K, w, N0, m, r, grad):
    h = 1e-4; e = np.array([h, 0, 0])
    Ekk = (disp(K + e, N0, m, r) + disp(K - e, N0, m, r) - 2 * disp(K, N0, m, r)) / h ** 2
    EN = (disp(K, N0 + h, m, r) - disp(K, N0 - h, m, r)) / (2 * h)
    return float(-np.sum(w * Ekk * EN) * grad), float(1 / np.sum(w * Ekk))


def semi_vterm(K, w, N0, m, r, grad):
    """POST-RUN DIAGNOSTIC (not registered): the term the A1 comparator omits, <E_k E_kN> dN/dx, from d(E_k)/dt = E_kk kdot + E_kN xdot dN/dx."""
    h = 1e-4; e = np.array([h, 0, 0]); d = lambda k, N: disp(k, N, m, r)
    Ek = (d(K + e, N0) - d(K - e, N0)) / (2 * h)
    EkN = ((d(K + e, N0 + h) - d(K - e, N0 + h)) - (d(K + e, N0 - h) - d(K - e, N0 - h))) / (4 * h * h)
    return float(np.sum(w * Ek * EkN) * grad)


# ------------------------------------------------------------------ 51C
def part_C(Bfield=None, u=0.035, sub=96, d0=16, T=12.0, dt=0.1, cases=None, nB=129):
    if Bfield is None:
        Bfield = np.load(OUT / "partial_field_B.npy")[U12.index(u)]
    c = nB // 2; lo = c + d0 - sub // 2; Nsub = Bfield[lo:lo + sub, c - sub // 2:c + sub // 2, c - sub // 2:c + sub // 2]
    ci = (c + d0 - lo, sub // 2, sub // 2); N0 = float(Nsub[ci])
    cases = cases or [(m, r) for m in (0.25, 0.5) for r in (1.0, 0.5, 0.0)]; res = {}
    for m, r in cases:
        sigma = 3.0 / m; psi, K, w = packet(sub, m, r, N0, ci, sigma)
        rho = np.sum(np.abs(psi) ** 2, -1); grad = float(np.sum(rho * (np.roll(Nsub, -1, 0) - np.roll(Nsub, 1, 0)) / 2))
        Nloc = float(np.sum(rho * Nsub)); aGR = -Nloc ** 3 * grad; asemi, mi = semi(K, w, Nloc, m, r, grad)
        t0 = time.time(); runs = {}
        for lab, Nf_ in (("field", Nsub), ("flat", np.full_like(Nsub, N0))):
            E0 = energy(psi, Nf_, m, r)
            _, _, hist = evolve([psi], [m], r, fN(Nf_), 0.0, T, dt, False, lambda ps, f, t: (t, centroid(ps[0]), energy(ps[0], Nf_, m, r)))
            runs[lab] = dict(t=[h[0] for h in hist], x=[h[1] for h in hist], dE_max=float(max(abs(h[2] - E0) for h in hist)))
        dx = np.array(runs["field"]["x"]) - np.array(runs["flat"]["x"]); a = fit_acc(runs["field"]["t"], dx)
        av = semi_vterm(K, w, Nloc, m, r, grad); tt = np.asarray(runs["field"]["t"]); late = tt >= T / 2
        a_plain = float(np.mean(dx[late] / (tt[late] ** 2 / 2)))   # diagnostic: mean of 2 dx/t^2 over the second half
        res[f"{m}_{r}"] = dict(m=m, r=r, sigma=sigma, N0=N0, Nloc=Nloc, grad=grad, a=a, a_GR=aGR, a_semi=asemi, m_i=mi,
                               diag_a_semi_full=asemi + av, diag_a_plain=a_plain,
                               a_over_GR=a / aGR, a_over_semi=a / asemi, pred_inf=1 + r * m, dE_max={k: v["dE_max"] for k, v in runs.items()},
                               t=runs["field"]["t"], dx=dx.tolist(), wall=time.time() - t0)
        print("51C m %.2f r %.1f: a/aGR %.5f a/asemi %.5f semi/GR %.5f (1+rm %.3f) | diag a/semi_full %.5f plain/semi_full %.5f | dEmax %.1e %.0fs" % (m, r, a / aGR, a / asemi, asemi / aGR, 1 + r * m, a / (asemi + av), a_plain / (asemi + av),
              runs["field"]["dE_max"], time.time() - t0), flush=True)
    for r in (1.0, 0.5, 0.0):
        if f"0.25_{r}" in res and f"0.5_{r}" in res:
            A1, A2 = res[f"0.25_{r}"], res[f"0.5_{r}"]
            res[f"ratio_r{r}"] = dict(meas=A1["a"] / A2["a"], semi=A1["a_semi"] / A2["a_semi"], rm=(1 + r * 0.25) / (1 + r * 0.5))
    return res


# ------------------------------------------------------------------ 51D
def part_D(n=97, sep=20, ms=(0.4, 0.6), r=1.0, u=0.01, T=20.0, dt=0.1):
    c = n // 2; x1 = (c - sep // 2, c, c); x2 = (c + sep // 2, c, c)
    p1, K, w1 = packet(n, ms[0], r, 1.0, x1, 3.0 / ms[0]); p2, _, w2 = packet(n, ms[1], r, 1.0, x2, 3.0 / ms[1])
    one = np.ones((n, n, n)); Qs = [float(charge(p, one, mm, r).sum()) for p, mm in ((p1, ms[0]), (p2, ms[1]))]
    kap = 4 * np.pi * sep * u / max(Qs); res = dict(kappa=kap, Q=Qs)
    def rec(ps, f, t):
        return (t, [centroid(p) for p in ps], [momentum_x(p) for p in ps])
    runs = {}
    for lab, ps, mm in (("pair", [p1, p2], list(ms)), ("one", [p1], [ms[0]]), ("two", [p2], [ms[1]])):
        t0 = time.time(); N, info = static(ps, mm, r, kap, n, (c, c, c))
        _, _, hist = evolve(ps, mm, r, fN(N), kap, T, dt, True, rec); runs[lab] = hist; print("51D", lab, "%.0fs" % (time.time() - t0), flush=True)
    t = np.array([h[0] for h in runs["pair"]])
    x1d = np.array([h[1][0] for h in runs["pair"]]) - np.array([h[1][0] for h in runs["one"]]); x2d = np.array([h[1][1] for h in runs["pair"]]) - np.array([h[1][0] for h in runs["two"]])
    P1 = np.array([h[2][0] for h in runs["pair"]]) - np.array([h[2][0] for h in runs["one"]]); P2 = np.array([h[2][1] for h in runs["pair"]]) - np.array([h[2][0] for h in runs["two"]])
    a1, a2 = fit_acc(t, x1d), fit_acc(t, x2d)
    F1 = float(np.polyfit(t, P1, 2)[1]); F2 = float(np.polyfit(t, P2, 2)[1])         # dP/dt at t = 0 from a quadratic fit
    _, mi1 = semi(K, w1, 1.0, ms[0], r, 0.0); _, mi2 = semi(K, w2, 1.0, ms[1], r, 0.0)
    res.update(a1=a1, a2=a2, F1=F1, F2=F2, mi1=mi1, mi2=mi2, momentum_balance=(F1 + F2) / abs(F1), inertia_balance=(mi1 * a1 + mi2 * a2) / abs(mi1 * a1),
               t=t.tolist(), x1d=x1d.tolist(), x2d=x2d.tolist(), P1=P1.tolist(), P2=P2.tolist())
    print("51D a1 %.4e a2 %.4e | F1 %.4e F2 %.4e | (F1+F2)/|F1| %.4f | (mi a) balance %.4f" % (a1, a2, F1, F2, res["momentum_balance"], res["inertia_balance"]), flush=True)
    return res


# ------------------------------------------------------------------ 51E
def part_E(Bfield=None, u=0.05, nB=129, bs=(8, 16, 24)):
    if Bfield is None:
        Bfield = np.load(OUT / "partial_field_B.npy")[U12.index(u)]
    c = nB // 2; res = {}
    for b in bs:
        Nl = Bfield[:, c + b, c][1:-1]
        two = float(np.sum(1 / Nl ** 2 - 1)); one = float(np.sum(1 / Nl - 1))
        res[str(b)] = dict(eikonal_two=two, eikonal_one=one, ratio_eikonal=two / one)
        # 1D chain waves: coincidence (speed N^2) vs one-site control (speed N), as 49B
        x = np.arange(len(Nl)) - len(Nl) // 2; k0, wd, x0, x1, dt = 1.0, 6.0, -45.0, 45.0, 0.05
        def run(Nb, kind):
            Nl_, Nr_ = np.roll(Nb, 1), np.roll(Nb, -1)
            if kind == "two":
                acc = lambda q: Nb ** 2 * (Nb * Nl_ * (np.roll(q, 1) - q) + Nb * Nr_ * (np.roll(q, -1) - q))
            else:
                wl = 2 * Nb * Nl_ / (Nb + Nl_) / Nb; wr = 2 * Nb * Nr_ / (Nb + Nr_) / Nb
                acc = lambda q: Nb ** 2 * (wl * (np.roll(q, 1) - q) + wr * (np.roll(q, -1) - q))
            i0 = np.argmin(np.abs(x - x0)); N0 = Nb[i0]; sp = N0 ** 2 if kind == "two" else N0
            om = sp * 2 * np.sin(k0 / 2); env = np.exp(-(x - x0) ** 2 / (2 * wd ** 2)); q = env * np.cos(k0 * (x - x0)); v = om * env * np.sin(k0 * (x - x0))
            t = 0.0; prev = None
            while t < 400:
                k1q, k1v = v, acc(q); k2q, k2v = v + dt / 2 * k1v, acc(q + dt / 2 * k1q); k3q, k3v = v + dt / 2 * k2v, acc(q + dt / 2 * k2q); k4q, k4v = v + dt * k3v, acc(q + dt * k3q)
                q = q + dt / 6 * (k1q + 2 * k2q + 2 * k3q + k4q); v = v + dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v); t += dt
                e = q ** 2 + (v / (sp * 2 * np.sin(k0 / 2))) ** 2; cen = float(np.sum(e * x) / e.sum())
                if prev is not None and prev < x1 <= cen:
                    return t - dt * (cen - x1) / (cen - prev)
                prev = cen
            return None
        flat = np.ones_like(Nl); d2 = run(Nl, "two") - run(flat, "two"); d1 = run(Nl, "one") - run(flat, "one")
        res[str(b)].update(wave_two=d2, wave_one=d1, ratio_wave=d2 / d1)
        print("51E b", b, res[str(b)], flush=True)
    return res


def save(name, obj, t0):
    OUT.mkdir(parents=True, exist_ok=True)
    rec = dict(study="EX-GC51", part=name, label="EXPLORATORY", platform=platform.platform(), python=platform.python_version(), numpy=np.__version__,
               wall_s=time.time() - t0, result=obj)
    (OUT / f"{name}.json").write_text(json.dumps(rec, indent=1, default=float))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for p in sys.argv[1:]:
        if p.startswith("Cm"):   # one probe mass, all three r (lets the two masses run in parallel); merged by exgc51_recompute.py merge
            t0 = time.time(); mm = float(p[2:]); save("C_m" + p[2:], part_C(cases=[(mm, r) for r in (1.0, 0.5, 0.0)]), t0); print(p, "done", flush=True); continue
        t0 = time.time(); save(p, {"A": part_A, "B": part_B, "C": part_C, "D": part_D, "E": part_E}[p](), t0); print(p, "done %.0fs" % (time.time() - t0), flush=True)
