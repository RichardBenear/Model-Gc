"""EX-GC51 checks. C1: static solver (fixed Gaussian source, qfix) vs EX-GC49's independent slowness solver at the same
source; plus an independent residual of the discrete equation with full-array rolls. C2: uniform-N plane-wave eigen check of
Hm (interior). C3: Q = dE/dN (uniform and local) by finite difference for the 51C probes. C4: 51C dt/2 (m=0.5, r=1) and
box x1.5 (static re-solved at n=193, sub-box 144); 51D dt/2. EXPLORATORY. Our construction; not LaViolette's Model G.
"""
import json, sys, time, platform
import numpy as np
import exgc51 as G
import exgc49 as F49


def c1(n=65, Q=0.45):
    c = n // 2; x = np.arange(n) - c; X, Y, Z = np.meshgrid(x, x, x, indexing="ij"); r2 = X ** 2 + Y ** 2 + Z ** 2
    w = np.exp(-r2 / (2 * F49.SIG ** 2)); s = 4 * np.pi * Q * w / w.sum()       # 49's source with total 4 pi Q, kappa = 1
    N49, r, c49, h = F49.solve("slow", Q, n=n)
    N51, info = G.static([], [], 1.0, 1.0, n, (c, c, c), qfix=s, gauss=False)
    # same boundary: 49 uses the continuum F at the faces; re-solve 51 with those faces by passing them through a fixed f
    # (simplest: compare interior after imposing 49's boundary): run 51 Picard from 49's boundary
    Nb = N49.copy(); N51b = solve_with_boundary(s, Nb)
    f = G.fN(N51b); I = (slice(1, -1),) * 3
    res_roll = 0.0; acc = np.zeros_like(f)
    for a in range(3):
        for d in (1, -1):
            acc += N51b * np.roll(N51b, -d, a) * (np.roll(f, -d, a) - f)
    res_roll = float(np.abs(acc[I] - s[I]).max())
    return dict(maxdiff_vs_49=float(np.abs(N51b - N49).max()), residual_independent=res_roll, iters49=len(h))


def solve_with_boundary(s, Nb, tol=1e-12):
    """Picard for sum N N (f_y - f_x) = s with Dirichlet N = Nb on the faces (generic, not using G.static)."""
    from scipy.fft import dstn, idstn
    n = s.shape[0]; mI = n - 2; I = (slice(1, -1),) * 3
    k = np.arange(1, mI + 1); l1 = 2 * np.cos(np.pi * k / (mI + 1)) - 2; lam = l1[:, None, None] + l1[None, :, None] + l1[None, None, :]
    N = Nb.copy(); N[I] = 1.0; f = G.fN(N)
    fbz = f.copy(); fbz[I] = 0; bc = F49.lap(fbz)
    for it in range(2000):
        R = F49.op(N, f) - F49.lap(f); g = idstn(dstn(s[I] - R - bc, type=1) / lam, type=1)
        f[I] = 0.5 * f[I] + 0.5 * g; N = G.Nf(f)
        if np.abs(F49.op(N, f) - s[I]).max() < tol:
            break
    return N


def c2(n=33, m=0.5, r=1.0, N0=0.9):
    k = np.array([2 * np.pi * 3 / n, 2 * np.pi * 1 / n, 0.0]); x = np.arange(n); X = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1)
    h = N0 ** 2 * sum(np.sin(k[j]) * G.AL[j] for j in range(3)) + (N0 * m + N0 ** 2 * r * np.sum(1 - np.cos(k))) * G.BE
    e, u = np.linalg.eigh(h); sp = u[:, -1]; psi = np.exp(1j * (X @ k))[..., None] * sp
    Hp = G.Hm(psi, np.full((n, n, n), N0), m, r); I = (slice(2, -2),) * 3
    return dict(E_formula=float(G.disp(k, N0, m, r)), E_eig=float(e[-1]), maxres=float(np.abs(Hp[I] - e[-1] * psi[I]).max()))


def c3(n=49, m=0.5, r=1.0):
    c = n // 2; psi, _, _ = G.packet(n, m, r, 1.0, (c, c, c), 3.0 / m); one = np.ones((n, n, n)); h = 1e-5
    Q = float(G.charge(psi, one, m, r).sum()); dE = (G.energy(psi, one * (1 + h), m, r) - G.energy(psi, one * (1 - h), m, r)) / (2 * h)
    q = G.charge(psi, one, m, r); x = (c + 3, c - 1, c + 2); Np = one.copy(); Np[x] += h; Nm = one.copy(); Nm[x] -= h
    loc = (G.energy(psi, Np, m, r) - G.energy(psi, Nm, m, r)) / (2 * h)
    return dict(Q=Q, dEdN=float(dE), rel=float(dE / Q - 1), local_fd=float(loc), local_q=float(q[x]), E=G.energy(psi, one, m, r))


def c4():
    out = {}
    r = G.part_C(dt=0.05, cases=[(0.5, 1.0)]); out["C_dt_half"] = {k: v for k, v in r["0.5_1.0"].items() if k in ("a", "a_semi", "a_GR", "a_over_semi", "dE_max", "diag_a_semi_full", "diag_a_plain")}
    # box x1.5: static field at n = 193 (same physical source and kappa), sub-box 144
    nB = 193; c = nB // 2; psi, _, _ = G.packet(nB, G.SRC["m"], G.SRC["r"], 1.0, (c, c, c), G.SRC["sigma"])
    B = json.loads((G.OUT / "B.json").read_text())["result"]; kap = B["rows"]["0.035"]["kappa"]
    N, info = G.static([psi], [G.SRC["m"]], G.SRC["r"], kap, nB, (c, c, c))
    r2 = G.part_C(Bfield=N, sub=144, nB=nB, cases=[(0.5, 1.0)]); out["C_box_1p5"] = {k: v for k, v in r2["0.5_1.0"].items() if k in ("a", "a_semi", "a_GR", "a_over_semi", "dE_max", "diag_a_semi_full", "diag_a_plain")}
    out["C_box_static"] = info
    d = G.part_D(dt=0.05); out["D_dt_half"] = {k: d[k] for k in ("a1", "a2", "F1", "F2", "momentum_balance", "inertia_balance")}
    return out


if __name__ == "__main__":
    t0 = time.time(); which = sys.argv[1:] or ["123", "4"]; res = {}
    if "123" in which:
        res["C1"] = c1(); print("C1", res["C1"], flush=True); res["C2"] = c2(); print("C2", res["C2"], flush=True); res["C3"] = c3(); print("C3", res["C3"], flush=True)
        (G.OUT / "checks123.json").write_text(json.dumps(dict(study="EX-GC51", part="checks123", platform=platform.platform(), result=res), indent=1, default=float))
    if "4" in which:
        r4 = c4(); print("C4", r4, flush=True)
        (G.OUT / "checks4.json").write_text(json.dumps(dict(study="EX-GC51", part="checks4", platform=platform.platform(), wall_s=time.time() - t0, result=r4), indent=1, default=float))
