"""EX-GC57: frame dragging read by light, by orbiting matter and by spin; and 56's light-drift readout closed.
57A exact Fourier light drift (one-way eigenpacket vs 56's launch; energy vs envelope centroid; |u|^2 vs mode-energy predictor).
57B 53's rotating ring (momentum rule, frozen static fields): local Omega = -1/2 curl beta read by a synchronised light loop, an orbiting
probe ring (energy splitting, torque) and a resting Dirac spin, with a pure-gradient control.
Plan: spec/EX-GC57_plan.md. EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory.

Implementation choices fixed before the run (not stated in the plan):
- Probes (ring and spin packets) are built in flat N in a 48^3 box centred on the probe position; the global shift is cut to that box.
- Omega_light = (dt_wave / readout) * Omega_loop, Omega_loop = -sum(bl) / (2 A): the wave/readout ratio times the coordinate loop curl,
  so the clock's 1/N^4 conversion to local time is removed and light is compared with the coordinate curl like the matter probes.
- The point curl uses the 4 plaquettes normal to the component that share the probe site (the plan said 8 cells; 4 is what touches the site).
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G
import exgc52 as S
import exgc53 as T53
import exgc55_posthoc2 as SY

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC57"
UP, DOWN = T53.UP, T53.DOWN; M = 0.6; NB = 48


# ------------------------------------------------------------------ 57A
def light_exact(k0, width, launch, beta, L=512, T=100.0):
    s = np.arange(L); d = (s - L // 2 + L) % L - L // 2; env = np.exp(-d ** 2 / (2 * width ** 2)); k = 2 * np.pi * np.fft.fftfreq(L)
    Sk = np.sin(k); Om = 2 * np.abs(np.sin(k / 2)); O0 = 2 * np.sin(k / 2)
    grad = lambda a: 0.5 * (np.roll(a, -1) - np.roll(a, 1))
    if launch == "a":                                                            # exact right-moving eigenpacket
        A = np.fft.fft(env * np.exp(1j * k0 * d)) * (k > 0)
        u = 2 * np.real(np.fft.ifft(A)); w = 2 * np.real(np.fft.ifft(-1j * O0 * A))
    else:                                                                        # 56's launch
        u = env * np.cos(k0 * d); ut = 2 * np.sin(k0 / 2) * env * np.sin(k0 * d); w = ut - beta * grad(u)
    uh, wh = np.fft.fft(u), np.fft.fft(w); ph = np.exp(1j * beta * Sk * T); c, sn = np.cos(Om * T), np.sin(Om * T)
    sinc = np.where(Om > 0, sn / np.where(Om > 0, Om, 1), T)
    uT = np.real(np.fft.ifft(ph * (c * uh + sinc * wh))); wT = np.real(np.fft.ifft(ph * (-Om * sn * uh + c * wh)))
    ang = 2 * np.pi * s / L
    def cen_E(u, w):
        gl = (np.roll(u, -1) - u) ** 2; e = 0.5 * w ** 2 + 0.25 * (gl + np.roll(gl, 1)) + beta * w * grad(u)
        return np.angle(np.sum(e * np.exp(1j * ang))) * L / (2 * np.pi)
    def cen_env(u):
        z = np.fft.ifft(2 * np.fft.fft(u) * (k > 0)); return np.angle(np.sum(np.abs(z) ** 2 * np.exp(1j * ang))) * L / (2 * np.pi)
    wrap = lambda x: (x + L / 2) % L - L / 2
    return wrap(cen_E(uT, wT) - cen_E(u, w)), wrap(cen_env(uT) - cen_env(u))


def predictors(k0, width, beta, L=512):
    s = np.arange(L); d = (s - L // 2 + L) % L - L // 2; env = np.exp(-d ** 2 / (2 * width ** 2)); k = 2 * np.pi * np.fft.fftfreq(L)
    full = np.abs(np.fft.fft(env * np.exp(1j * k0 * d))) ** 2; A2 = full * (k > 0); O0 = 2 * np.sin(k / 2)
    p_u2 = -float(np.sum(A2 * np.cos(k)) / A2.sum()); p_56 = -float(np.sum(full * np.cos(k)) / full.sum())
    eps = lambda b: O0 ** 2 * A2 - b * O0 * np.sin(k) * A2; vg = lambda b: np.cos(k / 2) - b * np.cos(k)
    p_eps = float((np.sum(eps(beta) * vg(beta)) / eps(beta).sum() - np.sum(eps(0) * vg(0)) / eps(0).sum()) / beta)
    p_eps0 = -float(np.sum(eps(0) * np.cos(k)) / eps(0).sum())
    return dict(u2=p_u2, pred56=p_56, eps=p_eps, eps0=p_eps0, point=-float(np.cos(k0)))


def part_A(beta=1e-3, T=100.0):
    rows = []
    for k0 in (0.25, 0.5):
        for width in (8.0, 16.0, 32.0):
            pr = predictors(k0, width, beta); row = dict(k0=k0, width=width, pred=pr)
            for la in ("a", "b"):
                eb, nb = light_exact(k0, width, la, beta, T=T); e0, n0 = light_exact(k0, width, la, 0.0, T=T)
                dE = (eb - e0) / (beta * T); dN = (nb - n0) / (beta * T)
                row[la] = dict(drift_E=dE, drift_env=dN, E_over_eps=dE / pr["eps"], E_over_u2=dE / pr["u2"], env_over_u2=dN / pr["u2"], E_over_56pred=dE / pr["pred56"])
            rows.append(row); print("57A k0 %.2f w %2.0f  a: E/eps %.6f env/u2 %.6f | b: E/pred56 %.6f" % (k0, width, row["a"]["E_over_eps"], row["a"]["env_over_u2"], row["b"]["E_over_56pred"]), flush=True)
    z = light_exact(0.25, 16.0, "a", 0.0, T=T); z2 = light_exact(0.25, 16.0, "a", 0.0, T=T); c1 = dict(zero_drift=float(abs(z[0] - z2[0]) + abs(z[1] - z2[1])))
    import exgc56 as X56
    rk = X56.light_drift(0.25); fb = [r for r in rows if r["k0"] == 0.25 and r["width"] == 16.0][0]["b"]["drift_E"]
    c1["rk4_drift"] = rk["drift"]; c1["fourier_drift"] = fb; c1["rel"] = abs(fb / rk["drift"] - 1); print("57A C1", c1, flush=True)
    return dict(rows=rows, C1=c1)


# ------------------------------------------------------------------ 57B helpers
def curl_site(B, x):
    """(curl beta) at site x from the 4 plaquettes normal to each component that share the site."""
    out = np.zeros(3)
    for comp in range(3):
        a, b = [i for i in range(3) if i != comp]
        if (a, b) != ((comp + 1) % 3, (comp + 2) % 3):
            a, b = b, a                                                          # right-handed (a, b, comp)
        vals = []
        for da in (-1, 0):
            for db in (-1, 0):
                p = list(x); p[a] += da; p[b] += db; pa = list(p); pa[a] += 1; pb = list(p); pb[b] += 1
                vals.append(B[a][tuple(p)] + B[b][tuple(pa)] - B[a][tuple(pb)] - B[b][tuple(p)])
        out[comp] = np.mean(vals)
    return out


def sub(B, x0):
    h = NB // 2; sl = tuple(slice(x0[i] - h, x0[i] + h) for i in range(3)); return np.stack([B[a][sl] for a in range(3)])


def embedded_link_pos():
    return T53.link_pos(NB)


def probe_ring(ell, chi, axis):
    return T53.ring_state(NB, M, 5.0, 1.5, ell, chi, axis)


def spin_packet(sig, chi):
    c = NB // 2; psi, _, _ = S.packet_k(NB, M, 0.0, (c, c, c), sig, (0, 0, 0), chi=chi); return psi


def Eint(Bs, P):
    return -float(np.sum(Bs * P))


def torque(psi, Bs):
    dpsi = -1j * T53.Hb(psi, Bs); Lp = T53.angmom(T53.link_P(psi + dpsi)); Lm = T53.angmom(T53.link_P(psi - dpsi)); return (Lp - Lm) / 2


def sigma_density(psi):
    Sz = -1j * G.AL[0] @ G.AL[1]; Sx = -1j * G.AL[1] @ G.AL[2]; Sy = -1j * G.AL[2] @ G.AL[0]
    return np.stack([np.real(np.einsum("...i,ij,...j->...", psi.conj(), Sm, psi)) for Sm in (Sx, Sy, Sz)])


def cB_fit(psi):
    s = sigma_density(psi); P = T53.link_P(psi); g = lambda f, a: 0.5 * (np.roll(f, -1, a) - np.roll(f, 1, a))
    cs = np.stack([g(s[2], 1) - g(s[1], 2), g(s[0], 2) - g(s[2], 0), g(s[1], 0) - g(s[0], 1)])
    C = np.stack([0.5 * 0.5 * (cs[a] + np.roll(cs[a], -1, a)) for a in range(3)])  # 1/2 curl s, averaged to the link midpoint
    C[:, -1, :, :] = 0; C[:, :, -1, :] = 0; C[:, :, :, -1] = 0
    c = float(np.sum(P * C) / np.sum(C * C)); return c, float(np.linalg.norm(P - c * C) / np.linalg.norm(P))


def light_loop(N, B, x0, h=8, width=4.0):
    st = T53.square_loop(x0, h); Ns, bl = T53.loop_fields(st, N, B); ro = T53.readout(Ns, bl)
    dt = SY.wave(Ns, bl, +1, -1, width=width) - SY.wave(Ns, bl, -1, -1, width=width); A = (2 * h) ** 2
    Oloop = -float(np.sum(bl)) / (2 * A)
    return dict(dt=dt, readout=ro, wave_over_readout=dt / ro if ro != 0 else None, Omega_loop=Oloop, Omega_light=(dt / ro) * Oloop if ro != 0 else None,
                Omega_light_raw=-dt / (4 * A))


def grad_field(n, x0, amp, w=8.0):
    x = np.arange(n); X = np.stack(np.meshgrid(x, x, x, indexing="ij"), 0).astype(float); c0 = np.asarray(x0, float) + np.array([3.0, 2.0, 1.0])
    chi = amp * np.exp(-np.sum((X - c0[:, None, None, None]) ** 2, 0) / (2 * w * w)); del X
    return np.stack([np.roll(chi, -1, a) - chi for a in range(3)])


# ------------------------------------------------------------------ 57B
def part_B(n=129, u=0.01, ds=(24, 32)):
    c = n // 2; res = dict(fields={}, positions={}); t0 = time.time()
    q0, P0 = T53.ring_sources(n, 0); Q0 = float(q0.sum()); kap = 4 * np.pi * 12 * u / Q0; Gk = kap / (4 * np.pi)
    B0 = T53.shift_static(P0, kap, T53.angmom(P0)); res["C4_B0_max"] = float(np.abs(B0).max()); del B0, P0
    F = {}
    for ell in (2, -2):
        q, P = T53.ring_sources(n, ell); J = T53.angmom(P); N, info = G.static([], [], 0.0, kap, n, (c, c, c), qfix=q); B = T53.shift_static(P, kap, J)
        F[ell] = (q, N, B, J); res["fields"][str(ell)] = dict(J=J.tolist(), Q=float(q.sum()), Bmax=float(np.abs(B).max()))
        print("57B field ell", ell, "J", np.round(J, 6), "%.0fs" % (time.time() - t0), flush=True)
    res["C4_q_diff"] = float(np.abs(F[2][0] - F[-2][0]).max() / np.abs(F[2][0]).max()); res["C4_N_diff"] = float(np.abs(F[2][1] - F[-2][1]).max())
    Jz = F[2][3][2]; res["kappa"] = kap; res["G"] = Gk
    # probes (built once, flat N, box NB centred on the probe)
    probes = {}
    for ax_name, axis in (("z", (0, 0, 1)), ("x", (1, 0, 0))):
        for lp in (1, -1):
            probes[(ax_name, lp)] = [probe_ring(lp, chi, axis) for chi in (UP, DOWN)]
    Pz = {lp: [T53.link_P(p) for p in probes[("z", lp)]] for lp in (1, -1)}
    Lz = {lp: np.mean([T53.angmom(P) for P in Pz[lp]], 0) for lp in (1, -1)}
    Lx = {lp: np.mean([T53.angmom(T53.link_P(p)) for p in probes[("x", lp)]], 0) for lp in (1, -1)}
    res["probe_L"] = dict(z={str(k): v.tolist() for k, v in Lz.items()}, x={str(k): v.tolist() for k, v in Lx.items()})
    spins = {sig: {nm: spin_packet(sig, chi) for nm, chi in (("up", UP), ("down", DOWN))} for sig in (3.0, 5.0)}
    sz = {sig: {nm: float(np.sum(sigma_density(p)[2])) for nm, p in d.items()} for sig, d in spins.items()}
    res["spin_Sz"] = {str(k): v for k, v in sz.items()}; res["cB"] = {str(sig): cB_fit(spins[sig]["up"]) for sig in spins}
    print("57B cB", res["cB"], "Sz", res["spin_Sz"], flush=True)
    # C2 rigid rotation identity
    om = np.array([1e-4, -2e-4, 3e-4]); LP = embedded_link_pos(); Brot = np.stack([np.cross(om[:, None], LP[a].reshape(3, -1), axis=0).reshape(LP[a].shape)[a] for a in range(3)])
    e = Eint(Brot, Pz[1][0]); Lr = T53.angmom(Pz[1][0]); res["C2"] = dict(E=e, minus_omega_L=-float(om @ Lr), rel=abs(e + float(om @ Lr)) / abs(float(om @ Lr)))
    print("57B C2", res["C2"], flush=True)
    # C3 wave-clock instrument on L = 64
    one = np.ones(64); ba = 1e-3 * one; dtu = SY.wave(one, ba, +1, -1, width=4.0) - SY.wave(one, ba, -1, -1, width=4.0); res["C3"] = dict(uniform_over_2sum=dtu / (2 * ba.sum()))
    print("57B C3", res["C3"], flush=True)
    for d in ds:
        for key, x0 in (("axis", (c, c, c + d)), ("plane", (c + d, c, c))):
            row = {}; r = np.array(x0, float) - c; rh = r / np.linalg.norm(r); Jv = np.array([0, 0, Jz])
            row["Omega_LT"] = (Gk * (3 * (Jv @ rh) * rh - Jv) / np.linalg.norm(r) ** 3).tolist()
            for ell in (2, -2):
                q, N, B, J = F[ell]; Bs = sub(B, x0); o = {}
                o["Omega_curl"] = (-0.5 * curl_site(B, x0)).tolist()
                E = {lp: np.mean([Eint(Bs, P) for P in Pz[lp]]) for lp in (1, -1)}; o["Omega_orb"] = float((E[1] - E[-1]) / (Lz[1][2] - Lz[-1][2]))
                if key == "axis":
                    tq = {lp: np.mean([torque(p, Bs) for p in probes[("x", lp)]], 0) for lp in (1, -1)}
                    o["torque"] = {str(lp): tq[lp].tolist() for lp in (1, -1)}; o["Omega_torque"] = float((tq[1][1] - tq[-1][1]) / (Lx[1][0] - Lx[-1][0]))
                    for sig in spins:
                        Eu = Eint(Bs, T53.link_P(spins[sig]["up"])); Ed = Eint(Bs, T53.link_P(spins[sig]["down"]))
                        o[f"Omega_spin_{sig:g}"] = float(2 * (Eu - Ed) / (sz[sig]["up"] - sz[sig]["down"]))
                o["light"] = light_loop(N, B, x0)
                row[str(ell)] = o
            # pure-gradient control at d = 24 (ell = +2 amplitude)
            if d == ds[0]:
                Bs2 = sub(F[2][2], x0); amp = float(np.abs(Bs2[:, 12:36, 12:36, 12:36]).max()) * 8.0 * np.exp(0.5)
                Bg = grad_field(n, x0, amp); Bgs = sub(Bg, x0); gc = dict(amp=amp, Bmax=float(np.abs(Bgs).max()))
                Eg = {lp: np.mean([Eint(Bgs, P) for P in Pz[lp]]) for lp in (1, -1)}; gc["Omega_orb"] = float((Eg[1] - Eg[-1]) / (Lz[1][2] - Lz[-1][2]))
                if key == "axis":
                    for sig in spins:
                        Eu = Eint(Bgs, T53.link_P(spins[sig]["up"])); Ed = Eint(Bgs, T53.link_P(spins[sig]["down"]))
                        gc[f"Omega_spin_{sig:g}"] = float(2 * (Eu - Ed) / (sz[sig]["up"] - sz[sig]["down"]))
                gc["light"] = light_loop(np.ones((n, n, n)), Bg, x0); row["grad"] = gc; del Bg
            res["positions"][f"{key}_{d}"] = row
            print("57B", key, d, json.dumps({k: (v if k in ("Omega_LT",) else {kk: vv for kk, vv in v.items() if kk not in ("torque",)}) for k, v in row.items()}, default=float), flush=True)
    return res


def save(name, res, t0):
    res["_meta"] = dict(wall_s=time.time() - t0, platform=platform.platform(), python=platform.python_version(), numpy=np.__version__)
    (OUT / f"EX-GC57_{name}.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True); job = sys.argv[1]; t0 = time.time()
    save(job, {"A": part_A, "B": part_B}[job](), t0); print("done", job, "%.0fs" % (time.time() - t0))
