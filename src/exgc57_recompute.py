"""EX-GC57 recomputation through separate paths -> results/EX-GC57/recompute.json.
57A: the one-way launch (a) time-stepped by RK4 with 56's right-hand side (not the exact Fourier propagator) at width 16, both k0;
     predictor recomputed from the analytic mode energy on a fine k grid (continuous Fourier transform of the Gaussian).
57B: matter energies from <psi|H_beta psi> (53's Hb operator) instead of -sum beta.P; the curl from 53's loop machinery on four
     2x2 plaquette loops (loop_fields), not curl_site. Shared dependencies: 53's ring sources, static N and shift, probe builders.
EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G
import exgc53 as T53
import exgc57 as E57
OUT = E57.OUT


def rk4_drift(k0, width, beta, L=512, T=100.0, dt=0.05):
    s = np.arange(L); d = (s - L // 2 + L) % L - L // 2; env = np.exp(-d ** 2 / (2 * width ** 2)); k = 2 * np.pi * np.fft.fftfreq(L)
    A = np.fft.fft(env * np.exp(1j * k0 * d)) * (k > 0); u0 = 2 * np.real(np.fft.ifft(A)); w0 = 2 * np.real(np.fft.ifft(-2j * np.sin(k / 2) * A))
    grad = lambda a: 0.5 * (np.roll(a, -1) - np.roll(a, 1)); ang = 2 * np.pi * s / L
    def run(b):
        rhs = lambda u, w: (w + b * grad(u), (np.roll(u, -1) - 2 * u + np.roll(u, 1)) + b * grad(w))
        def cen(u, w):
            gl = (np.roll(u, -1) - u) ** 2; e = 0.5 * w ** 2 + 0.25 * (gl + np.roll(gl, 1)) + b * w * grad(u); return np.angle(np.sum(e * np.exp(1j * ang))) * L / (2 * np.pi)
        u, w = u0.copy(), w0.copy(); prev = cen(u, w); un = 0.0
        for _ in range(int(round(T / dt))):
            k1 = rhs(u, w); k2 = rhs(u + dt / 2 * k1[0], w + dt / 2 * k1[1]); k3 = rhs(u + dt / 2 * k2[0], w + dt / 2 * k2[1]); k4 = rhs(u + dt * k3[0], w + dt * k3[1])
            u = u + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); w = w + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]); c = cen(u, w); un += (c - prev + L / 2) % L - L / 2; prev = c
        return un
    return (run(beta) - run(0.0)) / (beta * T)


def pred_cont(k0, width, beta):
    k = np.linspace(1e-6, np.pi, 200001); A2 = np.exp(-(k - k0) ** 2 * width ** 2) + 0 * k; O = 2 * np.sin(k / 2)
    eps = lambda b: O ** 2 * A2 - b * O * np.sin(k) * A2; vg = lambda b: np.cos(k / 2) - b * np.cos(k)
    return float((np.sum(eps(beta) * vg(beta)) / eps(beta).sum() - np.sum(eps(0) * vg(0)) / eps(0).sum()) / beta)


def curl_loops(B, x0):
    one = np.ones(B[0].shape); vals = []
    for da in (-1, 0):
        for db in (-1, 0):
            p = (x0[0] + da, x0[1] + db, x0[2]); px = (p[0] + 1, p[1], p[2]); pxy = (p[0] + 1, p[1] + 1, p[2]); py = (p[0], p[1] + 1, p[2])
            st = [(p, 0, 1), (px, 1, 1), (pxy, 0, -1), (py, 1, -1)]
            _, bl = T53.loop_fields(st, one, B); vals.append(np.sum(bl))
    return -0.5 * float(np.mean(vals))


if __name__ == "__main__":
    res = dict(A={})
    for k0 in (0.25, 0.5):
        dr = rk4_drift(k0, 16.0, 1e-3); pc = pred_cont(k0, 16.0, 1e-3); res["A"][str(k0)] = dict(rk4=dr, pred_cont=pc, ratio=dr / pc); print("A", k0, res["A"][str(k0)], flush=True)
    n = 129; c = n // 2; q0, _ = T53.ring_sources(n, 0); kap = 4 * np.pi * 12 * 0.01 / float(q0.sum())
    q, P = T53.ring_sources(n, 2); B = T53.shift_static(P, kap, T53.angmom(P))
    pz = {lp: [E57.probe_ring(lp, chi, (0, 0, 1)) for chi in (E57.UP, E57.DOWN)] for lp in (1, -1)}
    Lz = {lp: np.mean([T53.angmom(T53.link_P(p))[2] for p in pz[lp]]) for lp in (1, -1)}
    sp = {s: [E57.spin_packet(s, chi) for chi in (E57.UP, E57.DOWN)] for s in (3.0, 5.0)}
    Sz = {s: [float(np.sum(np.real(np.conj(p) * (p @ (-1j * G.AL[0] @ G.AL[1]).T)))) for p in sp[s]] for s in sp}
    res["B"] = {}
    for key, x0 in (("axis_24", (c, c, c + 24)), ("plane_24", (c + 24, c, c)), ("axis_32", (c, c, c + 32)), ("plane_32", (c + 32, c, c))):
        Bs = E57.sub(B, x0); H = lambda p: float(np.real(np.vdot(p, T53.Hb(p, Bs))))
        E = {lp: np.mean([H(p) for p in pz[lp]]) for lp in (1, -1)}; oc = curl_loops(B, x0); r = dict(Omega_curl=oc, orb_over_curl=(E[1] - E[-1]) / (Lz[1] - Lz[-1]) / oc)
        if key.startswith("axis"):
            for s in sp:
                r[f"spin_{s:g}_over_curl"] = 2 * (H(sp[s][0]) - H(sp[s][1])) / (Sz[s][0] - Sz[s][1]) / oc
        res["B"][key] = r; print("B", key, r, flush=True)
    (OUT / "recompute.json").write_text(json.dumps(res, indent=1, default=float))
