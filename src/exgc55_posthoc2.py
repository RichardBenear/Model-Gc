"""EX-GC55 post-hoc diagnostic (not scored): prepare the clock's pulse in the clock's locally synchronised frame.
u(s,0) = u0(s) + sg*chi_loc(s) u0_t(s), u_t(s,0) = u0_t + sg*chi_loc u0_tt, with chi_loc(d) = sum of the link shift from the clock
(s = 0) to d along the loop (Einstein synchronisation near the clock; for a pure-gauge pattern chi_loc = chi - chi(0)).
The sign sg is fixed on the gauge cases; uniform, one-sided and the rotating ring are then positive controls."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc55 as X
import exgc53 as T53
import exgc51 as G
OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC55"


def wave(Ns, bl, direction, sg, width=16.0, dt=0.05, k0=X.K0):
    L = len(Ns); s = np.arange(L); bs = 0.5 * (bl + np.roll(bl, 1)); Nl = Ns * np.roll(Ns, -1)
    grad = lambda a: 0.5 * (np.roll(a, -1) - np.roll(a, 1))
    lapf = lambda u: Nl * (np.roll(u, -1) - u) + np.roll(Nl, 1) * (np.roll(u, 1) - u)
    def rhs(u, w):
        return Ns * w + bs * grad(u), Ns * lapf(u) + grad(bs * w)
    d = (s - L // 2 + L) % L - L // 2; env = np.exp(-d ** 2 / (2 * width ** 2)); kk = direction * k0; om = 2 * np.sin(k0 / 2) * Ns[0] ** 2
    cs = np.cumsum(bl); chi = np.where(d >= 0, np.concatenate([[0.0], cs[:-1]])[s % L], 0.0)
    neg = d < 0; chi[neg] = -np.array([bl[(L + dd) % L:].sum() if dd < 0 else 0.0 for dd in d[neg]])
    u0 = env * np.cos(kk * d); u0t = om * env * np.sin(kk * d)
    u0tt = Ns * lapf(u0) * Ns                                                    # beta = 0 wave equation, N uniform at the clock
    u = u0 + sg * chi * u0t; ut = u0t + sg * chi * u0tt; w = (ut - bs * grad(u)) / Ns
    ang = 2 * np.pi * s / L; t = 0.0; unw = 0.0
    def cen(u, w):
        gl = Nl * (np.roll(u, -1) - u) ** 2; e = 0.5 * Ns * w ** 2 + 0.25 * (gl + np.roll(gl, 1)) + bs * w * grad(u)
        return np.angle(np.sum(e * np.exp(1j * ang))) * L / (2 * np.pi)
    prev = cen(u, w)
    while t < 10 * L:
        k1 = rhs(u, w); k2 = rhs(u + dt / 2 * k1[0], w + dt / 2 * k1[1]); k3 = rhs(u + dt / 2 * k2[0], w + dt / 2 * k2[1]); k4 = rhs(u + dt * k3[0], w + dt * k3[1])
        u = u + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); w = w + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]); t += dt
        c = cen(u, w); dc = (c - prev + L / 2) % L - L / 2; un = unw + dc
        if direction * un >= L:
            return t - dt * (direction * un - L) / (direction * dc)
        unw, prev = un, c


def ldt(Ns, bl, sg):
    return wave(Ns, bl, +1, sg) - wave(Ns, bl, -1, sg)


if __name__ == "__main__":
    L = 256; one = np.ones(L); s = np.arange(L); res = {}
    ba = 1e-3 * one
    for sg in (+1, -1):
        da = ldt(one, ba, sg); r = dict(uniform_over_LF=da / T53.readout(one, ba) / X.LF, one_side_over_quarter=ldt(one, np.where(s < 64, 1e-3, 0.0), sg) / da / 0.25, gauge={})
        for Lam, A, ph in [(L_, 1e-3, 0.0) for L_ in (16, 32, 64, 128)] + [(64, 1e-3, np.pi / 4), (64, 1e-3, 3 * np.pi / 4), (64, 1e-2, 0.0)]:
            chi = A * np.sin(2 * np.pi * s / Lam + ph); bl = np.roll(chi, -1) - chi; norm = abs(da) * np.abs(bl).sum() / np.abs(ba).sum()
            r["gauge"][f"{Lam}_{A:g}_{ph:.3f}"] = ldt(one, bl, sg) / norm
        res[str(sg)] = r; print(sg, json.dumps(r), flush=True)
    if "ring" in sys.argv:
        sg = int(sys.argv[sys.argv.index("ring") + 1]); n = 129; c = n // 2; q0, _ = T53.ring_sources(n, 0); kap = 4 * np.pi * 12 * 0.01 / float(q0.sum()); rr = {}
        for ell in (2, -2):
            q, P = T53.ring_sources(n, ell); J = T53.angmom(P); N, _ = G.static([], [], 0.0, kap, n, (c, c, c), qfix=q); B = T53.shift_static(P, kap, J)
            for h in (20, 28, 36):
                Ns, bl = T53.loop_fields(T53.square_loop((c, c, c), h), N, B); ro = T53.readout(Ns, bl); w = ldt(Ns, bl, sg)
                rr[f"{ell}_{h}"] = w / ro / X.LF; print("ring", ell, h, "%.5f" % rr[f"{ell}_{h}"], flush=True)
        res["ring_sg"] = sg; res["ring"] = rr
    (OUT / "posthoc_sync.json").write_text(json.dumps(res, indent=1, default=float))
