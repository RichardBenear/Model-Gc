"""EX-GC55 post-hoc diagnostic (not scored): emit both pulses at the same coordinate frequency as the clock.
55A's registered readout launches both directions with the same wavenumber k0. Where the shift b0 at the clock is nonzero,
the local lattice dispersion omega = 2 sin(k/2) N^2 -/+ b0 sin k then gives the two pulses different frequencies, and lattice
dispersion turns that into a delay of order L * b0. A clock fixes its frequency, not its wavenumber: here each direction's k
solves omega(k; b0) = omega0 = 2 sin(k0/2) at the clock (b0 = site shift at s = 0)."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc55 as X
import exgc53 as T53
import exgc51 as G

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC55"


def k_match(om0, b0, N0, sgn):
    lo, hi = 1e-6, 1.5
    f = lambda k: 2 * np.sin(k / 2) * N0 ** 2 - sgn * b0 * np.sin(k) - om0
    for _ in range(100):
        mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return 0.5 * (lo + hi)


def wave(Ns, bl, direction, rule="cov", width=16.0, dt=0.05):
    L = len(Ns); s = np.arange(L); bs = 0.5 * (bl + np.roll(bl, 1)); Nl = Ns * np.roll(Ns, -1)
    grad = lambda a: 0.5 * (np.roll(a, -1) - np.roll(a, 1))
    def rhs(u, w):
        lap = Nl * (np.roll(u, -1) - u) + np.roll(Nl, 1) * (np.roll(u, 1) - u)
        return Ns * w + bs * grad(u), Ns * lap + (grad(bs * w) if rule == "cov" else bs * grad(w))
    om0 = 2 * np.sin(X.K0 / 2) * Ns[0] ** 2; k = k_match(om0, bs[0], Ns[0], direction)
    d = (s - L // 2 + L) % L - L // 2; env = np.exp(-d ** 2 / (2 * width ** 2)); kk = direction * k
    u = env * np.cos(kk * d); ut = om0 * env * np.sin(kk * d); w = (ut - bs * grad(u)) / Ns
    ang = 2 * np.pi * s / L; t = 0.0; unw = 0.0; target = L
    def cen(u, w):
        gl = Nl * (np.roll(u, -1) - u) ** 2; e = 0.5 * Ns * w ** 2 + 0.25 * (gl + np.roll(gl, 1)) + bs * w * grad(u)
        return np.angle(np.sum(e * np.exp(1j * ang))) * L / (2 * np.pi)
    prev = cen(u, w)
    while t < 10 * L:
        k1 = rhs(u, w); k2 = rhs(u + dt / 2 * k1[0], w + dt / 2 * k1[1]); k3 = rhs(u + dt / 2 * k2[0], w + dt / 2 * k2[1]); k4 = rhs(u + dt * k3[0], w + dt * k3[1])
        u = u + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); w = w + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]); t += dt
        c = cen(u, w); dc = (c - prev + L / 2) % L - L / 2; un = unw + dc
        if direction * un >= target:
            return t - dt * (direction * un - target) / (direction * dc)
        unw, prev = un, c


def ldt(Ns, bl, **kw):
    return wave(Ns, bl, +1, **kw) - wave(Ns, bl, -1, **kw)


if __name__ == "__main__":
    L = 256; one = np.ones(L); s = np.arange(L); res = {}
    ba = 1e-3 * one; da = ldt(one, ba); res["uniform_ratio_over_LF"] = da / T53.readout(one, ba) / X.LF
    bc = np.where(s < 64, 1e-3, 0.0); res["one_side_over_quarter"] = ldt(one, bc) / da / 0.25
    for rule in ("cov", "adv"):
        for Lam, A, ph in [(L_, 1e-3, 0.0) for L_ in (16, 32, 64, 128)] + [(64, 1e-3, np.pi / 4), (64, 1e-3, 3 * np.pi / 4), (64, 1e-2, 0.0)]:
            chi = A * np.sin(2 * np.pi * s / Lam + ph); bl = np.roll(chi, -1) - chi; norm = abs(da) * np.abs(bl).sum() / np.abs(ba).sum()
            r = ldt(one, bl, rule=rule) / norm; res[f"{rule}_{Lam}_{A:g}_{ph:.3f}"] = r; print(rule, Lam, A, round(ph, 3), "rel %.2e" % r, flush=True)
    print("uniform/LF %.5f one-side/0.25 %.4f" % (res["uniform_ratio_over_LF"], res["one_side_over_quarter"]), flush=True)
    if "ring" in sys.argv:
        n = 129; c = n // 2; q0, _ = T53.ring_sources(n, 0); kap = 4 * np.pi * 12 * 0.01 / float(q0.sum())
        for ell in (2, -2):
            q, P = T53.ring_sources(n, ell); J = T53.angmom(P); N, _ = G.static([], [], 0.0, kap, n, (c, c, c), qfix=q); B = T53.shift_static(P, kap, J)
            for h in (20, 28, 36):
                Ns, bl = T53.loop_fields(T53.square_loop((c, c, c), h), N, B); ro = T53.readout(Ns, bl); w = ldt(Ns, bl)
                res[f"ring_{ell}_{h}"] = w / ro / X.LF; print("ring", ell, h, "wave/readout/LF %.5f" % (w / ro / X.LF), flush=True)
    (OUT / "posthoc_freq_matched.json").write_text(json.dumps(res, indent=1, default=float))
