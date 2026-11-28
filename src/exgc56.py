"""EX-GC56: the link shift as a substrate drift. 56A drift universality in uniform beta (Dirac matter vs light; P vs J_E coupling);
56B conveyor identity (-i H_beta psi = (beta/2)(T - T^+) psi for P); 56C Codex generalisation of the momentum rule (kappa/2, V = 0.15,
sigma = 4) in the 55C pair. Plan: spec/EX-GC56_plan.md. EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G
import exgc52 as S
import exgc53 as T53
import exgc54 as X54
import exgc55 as X55

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC56"
shp = G.shp; BETA = 1e-3


# ------------------------------------------------------------------ 56A matter drift
def je_r(psi, N, m, r, B):
    """J_E/J_E coupling for general r: O_a = (1/2){H, v_a} with v_a = i[H, X_a] built numerically from H and the coordinate."""
    n = N.shape[0]; X = np.arange(n, dtype=float); b = X55.site_b(B); out = np.zeros_like(psi)
    for a in range(3):
        if not np.any(b[a]):
            continue
        Xa = X.reshape([-1 if i == a else 1 for i in range(3)])[..., None]
        v = lambda ph: 1j * (G.Hm(Xa * ph, N, m, r) - Xa * G.Hm(ph, N, m, r))
        O = lambda ph: 0.5 * (G.Hm(v(ph), N, m, r) + v(G.Hm(ph, N, m, r)))
        out -= 0.5 * (b[a][..., None] * O(psi) + O(b[a][..., None] * psi))
    return out


def centroid(psi):
    w = np.sum(np.abs(psi) ** 2, -1); x = np.arange(psi.shape[0]); return float(np.sum(w * x[:, None, None]) / w.sum())


def evolve_static(psi, Hf, T, dt):
    t = 0.0
    while t < T - 1e-9:
        k1 = -1j * Hf(psi); k2 = -1j * Hf(psi + dt / 2 * k1); k3 = -1j * Hf(psi + dt / 2 * k2); k4 = -1j * Hf(psi + dt * k3)
        psi = psi + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt
    return psi


def matter_case(m, k0, r, rule, chi=(1, 0, 0, 0), n=48, T=10.0, dt=0.1, beta=BETA):
    c = n // 2; one = np.ones((n, n, n)); psi, K, w = S.packet_k(n, m, r, (c, c, c), S.SIG, (k0, 0, 0), chi=chi)
    B = np.zeros((3, n, n, n)); B[0] = beta
    Hf = lambda p, BB: (G.Hm(p, one, m, r) + (T53.Hb(p, BB) if rule == "P" else je_r(p, one, m, r, BB))) if BB is not None else G.Hm(p, one, m, r)
    p0 = evolve_static(psi, lambda p: Hf(p, None), T, dt); pb = evolve_static(psi, lambda p: Hf(p, B), T, dt)
    drift = (centroid(pb) - centroid(p0)) / (beta * T)
    h = 1e-4; e1 = np.array([h, 0, 0]); E2 = lambda KK: G.disp(KK, 1.0, m, r) ** 2
    pred_P = -float(np.sum(w * np.cos(K[..., 0])))
    pred_JE = -float(np.sum(w * 0.5 * (E2(K + e1) + E2(K - e1) - 2 * E2(K)) / h ** 2))
    return dict(m=m, k0=k0, r=r, rule=rule, chi=list(chi), drift=drift, pred_P=pred_P, pred_JE=pred_JE,
                ratio_P=drift / pred_P, ratio_JE=drift / pred_JE)


# ------------------------------------------------------------------ 56A light drift
def light_drift(k0, L=512, width=16.0, T=100.0, dt=0.05, beta=BETA):
    s = np.arange(L); one = np.ones(L)
    def run(b):
        bl = b * one; bs = 0.5 * (bl + np.roll(bl, 1)); grad = lambda a: 0.5 * (np.roll(a, -1) - np.roll(a, 1))
        def rhs(u, w):
            lap = (np.roll(u, -1) - u) + (np.roll(u, 1) - u); return w + bs * grad(u), lap + grad(bs * w)
        d = (s - L // 2 + L) % L - L // 2; env = np.exp(-d ** 2 / (2 * width ** 2)); om = 2 * np.sin(k0 / 2)
        u = env * np.cos(k0 * d); ut = om * env * np.sin(k0 * d); w = ut - bs * grad(u); ang = 2 * np.pi * s / L
        def cen(u, w):
            gl = (np.roll(u, -1) - u) ** 2; e = 0.5 * w ** 2 + 0.25 * (gl + np.roll(gl, 1)) + bs * w * grad(u)
            return np.angle(np.sum(e * np.exp(1j * ang))) * L / (2 * np.pi)
        t = 0.0; prev = cen(u, w); unw = 0.0
        while t < T - 1e-9:
            k1 = rhs(u, w); k2 = rhs(u + dt / 2 * k1[0], w + dt / 2 * k1[1]); k3 = rhs(u + dt / 2 * k2[0], w + dt / 2 * k2[1]); k4 = rhs(u + dt * k3[0], w + dt * k3[1])
            u = u + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); w = w + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]); t += dt
            c = cen(u, w); unw += (c - prev + L / 2) % L - L / 2; prev = c
        return unw, u
    x0, u0 = run(0.0); xb, _ = run(beta); drift = (xb - x0) / (beta * T)
    d = (s - L // 2 + L) % L - L // 2; env = np.exp(-d ** 2 / (2 * width ** 2)); uk = np.fft.fft(env * np.exp(1j * k0 * d)); kk = 2 * np.pi * np.fft.fftfreq(L)
    wgt = np.abs(uk) ** 2; pred = -float(np.sum(wgt * np.cos(kk)) / wgt.sum())
    return dict(k0=k0, drift=drift, pred=pred, ratio=drift / pred, pred_point=-float(np.cos(k0)))


def part_A():
    res = dict(matter=[], light=[])
    for rule in ("P", "JE"):
        for r in (0.0, 1.0):
            for m in (0.3, 0.6, 1.2):
                for k0 in (0.0, 0.15, 0.3):
                    row = matter_case(m, k0, r, rule); res["matter"].append(row)
                    print("56A %s r %.0f m %.1f k0 %.2f drift %.5f predP %.5f predJE %.5f  rP %.5f rJE %.5f" % (rule, r, m, k0, row["drift"], row["pred_P"], row["pred_JE"], row["ratio_P"], row["ratio_JE"]), flush=True)
        row = matter_case(0.6, 0.15, 0.0, rule, chi=(0, 1, 0, 0)); res["matter"].append(row); print("56A spin down", rule, row["ratio_P"], row["ratio_JE"], flush=True)
    for k0 in (0.25, 0.5):
        row = light_drift(k0); res["light"].append(row); print("56A light", row, flush=True)
    # C1 zero-beta, C2 dt/2
    c = 24; one = np.ones((48, 48, 48)); psi, _, _ = S.packet_k(48, 0.6, 0.0, (c, c, c), S.SIG, (0.15, 0, 0))
    pa = evolve_static(psi, lambda p: G.Hm(p, one, 0.6, 0.0), 10.0, 0.1); pb = evolve_static(psi, lambda p: G.Hm(p, one, 0.6, 0.0), 10.0, 0.1)
    res["C1"] = abs(centroid(pa) - centroid(pb))
    d1 = matter_case(0.6, 0.15, 0.0, "P")["drift"]; d2 = matter_case(0.6, 0.15, 0.0, "P", dt=0.05)["drift"]; res["C2"] = dict(dt=d1, half=d2, rel=abs(d2 / d1 - 1))
    print("56A C1", res["C1"], "C2", res["C2"], flush=True)
    return res


# ------------------------------------------------------------------ 56B conveyor identity
def part_B():
    rng = np.random.default_rng(56); n = 12; m = 0.6; one = np.ones((n, n, n))
    psi = rng.standard_normal((n, n, n, 4)) + 1j * rng.standard_normal((n, n, n, 4)); B = rng.standard_normal((3, n, n, n))
    lhs = -1j * T53.Hb(psi, B)
    conv = sum(0.5 * (B[a][..., None] * shp(psi, a, 1) - shp(B[a], a, -1)[..., None] * shp(psi, a, -1)) for a in range(3))
    res = dict(P_identity=float(np.linalg.norm(lhs - conv) / np.linalg.norm(lhs)), P_real_transfer=True)
    Bu = np.zeros((3, n, n, n)); Bu[0] = 1.0; je = -1j * X55.je_terms(psi, one, m, Bu)[0]
    adv = 0.5 * (shp(psi, 0, 1) - shp(psi, 0, -1)); cfit = np.vdot(adv, je) / np.vdot(adv, adv)
    res["JE_fit_coef"] = [float(np.real(cfit)), float(np.imag(cfit))]; res["JE_residual"] = float(np.linalg.norm(je - cfit * adv) / np.linalg.norm(je))
    Bs = 1e-3 * rng.standard_normal((3, n, n, n)); res["je_r_vs_55"] = float(np.linalg.norm(je_r(psi, one, m, 0.0, Bs) - X55.je_terms(psi, one, m, Bs)[0]) / np.linalg.norm(X55.je_terms(psi, one, m, Bs)[0]))
    print("56B", res, flush=True)
    return res


# ------------------------------------------------------------------ 56C generalisation (P/P)
def run_pair(V, kscale=1.0, sigma=None, n=49, sep=14, T=8.0, dt=0.1, lam=4.0):
    sig = S.SIG if sigma is None else sigma; c = n // 2; m = S.M; k0 = S.k_for_v(V, m, 0.0)
    p1, _, _ = S.packet_k(n, m, 0.0, (c, c - sep // 2, c), sig, (k0, 0, 0)); p2, _, _ = S.packet_k(n, m, 0.0, (c, c + sep // 2, c), sig, (k0, 0, 0))
    p0, _, _ = S.packet_k(n, m, 0.0, (c, c - sep // 2, c), S.SIG, (0, 0, 0)); kap = kscale * 4 * np.pi * sep * 0.01 / float(G.charge(p0, np.ones((n, n, n)), m, 0.0).sum())
    out = dict(V=V, kscale=kscale, sigma=sig, n=n, T=T, dt=dt, kappa=kap, runs={})
    for lab, ps in (("pair", [p1, p2]), ("one", [p1])):
        t0 = time.time(); f, g, info = S.initial_field(ps, [m] * len(ps), 0.0, kap, n, c, V)
        j0 = sum(T53.link_P(p) for p in ps); B, gB = X54.initial_shift(j0, kap, lam, n, V)
        _, hist = X55.evolve(ps, [m] * len(ps), f, g, B, gB, kap, lam, T, dt, "P")
        out["runs"][lab] = dict(t=[h["t"] for h in hist], Py=[h["Py"] for h in hist], wall=time.time() - t0)
        print("56C V %.2f k %.2f sig %.1f %s %.0fs" % (V, kscale, sig, lab, time.time() - t0), flush=True)
    t = np.array(out["runs"]["pair"]["t"]); P1 = np.array([p[0] for p in out["runs"]["pair"]["Py"]]) - np.array([p[0] for p in out["runs"]["one"]["Py"]])
    out["F1"] = float(np.polyfit(t, P1, 2)[1]); print("56C V %.2f k %.2f sig %.1f F1 %.6e" % (V, kscale, sig, out["F1"]), flush=True)
    return out


def save(name, res, t0):
    res["_meta"] = dict(wall_s=time.time() - t0, platform=platform.platform(), python=platform.python_version(), numpy=np.__version__)
    (OUT / f"EX-GC56_{name}.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True); job = sys.argv[1]; t0 = time.time()
    jobs = {"A": part_A, "B": part_B, "Ck_0.0": lambda: run_pair(0.0, kscale=0.5), "Ck_0.2": lambda: run_pair(0.2, kscale=0.5),
            "Cv_0.15": lambda: run_pair(0.15), "Cs_0.0": lambda: run_pair(0.0, sigma=4.0), "Cs_0.2": lambda: run_pair(0.2, sigma=4.0),
            "C3_0.2": lambda: run_pair(0.2)}
    save(job, jobs[job](), t0); print(job, "done %.0fs" % (time.time() - t0), flush=True)
