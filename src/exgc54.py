"""EX-GC54: where the link shift comes from. Two event-counting source rules with no free strength (Doppler-squared two-site
counting, lambda = 2 x 2 = 4): D2-q, source 4 kappa J_q with J_q = (1/2){Q, v}; D2-E, source 4 kappa J_E with J_E = (1/2){H, v}.
54A no-go (time-ordered N-history offset is a gradient), 54B both rules in the 53E co-moving pair, 54C lambda = 2 with 53's law,
54D smooth gauge patterns on the one-clock loop, operator checks (P2), C1/C3.
Plan: spec/EX-GC54_plan.md. EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G
import exgc52 as S
import exgc53 as T53
import exgc49 as F49

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC54"
fN, Nf, shp, mv = G.fN, G.Nf, G.shp, G.mv
UP = (1, 0, 0, 0)


# ------------------------------------------------------------------ operators (r = 0)
def Hgen(psi, msite, bond):
    """sum_x msite_x beta psi_x + sum_bonds W_{x,x+a} hop, with W = bond(N_x, N_{x+a}) given as arrays per axis (forward)."""
    out = msite[..., None] * mv(G.BE, psi)
    for a in range(3):
        Wf = bond[a]; Wb = np.roll(Wf, 1, axis=a)
        sl = [slice(None)] * 3; sl[a] = 0; Wb[tuple(sl)] = 0
        out += Wf[..., None] * mv(G.AL2I[a], shp(psi, a, 1)) - Wb[..., None] * mv(G.AL2I[a], shp(psi, a, -1))
    return out


def bonds_H(N):
    return [N * shp(N, a, 1) for a in range(3)]


def bonds_Q(N):
    return [(N + shp(N, a, 1)) * (shp(np.ones_like(N), a, 1)) for a in range(3)]


def Hop(psi, N, m):
    return Hgen(psi, m * N, bonds_H(N))


def Qop(psi, N, m):
    """Q = sum_x dH/dN_x: mass weight m, bond weight N_x + N_{x+a}."""
    return Hgen(psi, m * np.ones_like(N), bonds_Q(N))


def vop(psi, N, a):
    """v_a psi = i[H, X_a] psi = i (W_f AL2I psi(x+a) + W_b AL2I psi(x-a))."""
    Wf = N * shp(N, a, 1); Wb = np.roll(Wf, 1, axis=a); sl = [slice(None)] * 3; sl[a] = 0; Wb[tuple(sl)] = 0
    return 1j * (Wf[..., None] * mv(G.AL2I[a], shp(psi, a, 1)) + Wb[..., None] * mv(G.AL2I[a], shp(psi, a, -1)))


def rule_current(psi, N, m, rule):
    """link-local current s_{x,a} = (1/2) Re[psi_x^+ (O psi)_x + psi_{x+a}^+ (O psi)_{x+a}], O = (1/2){A, v_a}, A = Q (D2q) or H (D2E)."""
    A = Qop if rule == "D2q" else Hop
    Ap = A(psi, N, m); s = np.zeros((3,) + N.shape)
    for a in range(3):
        Op = 0.5 * (A(vop(psi, N, a), N, m) + vop(Ap, N, a))
        d = np.real(np.sum(psi.conj() * Op, -1)); s[a] = 0.5 * (d + shp(d, a, 1))
    return s


def source(ps, N, ms, rule):
    if rule == "P":
        return sum(T53.link_P(p) for p in ps)
    return sum(rule_current(p, N, mm, rule) for p, mm in zip(ps, ms))


# ------------------------------------------------------------------ dynamics (53E with a general source)
def initial_shift(j, kap, lam, n, V):
    Ptot = j.sum(axis=(1, 2, 3)); gam = 1 / np.sqrt(1 - V * V); LP = T53.link_pos(n); B = np.zeros_like(j)
    for a in range(3):
        X = LP[a]; Rg = np.sqrt((gam * X[0]) ** 2 + X[1] ** 2 + X[2] ** 2); Rg = np.where(Rg == 0, 1e-9, Rg)
        B[a] = T53.poisson(lam * kap * j[a], -lam * kap * Ptot[a] * gam / (4 * np.pi * Rg), V)
    gB = np.zeros_like(B)
    for a in range(3):
        gB[a][1:-1] = -V * (B[a][2:] - B[a][:-2]) / 2
    return B, gB


def evolve(psis, ms, f, g, B, gB, kappa, lam, T, dt, rule):
    """RK4 as 53E. For rule P the source is evaluated at every stage (identical to 53E). For D2q/D2E the rule source is
    refreshed once per step from the step-start state and held across the four stages (amendment A1)."""
    I = (slice(1, -1),) * 3
    def faces0(a):
        for ax in range(3):
            sl = [slice(None)] * 3; sl[ax] = 0; a[tuple(sl)] = 0; sl[ax] = -1; a[tuple(sl)] = 0
        return a
    def rhs(ps, f, g, B, gB, jfix):
        N = Nf(f); dps = [-1j * (G.Hm(p, N, mm, 0.0) + T53.Hb(p, B)) for p, mm in zip(ps, ms)]
        q = sum(G.charge(p, N, mm, 0.0) for p, mm in zip(ps, ms)); df = faces0(N * g); dg = np.zeros_like(g); dg[I] = N[I] * (F49.op(N, f) - kappa * q[I])
        j = sum(T53.link_P(p) for p in ps) if jfix is None else jfix
        dB = np.zeros_like(B); dgB = np.zeros_like(gB)
        for a in range(3):
            Nb = 0.5 * (N + np.roll(N, -1, a)); dB[a] = faces0(Nb * gB[a]); dgB[a][I] = Nb[I] * (F49.op(Nb, B[a]) - lam * kappa * j[a][I])
        return dps, df, dg, dB, dgB
    def add(st, k, h):
        return [p + h * kp for p, kp in zip(st[0], k[0])], st[1] + h * k[1], st[2] + h * k[2], st[3] + h * k[3], st[4] + h * k[4]
    st = (psis, f, g, B, gB); t = 0.0; hist = [(0.0, [S.momentum(p, 1) for p in psis])]
    while t < T - 1e-9:
        jfix = None if rule == "P" else source(st[0], Nf(st[1]), ms, rule)
        k1 = rhs(*st, jfix); k2 = rhs(*add(st, k1, dt / 2), jfix); k3 = rhs(*add(st, k2, dt / 2), jfix); k4 = rhs(*add(st, k3, dt), jfix)
        st = ([p + dt / 6 * (a + 2 * b + 2 * c + d) for p, a, b, c, d in zip(st[0], k1[0], k2[0], k3[0], k4[0])],) + tuple(
            st[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(1, 5))
        t += dt; hist.append((t, [S.momentum(p, 1) for p in st[0]]))
    return st, hist


def packets(V, n=65, sep=14, chi=UP):
    c = n // 2; m = S.M; k0 = S.k_for_v(V, m, 0.0)
    p1, _, _ = S.packet_k(n, m, 0.0, (c, c - sep // 2, c), S.SIG, (k0, 0, 0), chi=chi); p2, _, _ = S.packet_k(n, m, 0.0, (c, c + sep // 2, c), S.SIG, (k0, 0, 0), chi=chi)
    p0, _, _ = S.packet_k(n, m, 0.0, (c, c - sep // 2, c), S.SIG, (0, 0, 0)); kap = 4 * np.pi * sep * 0.01 / float(G.charge(p0, np.ones((n, n, n)), m, 0.0).sum())
    return p1, p2, kap


def run_pair(V, rule, lam=4.0, n=65, T=16.0, dt=0.1, only=None):
    c = n // 2; m = S.M; p1, p2, kap = packets(V, n)
    out = dict(V=V, rule=rule, lam=lam, kappa=kap, n=n, T=T, dt=dt, runs={})
    for lab, ps in (("pair", [p1, p2]), ("one", [p1])):
        if only and lab not in only:
            continue
        t0 = time.time(); f, g, info = S.initial_field(ps, [m] * len(ps), 0.0, kap, n, c, V)
        j0 = source(ps, Nf(f), [m] * len(ps), rule); B, gB = initial_shift(j0, kap, lam, n, V)
        _, hist = evolve(ps, [m] * len(ps), f, g, B, gB, kap, lam, T, dt, rule)
        out["runs"][lab] = dict(t=[h[0] for h in hist], Py=[h[1] for h in hist], wall=time.time() - t0, j0_tot=j0.sum(axis=(1, 2, 3)).tolist())
        print("54 V %.1f rule %s lam %.1f %s %.0fs" % (V, rule, lam, lab, time.time() - t0), flush=True)
    if "pair" in out["runs"] and "one" in out["runs"]:
        t = np.array(out["runs"]["pair"]["t"]); P1 = np.array([p[0] for p in out["runs"]["pair"]["Py"]]) - np.array([p[0] for p in out["runs"]["one"]["Py"]])
        out["F1"] = float(np.polyfit(t, P1, 2)[1]); print("54 V %.1f rule %s lam %.1f F1 %.6e" % (V, rule, lam, out["F1"]), flush=True)
    return out


# ------------------------------------------------------------------ P2 operator checks, C3 ring
def part_ops():
    res = {}
    n = 65; c = n // 2; m = S.M
    for V in (0.0, 0.2):
        k0 = S.k_for_v(V, m, 0.0); psi, K, w = S.packet_k(n, m, 0.0, (c, c, c), S.SIG, (k0, 0, 0)); one = np.ones((n, n, n))
        h = 1e-5; e1 = np.array([h, 0, 0])
        q = (G.disp(K, 1 + h, m, 0) - G.disp(K, 1 - h, m, 0)) / (2 * h); v = (G.disp(K + e1, 1, m, 0) - G.disp(K - e1, 1, m, 0)) / (2 * h); E = G.disp(K, 1, m, 0)
        P = float(T53.link_P(psi).sum(axis=(1, 2, 3))[0]); row = dict(P=P, P_k=float(np.sum(w * np.sin(K[..., 0]))))
        for rule, A, kv in (("D2q", Qop, float(np.sum(w * q * v))), ("D2E", Hop, float(np.sum(w * E * v)))):
            s = rule_current(psi, one, m, rule); Ap = A(psi, one, m); Ov = 0.5 * (A(vop(psi, one, 0), one, m) + vop(Ap, one, 0))
            ex = float(np.real(np.vdot(psi, Ov))); tot = float(s[0].sum())
            row[rule] = dict(sum_s=tot, expect=ex, rel_sum=abs(tot - ex) / max(abs(ex), 1e-30), kspace=kv, rel_k=(tot / kv - 1) if V > 0 else None,
                             lam_eff=4 * tot / P if V > 0 else None, s_yz_max=float(np.abs(s[1:]).max()), s_x_max=float(np.abs(s[0]).max()))
        res[str(V)] = row; print("54 ops V", V, json.dumps(row), flush=True)
    # Hermiticity / Q identity: <Q> = sum charge at uniform N
    psi, _, _ = S.packet_k(n, m, 0.0, (c, c, c), S.SIG, (S.k_for_v(0.2, m, 0.0), 0, 0)); one = np.ones((n, n, n))
    res["Q_identity"] = abs(float(np.real(np.vdot(psi, Qop(psi, one, m)))) / float(G.charge(psi, one, m, 0.0).sum()) - 1)
    res["H_identity"] = abs(float(np.real(np.vdot(psi, Hop(psi, one, m)))) / G.energy(psi, one, m, 0.0) - 1)
    print("54 identities", res["Q_identity"], res["H_identity"], flush=True)
    return res


def part_ring(n=129):
    """C3 and the reported ring J: rule sources for the 53C spin-balanced ring, ell = +2, 0, -2 (uniform N; reported)."""
    res = {}; m = 0.6; one = np.ones((n, n, n))
    for ell in (2, 0, -2):
        t0 = time.time(); row = dict(J_P=T53.angmom(T53.ring_sources(n, ell)[1]).tolist()); acc = {"D2q": 0, "D2E": 0}
        for chi in (T53.UP, T53.DOWN):
            psi = T53.ring_state(n, m, 10.0, 3.0, ell, chi)
            for rule in acc:
                acc[rule] = acc[rule] + 0.5 * rule_current(psi, one, m, rule)
            del psi
        for rule, s in acc.items():
            row[rule] = dict(J=T53.angmom(s).tolist(), maxabs=float(np.abs(s).max()))
        res[str(ell)] = row; print("54 ring ell", ell, json.dumps(row), "%.0fs" % (time.time() - t0), flush=True)
    return res


# ------------------------------------------------------------------ 54A no-go
def part_A(n=65, V=0.2):
    c = n // 2; m = S.M; k0 = S.k_for_v(V, m, 0.0); psi, _, _ = S.packet_k(n, m, 0.0, (c, c, c), S.SIG, (k0, 0, 0))
    _, _, kap = packets(V, n)
    q = G.charge(psi, np.ones((n, n, n)), m, 0.0)
    for _ in range(3):
        N, info = S.moving([psi], [m], 0.0, kap, n, (c, c, c), V, qfix=q, tol=1e-11); q = G.charge(psi, N, m, 0.0)
    psiN = np.log(N); dpsi = np.zeros_like(psiN); dpsi[1:-1] = -V * (psiN[2:] - psiN[:-2]) / 2      # d_t psi = -V d_x psi
    b = np.stack([N * shp(N, a, 1) * (shp(dpsi, a, 1) - dpsi) for a in range(3)])
    B, _ = initial_shift(T53.link_P(psi), kap, 4.0, n, V)
    def curl_xy(F):            # plaquette at x in the (x, y) plane: F_x(x) + F_y(x+x) - F_x(x+y) - F_y(x)
        return F[0] + shp(F[1], 0, 1) - shp(F[0], 1, 1) - F[1]
    cb, cB = curl_xy(b), curl_xy(B); rows = []
    for rho in range(8, 25):
        rows.append(dict(rho=rho, curl_b=float(cb[c, c + rho, c]), curl_B=float(cB[c, c + rho, c]), ratio=float(abs(cb[c, c + rho, c]) / abs(cB[c, c + rho, c]))))
    # the same offset's lattice-gradient part: b - grad of (N^2-weighted) d_t psi is the nonlinear remainder
    res = dict(rows=rows, max_ratio=max(r["ratio"] for r in rows), b_max=float(np.abs(b).max()), B_max=float(np.abs(B).max()), info=str(info), kappa=kap)
    print("54A max curl ratio %.3e  |b|max %.3e |B|max %.3e" % (res["max_ratio"], res["b_max"], res["B_max"]), flush=True)
    return res


# ------------------------------------------------------------------ 54D smooth gauge
def part_D():
    L = 256; one = np.ones(L); s = np.arange(L); res = {}
    ba = 1e-3 * one; da = T53.loop_dt(one, ba); res["a_uniform"] = dict(dt=da, sum_abs=float(np.abs(ba).sum()))
    for Lam, A in ((16, 1e-3), (32, 1e-3), (64, 1e-3), (128, 1e-3), (64, 1e-4), (64, 1e-2)):
        chi = A * np.sin(2 * np.pi * s / Lam); bl = np.roll(chi, -1) - chi; d = T53.loop_dt(one, bl)
        rel = abs(d) / (abs(da) * np.abs(bl).sum() / np.abs(ba).sum())
        res[f"{Lam}_{A:g}"] = dict(Lam=Lam, A=A, dt=d, sum_abs=float(np.abs(bl).sum()), rel=float(rel), readout=T53.readout(one, bl)); print("54D", Lam, A, d, rel, flush=True)
    res["zero"] = T53.loop_dt(one, 0 * one)
    return res


def save(name, res, t0):
    res["_meta"] = dict(wall_s=time.time() - t0, platform=platform.platform(), python=platform.python_version(), numpy=np.__version__)
    (OUT / f"EX-GC54_{name}.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    job = sys.argv[1]; t0 = time.time()
    if job in ("A", "D", "ops", "ring"):
        save(job, {"A": part_A, "D": part_D, "ops": part_ops, "ring": part_ring}[job](), t0)
    else:                       # e.g. B_D2E_0.2  B_D2q_0.0  C_P_0.2_lam2  C1_P_0.2
        parts = job.split("_"); kind, rule, V = parts[0], parts[1], float(parts[2])
        lam = 2.0 if kind == "C" else 4.0
        T = 8.0 if kind in ("C1", "C2") else 16.0
        dt = 0.05 if kind == "C2" else 0.1
        save(job, run_pair(V, rule, lam=lam, T=T, dt=dt), t0)
    print(job, "done %.0fs" % (time.time() - t0), flush=True)
