"""EX-GC55: one interaction for the link shift, and a gauge-covariant clock.
55A covariant light rule on the one-clock loop (Hamiltonian ADM form: u_t = N pi + beta d_s u, pi_t = lap u + d_s(beta pi)),
    smooth-gauge / uniform / one-sided controls and the 53C/D rotating ring;
55B energy bookkeeping (N = 1) for P/P, mixed J_E/P (54) and reciprocal J_E/J_E couplings;
55C reciprocal energy-current rule (RJE) in the co-moving pair vs P/P at n = 49, T = 8.
Plan: spec/EX-GC55_plan.md. EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G
import exgc52 as S
import exgc53 as T53
import exgc54 as X54
import exgc49 as F49

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC55"
Nf, shp = G.Nf, G.shp
LAM = 4.0; K0 = 0.5; LF = float(np.cos(K0) / np.cos(K0 / 2) ** 2)


# ------------------------------------------------------------------ 55A covariant loop clock
def loop_wave(Ns, bl, rule="cov", init="new", k0=K0, width=16.0, dt=0.05, laps=1.0, direction=+1):
    L = len(Ns); s = np.arange(L); bs = 0.5 * (bl + np.roll(bl, 1)); Nl = Ns * np.roll(Ns, -1)
    grad = lambda a: 0.5 * (np.roll(a, -1) - np.roll(a, 1))
    def rhs(u, w):
        lap = Nl * (np.roll(u, -1) - u) + np.roll(Nl, 1) * (np.roll(u, 1) - u)
        return Ns * w + bs * grad(u), Ns * lap + (grad(bs * w) if rule == "cov" else bs * grad(w))
    d = (s - L // 2 + L) % L - L // 2; env = np.exp(-d ** 2 / (2 * width ** 2)); kk = direction * k0; om = 2 * np.sin(k0 / 2) * Ns[0] ** 2
    u = env * np.cos(kk * d); ut = om * env * np.sin(kk * d)
    w = (ut - bs * grad(u)) / Ns if init == "new" else ut / Ns
    ang = 2 * np.pi * s / L; t = 0.0; unw = 0.0; target = laps * L
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
    return None


def loop_dt(Ns, bl, **kw):
    return loop_wave(Ns, bl, direction=+1, **kw) - loop_wave(Ns, bl, direction=-1, **kw)


def part_A(ring=True):
    L = 256; one = np.ones(L); s = np.arange(L); res = {"lattice_factor": LF}
    ba = 1e-3 * one; da = loop_dt(one, ba); res["a_uniform"] = dict(dt=da, readout=T53.readout(one, ba), ratio=da / T53.readout(one, ba))
    bc = np.where(s < 64, 1e-3, 0.0); dc = loop_dt(one, bc); res["c_one_side"] = dict(dt=dc, ratio_to_a=dc / da)
    res["zero"] = dict(cov=loop_dt(one, 0 * one), adv=T53.loop_dt(one, 0 * one))
    cases = [(L_, 1e-3, 0.0) for L_ in (16, 32, 64, 128)] + [(64, 1e-3, ph) for ph in (np.pi / 4, np.pi / 2, 3 * np.pi / 4)] + [(64, 1e-4, 0.0), (64, 1e-2, 0.0)]
    g = {}
    for Lam, A, ph in cases:
        chi = A * np.sin(2 * np.pi * s / Lam + ph); bl = np.roll(chi, -1) - chi; norm = abs(da) * np.abs(bl).sum() / np.abs(ba).sum()
        row = dict(Lam=Lam, A=A, phi=ph, beta_clock=float(bl[0]), dt=loop_dt(one, bl))
        row["rel"] = abs(row["dt"]) / norm
        row["diag_cov_oldinit"] = loop_dt(one, bl, init="old") / norm; row["diag_adv_newinit"] = loop_dt(one, bl, rule="adv") / norm
        g[f"{Lam}_{A:g}_{ph:.3f}"] = row; print("55A", Lam, A, round(ph, 3), "rel %.2e  covOld %.2e advNew %.2e" % (row["rel"], row["diag_cov_oldinit"], row["diag_adv_newinit"]), flush=True)
    res["gauge"] = g; print("55A uniform ratio %.5f one-side %.4f" % (res["a_uniform"]["ratio"], res["c_one_side"]["ratio_to_a"]), flush=True)
    if ring:
        res["ring"] = part_ring()
    return res


def part_ring(n=129, u=0.01):
    c = n // 2; q0, P0 = T53.ring_sources(n, 0); kap = 4 * np.pi * 12 * u / float(q0.sum()); out = {}
    for ell in (2, -2):
        t0 = time.time(); q, P = T53.ring_sources(n, ell); J = T53.angmom(P); N, info = G.static([], [], 0.0, kap, n, (c, c, c), qfix=q); B = T53.shift_static(P, kap, J)
        for h in (20, 28, 36):
            Ns, bl = T53.loop_fields(T53.square_loop((c, c, c), h), N, B); ro = T53.readout(Ns, bl); w = loop_dt(Ns, bl)
            out[f"{ell}_{h}"] = dict(readout=ro, wave=w, ratio=w / ro, ratio_over_LF=w / ro / LF)
            print("55A ring ell %d h %d readout %.6e wave %.6e ratio %.5f (/LF %.5f)" % (ell, h, ro, w, w / ro, w / ro / LF), flush=True)
        print("55A ring ell", ell, "%.0fs" % (time.time() - t0), flush=True)
    return out


# ------------------------------------------------------------------ reciprocal energy-current coupling
def site_b(B):
    return [0.5 * (B[a] + shp(B[a], a, -1)) for a in range(3)]


def Oa(phi, N, m, a):
    return 0.5 * (X54.Hop(X54.vop(phi, N, a), N, m) + X54.vop(X54.Hop(phi, N, m), N, a))


def je_terms(psi, N, m, B):
    """returns (H_beta psi, link source s) for the J_E/J_E coupling H_beta = -sum_a (1/2){b_a, O_a}."""
    b = site_b(B); Hb = np.zeros_like(psi); s = np.zeros((3,) + N.shape)
    for a in range(3):
        Op = Oa(psi, N, m, a); d = np.real(np.sum(psi.conj() * Op, -1)); s[a] = 0.5 * (d + shp(d, a, 1))
        Hb -= 0.5 * (b[a][..., None] * Op + Oa(b[a][..., None] * psi, N, m, a))
    return Hb, s


def coupling(psi, N, m, B, rule):
    """(H_beta psi, source current, coupling current) for rule in P (P/P), MIX (J_E source, P coupling), RJE (J_E/J_E)."""
    if rule == "P":
        P = T53.link_P(psi); return T53.Hb(psi, B), P, P
    if rule == "MIX":
        P = T53.link_P(psi); return T53.Hb(psi, B), X54.rule_current(psi, N, m, "D2E"), P
    Hb, s = je_terms(psi, N, m, B); return Hb, s, s


def faces0(a):
    for ax in range(3):
        sl = [slice(None)] * 3; sl[ax] = 0; a[tuple(sl)] = 0; sl[ax] = -1; a[tuple(sl)] = 0
    return a


def evolve(psis, ms, f, g, B, gB, kappa, lam, T, dt, rule, frozenN=False, energy=False):
    I = (slice(1, -1),) * 3
    def rhs(ps, f, g, B, gB):
        N = Nf(f); terms = [coupling(p, N, mm, B, rule) for p, mm in zip(ps, ms)]
        dps = [-1j * (G.Hm(p, N, mm, 0.0) + tb[0]) for p, mm, tb in zip(ps, ms, terms)]
        j = sum(tb[1] for tb in terms)
        if frozenN:
            df = np.zeros_like(f); dg = np.zeros_like(g)
        else:
            q = sum(G.charge(p, N, mm, 0.0) for p, mm in zip(ps, ms)); df = faces0(N * g); dg = np.zeros_like(g); dg[I] = N[I] * (F49.op(N, f) - kappa * q[I])
        dB = np.zeros_like(B); dgB = np.zeros_like(gB)
        for a in range(3):
            Nb = 0.5 * (N + np.roll(N, -1, a)); dB[a] = faces0(Nb * gB[a]); dgB[a][I] = Nb[I] * (F49.op(Nb, B[a]) - lam * kappa * j[a][I])
        return dps, df, dg, dB, dgB
    def add(st, k, h):
        return [p + h * kp for p, kp in zip(st[0], k[0])], st[1] + h * k[1], st[2] + h * k[2], st[3] + h * k[3], st[4] + h * k[4]
    def rec(st, t):
        r = dict(t=t, Py=[S.momentum(p, 1) for p in st[0]])
        if energy:
            r.update(Etot(st, ms, kappa, lam, rule))
        return r
    st = (psis, f, g, B, gB); t = 0.0; hist = [rec(st, 0.0)]
    while t < T - 1e-9:
        k1 = rhs(*st); k2 = rhs(*add(st, k1, dt / 2)); k3 = rhs(*add(st, k2, dt / 2)); k4 = rhs(*add(st, k3, dt))
        st = ([p + dt / 6 * (a + 2 * b + 2 * c + d) for p, a, b, c, d in zip(st[0], k1[0], k2[0], k3[0], k4[0])],) + tuple(
            st[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(1, 5))
        t += dt; hist.append(rec(st, t))
    return st, hist


def Etot(st, ms, kappa, lam, rule):
    """N = 1: E = <H_m> - sum beta j_coupling (interior) - (1/2 lam kappa) sum[(d_t beta)^2 + (grad beta)^2]."""
    ps, f, g, B, gB = st; N = Nf(f); I = (slice(1, -1),) * 3
    Em = sum(float(np.real(np.vdot(p, G.Hm(p, N, mm, 0.0)))) for p, mm in zip(ps, ms))
    jc = sum(coupling(p, N, mm, B, rule)[2] for p, mm in zip(ps, ms)); Ei = -float(np.sum(B[(slice(None),) + I] * jc[(slice(None),) + I]))
    kin = float(np.sum(gB[(slice(None),) + I] ** 2)); grad2 = 0.0
    for a in range(3):
        for ax in range(3):
            dB = np.diff(B[a], axis=ax); grad2 += float(np.sum(dB ** 2))
    Eb = -(kin + grad2) / (2 * lam * kappa)
    return dict(Em=Em, Ei=Ei, Eb=Eb, E=Em + Ei + Eb)


# ------------------------------------------------------------------ 55B bookkeeping
def part_B(rule, dt=0.1, n=49, V=0.2, T=8.0):
    c = n // 2; m = S.M; k0 = S.k_for_v(V, m, 0.0); psi, _, _ = S.packet_k(n, m, 0.0, (c, c, c), S.SIG, (k0, 0, 0))
    _, _, kap = X54.packets(V, 65)                                    # kappa as in 53E (n = 65 normalisation)
    one = np.ones((n, n, n)); f = G.fN(one); g = np.zeros_like(f)
    j0 = coupling(psi, one, m, np.zeros((3, n, n, n)), rule)[1]; B, gB = X54.initial_shift(j0, kap, LAM, n, V)
    t0 = time.time(); _, hist = evolve([psi], [m], f, g, B, gB, kap, LAM, T, dt, rule, frozenN=True, energy=True)
    E = np.array([h["E"] for h in hist]); Ei0 = hist[0]["Ei"]
    res = dict(rule=rule, dt=dt, n=n, V=V, T=T, kappa=kap, hist=hist, drift=float(np.max(np.abs(E - E[0]))), Ei0=Ei0, drift_rel=float(np.max(np.abs(E - E[0])) / abs(Ei0)), wall=time.time() - t0)
    print("55B %s dt %.3f drift %.3e  |Ei0| %.3e  rel %.3e  %.0fs" % (rule, dt, res["drift"], abs(Ei0), res["drift_rel"], res["wall"]), flush=True)
    return res


# ------------------------------------------------------------------ 55C pair
def run_pair(V, rule, n=49, sep=14, T=8.0, dt=0.1, lam=LAM):
    c = n // 2; m = S.M; k0 = S.k_for_v(V, m, 0.0)
    p1, _, _ = S.packet_k(n, m, 0.0, (c, c - sep // 2, c), S.SIG, (k0, 0, 0)); p2, _, _ = S.packet_k(n, m, 0.0, (c, c + sep // 2, c), S.SIG, (k0, 0, 0))
    p0, _, _ = S.packet_k(n, m, 0.0, (c, c - sep // 2, c), S.SIG, (0, 0, 0)); kap = 4 * np.pi * sep * 0.01 / float(G.charge(p0, np.ones((n, n, n)), m, 0.0).sum())
    out = dict(V=V, rule=rule, n=n, T=T, dt=dt, kappa=kap, runs={})
    for lab, ps in (("pair", [p1, p2]), ("one", [p1])):
        t0 = time.time(); f, g, info = S.initial_field(ps, [m] * len(ps), 0.0, kap, n, c, V)
        N = Nf(f); j0 = sum(coupling(p, N, m, np.zeros((3, n, n, n)), rule)[1] for p in ps); B, gB = X54.initial_shift(j0, kap, lam, n, V)
        _, hist = evolve(ps, [m] * len(ps), f, g, B, gB, kap, lam, T, dt, rule)
        out["runs"][lab] = dict(t=[h["t"] for h in hist], Py=[h["Py"] for h in hist], wall=time.time() - t0, j0_tot=j0.sum(axis=(1, 2, 3)).tolist())
        print("55C V %.1f %s %s %.0fs" % (V, rule, lab, time.time() - t0), flush=True)
    t = np.array(out["runs"]["pair"]["t"]); P1 = np.array([p[0] for p in out["runs"]["pair"]["Py"]]) - np.array([p[0] for p in out["runs"]["one"]["Py"]])
    out["F1"] = float(np.polyfit(t, P1, 2)[1]); print("55C V %.1f %s F1 %.6e" % (V, rule, out["F1"]), flush=True)
    return out


# ------------------------------------------------------------------ C1, C3
def part_checks():
    res = {}; L = 256; one = np.ones(L)
    res["C1"] = dict(cov=loop_dt(one, 0 * one), adv=T53.loop_dt(one, 0 * one))
    res["C1"]["diff"] = abs(res["C1"]["cov"] - res["C1"]["adv"])
    rng = np.random.default_rng(55); n = 13; m = 0.6; N = 1 + 0.02 * rng.standard_normal((n, n, n)); B = 1e-3 * rng.standard_normal((3, n, n, n))
    ph = rng.standard_normal((n, n, n, 4)) + 1j * rng.standard_normal((n, n, n, 4)); ps = rng.standard_normal((n, n, n, 4)) + 1j * rng.standard_normal((n, n, n, 4))
    Hph, _ = je_terms(ph, N, m, B); Hps, s = je_terms(ps, N, m, B)
    res["C3"] = dict(herm=float(abs(np.vdot(ph, Hps) - np.vdot(Hph, ps)) / abs(np.vdot(ph, Hps))),
                     expect=float(abs(np.real(np.vdot(ps, Hps)) + np.sum(B * s)) / abs(np.sum(B * s))))
    print("55 checks", res, flush=True)
    return res


def save(name, res, t0):
    res["_meta"] = dict(wall_s=time.time() - t0, platform=platform.platform(), python=platform.python_version(), numpy=np.__version__)
    (OUT / f"EX-GC55_{name}.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True); job = sys.argv[1]; t0 = time.time()
    if job == "A":
        r = part_A()
    elif job == "checks":
        r = part_checks()
    elif job.startswith("B_"):                     # B_P, B_MIX, B_RJE, B_P_half
        p = job.split("_"); r = part_B(p[1], dt=0.05 if len(p) > 2 else 0.1)
    else:                                           # C_RJE_0.2, C_P_0.0
        p = job.split("_"); r = run_pair(float(p[2]), p[1])
    save(job, r, t0); print(job, "done %.0fs" % (time.time() - t0), flush=True)
