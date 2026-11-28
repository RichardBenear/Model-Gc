"""EX-GC61: the shift as a travelling field. Local, order-free events carry beta (drift, two-link coincidence, source); continuum limit is
53's wave law, static limit is 59's Poisson law. 61A events vs RK4 PDE and order; 61B arrival, settling, boxes; 61C probe J and Sigma in the
travelling field with nothing refitted; 61D detector add-on (matched transformed readout, NB = 64, tilted spin vector prediction).
Plan: spec/EX-GC61_plan.md (amendment A1: scored runs use flat weights). EXPLORATORY. Our construction; not LaViolette's Model G and
not Chain Theory.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np
import numba

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc49 as F49
import exgc51 as G
import exgc52 as S
import exgc53 as T53
import exgc57 as E57
import exgc58 as X
import exgc59 as E
import exgc60 as P60

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC61"; PART = OUT / "partial"
LAM = 4.0; TAU = 8.0; DT = 0.1; I3 = (slice(1, -1),) * 3
K60 = dict(S=dict(rot=-0.9975751952997015, sh=0.0927428947179676), J=dict(rot=-0.9975414696755891, sh=0.0019794638314299963))


def ramp(t):
    return float(np.sin(np.pi * t / (2 * TAU)) ** 2) if t < TAU else 1.0


def faces0(a):
    nd = a.ndim
    for ax in range(nd - 3, nd):
        sl = [slice(None)] * nd; sl[ax] = 0; a[tuple(sl)] = 0; sl[ax] = -1; a[tuple(sl)] = 0
    return a


def lap(f):
    out = -6.0 * f[I3]
    for ax in range(3):
        for d in (-1, 1):
            sl = [slice(1, -1)] * 3; sl[ax] = slice(1 + d, f.shape[ax] - 1 + d); out = out + f[tuple(sl)]
    return out


class Field:
    """beta on links, rate pi; events per tick: drift (dt/2), coincidence + source (dt), drift (dt/2). Nb None = flat weights."""
    def __init__(self, j, lamk, Nb=None):
        self.j = j; self.lamk = lamk; self.Nb = Nb; self.B = np.zeros_like(j); self.Pi = np.zeros_like(j); self.t = 0.0

    def op(self, a, f):
        return lap(f) if self.Nb is None else F49.op(self.Nb[a], f)

    def w(self):
        return 1.0 if self.Nb is None else self.Nb

    def step(self, dt):
        self.B += 0.5 * dt * self.w() * self.Pi; faces0(self.B); r = ramp(self.t + 0.5 * dt)
        for a in range(3):
            wa = 1.0 if self.Nb is None else self.Nb[a][I3]
            self.Pi[a][I3] += dt * wa * (self.op(a, self.B[a]) - self.lamk * r * self.j[a][I3])
        self.B += 0.5 * dt * self.w() * self.Pi; faces0(self.B); self.t += dt

    def energy(self):
        Kin = 0.5 * float(np.sum(self.Pi ** 2)); V = 0.0
        for a in range(3):
            for ax in range(3):
                V += 0.5 * float(np.sum(np.diff(self.B[a], axis=ax) ** 2))
        Vs = self.lamk * ramp(self.t) * float(np.sum(self.j * self.B)); return Kin + V + Vs, Kin + V + abs(Vs)


def rk4_field(j, lamk, T, dt):
    """53's evolve_b shift right-hand side (flat weights), matter held fixed."""
    B = np.zeros_like(j); Pi = np.zeros_like(j); t = 0.0
    def rhs(B, Pi, t):
        dB = faces0(Pi.copy()); dP = np.zeros_like(Pi); r = ramp(t)
        for a in range(3):
            dP[a][I3] = lap(B[a]) - lamk * r * j[a][I3]
        return dB, dP
    while t < T - 1e-9:
        k1 = rhs(B, Pi, t); k2 = rhs(B + dt / 2 * k1[0], Pi + dt / 2 * k1[1], t + dt / 2); k3 = rhs(B + dt / 2 * k2[0], Pi + dt / 2 * k2[1], t + dt / 2)
        k4 = rhs(B + dt * k3[0], Pi + dt * k3[1], t + dt)
        B = B + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); Pi = Pi + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]); t += dt
    return B


@numba.njit(cache=True)
def _kick_events(B, Pi, j, lamk_r, dt, perm, n):
    n3 = n * n * n
    for e in perm:
        if e < 9 * n3:
            a = e // (3 * n3); b = (e // n3) % 3; xl = e % n3; x = xl // (n * n); y = (xl // n) % n; z = xl % n
            x2, y2, z2 = x, y, z
            if b == 0:
                x2 += 1
            elif b == 1:
                y2 += 1
            else:
                z2 += 1
            if x2 >= n or y2 >= n or z2 >= n:
                continue
            d = B[a, x2, y2, z2] - B[a, x, y, z]
            Pi[a, x, y, z] += dt * d; Pi[a, x2, y2, z2] -= dt * d
        else:
            k = e - 9 * n3; a = k // n3; xl = k % n3; x = xl // (n * n); y = (xl // n) % n; z = xl % n
            Pi[a, x, y, z] -= dt * lamk_r * j[a, x, y, z]


def step_shuffled(F, dt, rng):
    n = F.B.shape[1]; F.B += 0.5 * dt * F.Pi; faces0(F.B); perm = rng.permutation(12 * n ** 3)
    _kick_events(F.B, F.Pi, F.j, F.lamk * ramp(F.t + 0.5 * dt), dt, perm, n); faces0(F.Pi)
    F.B += 0.5 * dt * F.Pi; faces0(F.B); F.t += dt


def embed(a, n_big):
    n = a.shape[-1]; o = (n_big - n) // 2; out = np.zeros(a.shape[:-3] + (n_big,) * 3); out[..., o:o + n, o:o + n, o:o + n] = a; return out


def site_omz(B, x):
    cr = B[:, x[0] - 2:x[0] + 3, x[1] - 2:x[1] + 3, x[2] - 2:x[2] + 3]; return -0.5 * float(X.curl_op(cr)[2][2, 2, 2])


def sites(n):
    c = n // 2; out = {}
    for d in (8, 12, 16, 24):
        out[f"axis_{d}"] = (c, c, c + d); out[f"plane_{d}"] = (c + d, c, c)
    return out


def t_half(ts, om, ref):
    om = np.asarray(om); k = np.argmax(om / ref >= 0.5); f = (0.5 * ref - om[k - 1]) / (om[k] - om[k - 1]); return float(ts[k - 1] + f * (ts[k] - ts[k - 1]))


# ------------------------------------------------------------------ 61A
def part_A():
    n = 65; src, kap, N, Nd = E.source_terms(n); j = src["up"]["PS"]; lk = LAM * kap; res = {}
    for dt in (0.1, 0.05):
        F = Field(j, lk)
        for _ in range(int(round(20 / dt))):
            F.step(dt)
        Bref = rk4_field(j, lk, 20.0, dt); res[f"diff_dt{dt}"] = float(np.abs(F.B - Bref).max() / np.abs(Bref).max())
        if dt == 0.1:
            Bv = F.B.copy()
    res["ratio"] = res["diff_dt0.1"] / res["diff_dt0.05"]
    for seed in (1, 2):
        F = Field(j, lk); rng = np.random.default_rng(seed)
        for _ in range(200):
            step_shuffled(F, 0.1, rng)
        res[f"shuffle_seed{seed}"] = float(np.abs(F.B - Bv).max() / np.abs(Bv).max())
    F = Field(np.zeros_like(j), lk)
    for _ in range(50):
        F.step(0.1)
    res["C1_zero"] = float(np.abs(F.B).max()); print("61A", res, flush=True); return res


# ------------------------------------------------------------------ 61B
def odd_source(n_box):
    src, kap, N, Nd = E.source_terms(97); jo = 0.5 * (src["up"]["PS"] - src["down"]["PS"])
    return (jo if n_box == 97 else embed(jo, n_box)), kap, N, src


def part_B(n_box, weighted=False, Tend=100.0):
    jo, kap, N, src = odd_source(n_box); lk = LAM * kap; Nb = None
    if weighted:
        Nb = np.stack([0.5 * (N + np.roll(N, -1, a)) for a in range(3)])
    F = Field(jo, lk, Nb); st = sites(n_box); rec = {k: [] for k in st}; ts = []; En = []
    def record():
        ts.append(F.t)
        for k, x in st.items():
            rec[k].append(site_omz(F.B, x))
    record(); B50 = None
    while F.t < Tend - 1e-9:
        F.step(DT); record()
        if abs(F.t - 50.0) < 1e-6:
            B50 = F.B.copy()
        if not weighted and n_box == 97 and abs(F.t / 10 - round(F.t / 10)) < 1e-6 and F.t > 9.99:
            En.append(F.energy())
    ts = np.array(ts); i50 = int(np.argmin(np.abs(ts - 50))); res = dict(n=n_box, weighted=weighted, ts_dt=DT, om=rec)
    res["t_half"] = {k: t_half(ts, rec[k], rec[k][i50]) for k in st}
    if n_box == 97:
        PART.mkdir(parents=True, exist_ok=True); np.save(PART / f"B50_{'w' if weighted else 'flat'}.npy", B50)
        j = src["up"]["PS"]; div = sum(j[a] - np.roll(j[a], 1, a) for a in range(3)); res["C4_div_rel"] = float(np.abs(div).max() / np.abs(j).max())
        # settled field vs static: flat -> n = 193 zero-box Poisson; weighted -> its own weighted static (CG) in n = 97
        LP = T53.link_pos(97); R = [np.linalg.norm(LP[a], axis=0) for a in range(3)]; msk = [(R[a] >= 8) & (R[a] <= 20) for a in range(3)]
        if not weighted:
            j193 = embed(jo, 193); Bs = np.stack([T53.poisson(lk * j193[a], np.zeros_like(j193[a])) for a in range(3)])[:, 48:145, 48:145, 48:145]
            Bs97 = np.stack([T53.poisson(lk * jo[a], np.zeros_like(jo[a])) for a in range(3)])
            res["rms_vs_static193"] = float(np.sqrt(sum(np.sum((B50[a] - Bs[a])[msk[a]] ** 2) for a in range(3)) / sum(np.sum(Bs[a][msk[a]] ** 2) for a in range(3))))
            res["rms_vs_static97"] = float(np.sqrt(sum(np.sum((B50[a] - Bs97[a])[msk[a]] ** 2) for a in range(3)) / sum(np.sum(Bs97[a][msk[a]] ** 2) for a in range(3))))
            Es = np.array(En); res["C2_energy"] = dict(E=[e[0] for e in En], spread_rel=float((Es[:, 0].max() - Es[:, 0].min()) / Es[:, 1].max()))
        else:
            from scipy.sparse.linalg import LinearOperator, cg
            m = 95; Bw = np.zeros_like(jo)
            for a in range(3):
                def mv(v, a=a):
                    f = np.zeros((97,) * 3); f[I3] = v.reshape((m,) * 3); return -F49.op(Nb[a], f).ravel()
                Aop = LinearOperator((m ** 3, m ** 3), matvec=mv); sol, info = cg(Aop, -(lk * jo[a][I3]).ravel(), rtol=1e-10, maxiter=3000)
                Bw[a][I3] = sol.reshape((m,) * 3)
            res["rms_vs_own_static"] = float(np.sqrt(sum(np.sum((B50[a] - Bw[a])[msk[a]] ** 2) for a in range(3)) / sum(np.sum(Bw[a][msk[a]] ** 2) for a in range(3))))
    print("61B", n_box, weighted, json.dumps(res["t_half"]), {k: v for k, v in res.items() if k.startswith(("rms", "C"))}, flush=True); return res


# ------------------------------------------------------------------ probes
def vecs(psi):
    s = np.array([float(np.real(np.vdot(psi, psi @ X.SIG[k].T))) for k in range(3)]); L = np.array(T53.angmom(T53.link_P(psi)), float)
    return s, L, L + 0.5 * s


class Probe:
    def __init__(self, Ns, psi):
        self.Ns = Ns; self.psi = psi.copy(); self.rec = []; self.add()

    def H(self, p, B, C):
        return G.Hm(p, self.Ns, P60.M, 0.0) + T53.Hb(p, B) + X.Hpl(p, C)

    def step(self, B0, B1, dt):
        C0 = X.curl_op(B0); C1 = X.curl_op(B1); Bm = 0.5 * (B0 + B1); Cm = 0.5 * (C0 + C1); p = self.psi
        k1 = -1j * self.H(p, B0, C0); k2 = -1j * self.H(p + dt / 2 * k1, Bm, Cm); k3 = -1j * self.H(p + dt / 2 * k2, Bm, Cm); k4 = -1j * self.H(p + dt * k3, B1, C1)
        self.psi = p + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); self.add()

    def add(self):
        s, L, J = vecs(self.psi); self.rec.append(np.concatenate([s, L, J]).tolist())


def ang(rec, off):
    r = np.array(rec); return np.arctan2(r[:, off + 1], r[:, off])


def odd_rate(ru, rd, off, dt=DT):
    return P60.slope((ang(ru, off) - ang(rd, off)) / 2, dt)


def wgrad(Bs, sig=3.0, chi="x"):
    """probe-density-weighted lattice gradient matrix M[a][b] = d_b beta_a on the 48^3 sub-box (60's stencils for a != b)."""
    p = P60.probe(sig, chi) if isinstance(chi, str) else S.packet_k(Bs.shape[1], P60.M, 0.0, (Bs.shape[1] // 2,) * 3, sig, (0, 0, 0), chi=chi)[0]
    w = np.sum(np.abs(p) ** 2, -1); w /= w.sum(); M = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            if a == b:
                g = Bs[a] - np.roll(Bs[a], 1, a)
            else:
                q = Bs[a]; g = 0.25 * (np.roll(q, -1, b) - np.roll(q, 1, b) + np.roll(np.roll(q, -1, b), 1, a) - np.roll(np.roll(q, 1, b), 1, a))
            M[a, b] = float(np.sum(w * g))
    return M


def pred_xy(M):
    sa, sb = M[0, 1], M[1, 0]; om = (sb - sa) / 2; sh = (sb + sa) / 2
    return {k: K60[k]["rot"] * om + K60[k]["sh"] * sh for k in ("S", "J")}, {k: K60[k]["rot"] * om for k in ("S", "J")}


def part_C(n_box, spin, windows):
    """field for one source spin in an n_box box; probes in the travelling field at the axis and diagonal sites over the windows."""
    src, kap, N, Nd = E.source_terms(97); j = src[{"dn": "down"}.get(spin, spin)]["PS"]; j = j if n_box == 97 else embed(j, n_box); o = (n_box - 97) // 2; c = 48
    F = Field(j, LAM * kap); xs97 = {"axis_16": (c, c, c + 16), "diag": (c + 12, c, c + 12)}; res = dict(n=n_box, spin=spin, probes={}, grad45={})
    jobs = [(k, t0) for (k, t0) in windows]; live = {}
    while F.t < max(t0 for _, t0 in jobs) + 10.0 - 1e-9:
        for (k, t0) in jobs:
            if abs(F.t - t0) < 1e-6:
                live[(k, t0)] = Probe(E.cut_site(N, xs97[k]), P60.probe(3.0, "x"))
        olds = {key: E57.sub(F.B, tuple(xs97[key[0]][i] + o for i in range(3))) for key in live}
        F.step(DT)
        for key, pr in list(live.items()):
            pr.step(olds[key], E57.sub(F.B, tuple(xs97[key[0]][i] + o for i in range(3))), DT)
            if abs(F.t - (key[1] + 10.0)) < 1e-6:
                res["probes"][f"{key[0]}@{key[1]:g}"] = pr.rec; del live[key]
        if abs(F.t - 45.0) < 1e-6:
            for k, x in xs97.items():
                np.save(PART / f"B45_{n_box}_{spin}_{k}.npy", E57.sub(F.B, tuple(x[i] + o for i in range(3))))
    print("61C", n_box, spin, "done", flush=True); return res


def part_C3():
    """C3: 61C probe code with 59's static zero-box field at axis 16 reproduces 59's Sigma rate."""
    src, kap, N, Nd = E.source_terms(97); c = 48; x0 = (c, c, c + 16); Ns = E.cut_site(N, x0); out = {}
    for sp in ("up", "down"):
        Bs = E57.sub(E.shift(src[sp]["PS"], kap, "zero"), x0); pr = Probe(Ns, P60.probe(3.0, "x"))
        for _ in range(100):
            pr.step(Bs, Bs, DT)
        out[sp] = pr.rec
    om = odd_rate(out["up"], out["down"], 0); ref = json.loads((ROOT / "results/EX-GC59/EX-GC59_main.json").read_text())["runs"]["axis_16"]
    res = dict(Omega=om, ref=ref["Omega_meas"], rel=om / ref["Omega_meas"] - 1, Omega_pred_59=ref["Omega_pred"],
               Omega_gr=E.omega_pred(src["up"]["PS"], kap, 97, x0)[0])
    print("C3", res, flush=True); return res


# ------------------------------------------------------------------ 61D
NT = 6


def gauge_ops(chi_s, NBx):
    one = np.ones((NBx,) * 3); x = np.arange(NBx) - NBx // 2; Xg, Yg = np.meshgrid(x, x, indexing="ij")
    chi = (-chi_s * Xg * Yg)[..., None][:, :, None, :] * np.ones((1, 1, NBx, 1)); Hm = lambda p: G.Hm(p, one, P60.M, 0.0)
    Hg = lambda p: -0.5j * (chi * Hm(Hm(p)) - Hm(Hm(chi * p))); Gop = lambda p: 0.5 * (chi * Hm(p) + Hm(chi * p))
    def expG(p, sgn):
        out = p.copy(); term = p.copy()
        for k in range(1, NT + 1):
            term = (sgn * -1j) / k * Gop(term); out = out + term
        return out
    return Hg, expG


def run_static(Ns, Bs, rule, psi, extra=None, tr=None, dt=DT, T=10.0):
    C = X.curl_op(Bs) if (Bs is not None and rule == "PS") else None; raw = []; trn = []
    def H(p):
        o = G.Hm(p, Ns, P60.M, 0.0)
        if Bs is not None:
            o = o + T53.Hb(p, Bs) + (X.Hpl(p, C) if C is not None else 0)
        return o + (extra(p) if extra is not None else 0)
    def rec(p):
        raw.append(np.concatenate(vecs(p)).tolist())
        if tr is not None:
            trn.append(np.concatenate(vecs(tr(p, -1))).tolist())
    rec(psi); t = 0.0
    while t < T - 1e-9:
        k1 = -1j * H(psi); k2 = -1j * H(psi + dt / 2 * k1); k3 = -1j * H(psi + dt / 2 * k2); k4 = -1j * H(psi + dt * k3)
        psi = psi + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt; rec(psi)
    return dict(raw=raw, tr=trn)


def part_D1(sig):
    one = np.ones((48,) * 3); s0 = P60.S0; psi = P60.probe(sig, "x"); out = dict(sigma=sig)
    Hg, expG = gauge_ops(s0, 48); out["free"] = run_static(one, None, "PS", psi)["raw"]
    out["Hg"] = run_static(one, None, "PS", psi, extra=Hg, tr=expG)
    out["P"] = run_static(one, P60.field("shear", "xy"), "PS", psi, tr=expG)
    if sig == 3.0:
        out["sum"] = run_static(one, P60.field("sum", "xy"), "PS", psi, tr=expG)
        out["rot"] = run_static(one, P60.field("rot", "xy"), "PS", psi)["raw"]
        Hg2, expG2 = gauge_ops(2 * s0, 48); out["Hg2"] = run_static(one, None, "PS", psi, extra=Hg2, tr=expG2)
    print("61D1", sig, "done", flush=True); return out


def part_D2():
    P60.NB = 64; one = np.ones((64,) * 3); r0 = P60.run(one, None, "PS", 5.0, "xy"); out = {}
    for rule in ("P", "PS"):
        out[rule] = P60.rates(P60.run(one, P60.field("shear", "xy"), rule, 5.0, "xy"), r0, P60.S0)
    print("61D2", out, flush=True); return out


GMAT = np.array([[0.3, -0.7, 0.4], [0.5, 0.2, -0.6], [-0.1, 0.8, -0.5]])
TH = np.arccos(1 / np.sqrt(3)); CHI_N = (np.cos(TH / 2), np.exp(1j * np.pi / 4) * np.sin(TH / 2), 0, 0)


def part_D3():
    out = {}; one = np.ones((48,) * 3); psi = S.packet_k(48, P60.M, 0.0, (24, 24, 24), 3.0, (0, 0, 0), chi=CHI_N)[0]
    LP = T53.link_pos(48); B = np.stack([P60.S0 * sum(GMAT[a, b] * LP[a][b] for b in range(3)) for a in range(3)])
    out["flat_free"] = run_static(one, None, "PS", psi)["raw"]; out["flat_G"] = run_static(one, B, "PS", psi)["raw"]
    src, kap, N, Nd = E.source_terms(97); c = 48; x0 = (c + 9, c + 5, c + 13); Ns = E.cut_site(N, x0); Bsub = {}
    for sp in ("up", "down"):
        Bsub[sp] = E57.sub(E.shift(src[sp]["PS"], kap, "zero"), x0); out[f"src_{sp}"] = run_static(Ns, Bsub[sp], "PS", psi)["raw"]
    out["M_src"] = wgrad(0.5 * (Bsub["up"] - Bsub["down"]), chi=CHI_N).tolist(); print("61D3 done", flush=True); return out


def save(name, res, t0):
    res["_meta"] = dict(wall_s=time.time() - t0, platform=platform.platform(), python=platform.python_version(), numpy=np.__version__)
    (OUT / f"EX-GC61_{name}.json").write_text(json.dumps(res, default=float))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True); PART.mkdir(parents=True, exist_ok=True); job = sys.argv[1]; t0 = time.time()
    if job == "A":
        r = part_A()
    elif job == "B97":
        r = part_B(97)
    elif job == "B97w":
        r = part_B(97, weighted=True, Tend=60.0)
    elif job == "B129":
        r = part_B(129)
    elif job.startswith("C97") or job.startswith("C129"):
        nb = 97 if job.startswith("C97") else 129; sp = job[len(f"C{nb}"):]
        w = [("axis_16", 40.0), ("diag", 40.0)] + ([("axis_16", 90.0)] if nb == 97 else [])
        r = part_C(nb, sp, w)
    elif job == "C3":
        r = part_C3()
    elif job.startswith("D1_"):
        r = part_D1(float(job[3:]))
    elif job == "D2":
        r = part_D2()
    else:
        r = part_D3()
    save(job, r, t0); print("done", job, "%.0fs" % (time.time() - t0), flush=True)
