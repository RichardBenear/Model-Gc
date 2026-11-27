"""EX-GC53: a link shift beta (GR's g_0i) added to the 51 model. 53A one-clock loop readout of an imposed shift; 53B matter
coupling H_beta = -sum beta_{x,a} Im(psi_x^+ psi_{x+a}); 53C/D rotating spin-balanced Dirac ring and stationary light loops;
53E 52D rerun with a dynamic shift (R -> 1?) and the strength lambda* Lorentz invariance needs.
Plan: spec/EX-GC53_plan.md (amendment A1). EXPLORATORY. Our construction; not LaViolette's Model G.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np
from scipy.fft import dstn, idstn

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G
import exgc52 as S
import exgc49 as F49

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC53"
fN, Nf, shp = G.fN, G.Nf, G.shp
LAM = 4.0
UP, DOWN = (1, 0, 0, 0), (0, 1, 0, 0)


# ------------------------------------------------------------------ matter coupling
def Hb(psi, B):
    """H_beta psi = -sum_a [B_a(x) psi(x+a) - B_a(x-a) psi(x-a)] / (2i); B[a] lives on link (x, x+a)."""
    out = np.zeros_like(psi)
    for a in range(3):
        out -= (B[a][..., None] * shp(psi, a, 1) - shp(B[a], a, -1)[..., None] * shp(psi, a, -1)) / 2j
    return out


def link_P(psi):
    """P_{x,a} = Im(psi_x^+ psi_{x+a}) (zero on the last layer)."""
    return np.stack([np.imag(np.sum(psi.conj() * shp(psi, a, 1), -1)) for a in range(3)])


# ------------------------------------------------------------------ 53A: 1D loop light with shift
def loop_wave(Ns, bl, k0=0.5, width=16.0, dt=0.05, laps=1.0, direction=+1, density="u2"):
    """1D periodic chain: (1/N) D_t[(1/N) D_t u] = sum N_s N_s'(u_s' - u_s), D_t = d_t - beta d_s.
    Ns: site N (L,); bl: link shift on (s, s+1) (L,). Returns the time for the energy centroid to travel laps*L in `direction`."""
    L = len(Ns); s = np.arange(L); bs = 0.5 * (bl + np.roll(bl, 1))           # beta at sites
    Nl = Ns * np.roll(Ns, -1)                                                    # N_s N_{s+1}
    def grad(a):
        return 0.5 * (np.roll(a, -1) - np.roll(a, 1))
    def rhs(u, w):
        lap = Nl * (np.roll(u, -1) - u) + np.roll(Nl, 1) * (np.roll(u, 1) - u)
        return Ns * w + bs * grad(u), Ns * lap + bs * grad(w)
    d = (s - L // 2 + L) % L - L // 2                                            # distance from site 0
    env = np.exp(-d ** 2 / (2 * width ** 2)); kk = direction * k0
    om = 2 * np.sin(k0 / 2) * Ns[0] ** 2
    u = env * np.cos(kk * d); w = env * om * np.sin(kk * d) / Ns                 # u_t = N w = om * env sin(k d) (right-mover for k>0)
    ang = 2 * np.pi * s / L; t = 0.0; prev = 0.0; unw = 0.0; target = laps * L
    def cen(u, w):
        if density == "u2":
            e = u ** 2 + w ** 2
        else:                                                                    # A1 diagnostic: conserved energy density (link gradient energy split to sites)
            gl = Nl * (np.roll(u, -1) - u) ** 2; e = w ** 2 + 0.5 * (gl + np.roll(gl, 1))
        z = np.sum(e * np.exp(1j * ang)); return np.angle(z) * L / (2 * np.pi)
    c0 = cen(u, w); prev = c0
    while t < 10 * L:
        k1 = rhs(u, w); k2 = rhs(u + dt / 2 * k1[0], w + dt / 2 * k1[1]); k3 = rhs(u + dt / 2 * k2[0], w + dt / 2 * k2[1]); k4 = rhs(u + dt * k3[0], w + dt * k3[1])
        u = u + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); w = w + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]); t += dt
        c = cen(u, w); dc = (c - prev + L / 2) % L - L / 2; unw_new = unw + dc
        if direction * unw_new >= target:
            return t - dt * (direction * unw_new - target) / (direction * dc)
        unw, prev = unw_new, c
    return None


def loop_dt(Ns, bl, **kw):
    """one-clock Sagnac delay: time(+) - time(-)."""
    return loop_wave(Ns, bl, direction=+1, **kw) - loop_wave(Ns, bl, direction=-1, **kw)


def readout(Ns, bl):
    Nlk = 0.5 * (Ns + np.roll(Ns, -1)); return float(2 * np.sum(bl / Nlk ** 4))


def part_A():
    L = 256; rng = np.random.default_rng(53); res = {}
    one = np.ones(L); b = 1e-3
    pats = {"a_uniform": b * one, "b_gauge": (lambda c: np.roll(c, -1) - c)(1e-3 * rng.standard_normal(L)),
            "c_one_side": np.where(np.arange(L) < 64, b, 0.0)}
    for k, bl in pats.items():
        d = loop_dt(one, bl); res[k] = dict(dt=d, readout=readout(one, bl)); print("53A", k, res[k], flush=True)
    Nw = 1 - 0.01 * np.cos(2 * np.pi * np.arange(L) / L); d = loop_dt(Nw, b * one); res["a_weakN"] = dict(dt=d, readout=readout(Nw, b * one)); print("53A a_weakN", res["a_weakN"], flush=True)
    z = loop_dt(one, 0 * one); res["zero"] = dict(dt=z); print("53A zero", z, flush=True)
    return res


# ------------------------------------------------------------------ 53B
def part_B(n=64, V=0.2):
    c = n // 2; m = S.M; k0 = S.k_for_v(V, m, 0.0); psi, K, w = S.packet_k(n, m, 0.0, (c, c, c), S.SIG, (k0, 0, 0))
    P = link_P(psi); tot = P.sum(axis=(1, 2, 3))
    pk = np.fft.fftn(psi, axes=(0, 1, 2)); wk = np.sum(np.abs(pk) ** 2, -1); wk /= wk.sum(); kk = 2 * np.pi * np.fft.fftfreq(n)
    sink = [float(np.sum(wk * np.sin(kk).reshape([-1 if i == a else 1 for i in range(3)]))) for a in range(3)]
    rng = np.random.default_rng(7); errs = []; Z = np.zeros((3, n, n, n)); h = 1e-4
    for _ in range(20):
        a = int(rng.integers(3)); x = tuple(int(v) for v in rng.integers(c - 6, c + 7, 3)); Bp = Z.copy(); Bp[(a,) + x] = h; Bm = Z.copy(); Bm[(a,) + x] = -h
        dH = (np.real(np.vdot(psi, Hb(psi, Bp))) - np.real(np.vdot(psi, Hb(psi, Bm)))) / (2 * h)
        errs.append(abs(dH + P[(a,) + x]) / max(abs(P).max(), 1e-30))
    # spin current of a spin-polarised packet at rest
    p0, _, _ = S.packet_k(n, m, 0.0, (c, c, c), S.SIG, (0, 0, 0)); P0 = link_P(p0)
    res = dict(total_P=tot.tolist(), sin_k=sink, rel_total=float(abs(tot[0] - sink[0]) / abs(sink[0])), dHdb_err_max=float(max(errs)),
               rest_spin_current=dict(max=float(np.abs(P0).max()), net=P0.sum(axis=(1, 2, 3)).tolist(), moving_max=float(np.abs(P).max())))
    print("53B", res, flush=True); return res


# ------------------------------------------------------------------ static Poisson (Dirichlet faces) for a link component
def poisson(src, bnd, V=0.0):
    """solve (1 - V^2) d_xx + d_yy + d_zz of u = src on the interior, u = bnd on the faces."""
    n = src.shape[0]; mI = n - 2; k = np.arange(1, mI + 1); l1 = 2 * np.cos(np.pi * k / (mI + 1)) - 2
    lam = (1 - V * V) * l1[:, None, None] + l1[None, :, None] + l1[None, None, :]
    u = bnd.copy(); I = (slice(1, -1),) * 3; u[I] = 0; bc = S.lap_aniso(u, V)
    u[I] = idstn(dstn(src[I] - bc, type=1) / lam, type=1); return u


# ------------------------------------------------------------------ 53C ring source
def ring_state(n, m, a_r, s_r, ell, chi, axis=(0, 0, 1)):
    c = n // 2; x = np.arange(n) - c; X = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1).astype(float)
    nz = np.asarray(axis, float); nz /= np.linalg.norm(nz); e1 = np.cross(nz, [1.0, 0, 0] if abs(nz[0]) < 0.9 else [0, 1.0, 0]); e1 /= np.linalg.norm(e1); e2 = np.cross(nz, e1)
    z = X @ nz; u1 = X @ e1; u2 = X @ e2; rc = np.hypot(u1, u2); phi = np.arctan2(u2, u1)
    A = np.exp(-((rc - a_r) ** 2 + z ** 2) / (2 * s_r ** 2)) * np.exp(1j * ell * phi)
    k = 2 * np.pi * np.fft.fftfreq(n); K = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
    E = G.disp(K, 1.0, m, 0.0); del K
    a = np.fft.fftn(A[..., None] * np.asarray(chi, complex), axes=(0, 1, 2)); kk = 2 * np.pi * np.fft.fftfreq(n)
    ha = m * (a @ G.BE.T)
    for j in range(3):
        sk = np.sin(kk).reshape([-1 if i == j else 1 for i in range(3)])[..., None]; ha += sk * (a @ G.AL[j].T)
    a = 0.5 * (a + ha / E[..., None])
    psi = np.fft.ifftn(a, axes=(0, 1, 2))
    for ax in range(3):
        sl = [slice(None)] * 4; sl[ax] = 0; psi[tuple(sl)] = 0; sl[ax] = -1; psi[tuple(sl)] = 0
    return psi / np.linalg.norm(psi)


def ring_sources(n, ell, axis=(0, 0, 1), m=0.6, a_r=10.0, s_r=3.0):
    one = np.ones((n, n, n)); q = 0; P = 0
    for chi in (UP, DOWN):
        psi = ring_state(n, m, a_r, s_r, ell, chi, axis); q = q + 0.5 * G.charge(psi, one, m, 0.0); P = P + 0.5 * link_P(psi)
    return q, P


def link_pos(n):
    c = n // 2; x = np.arange(n) - c; X = np.stack(np.meshgrid(x, x, x, indexing="ij"), 0).astype(float)
    return [X + 0.5 * np.eye(3)[a][:, None, None, None] for a in range(3)]       # positions of link midpoints, each (3, n, n, n)


def dipole_beta(J, Rpos, Gk, lam=LAM):
    """beta_dip = -(lam G / 2) (J x r) / r^3 at positions Rpos (3, ...)."""
    r = np.sqrt(np.sum(Rpos ** 2, 0)); r = np.where(r == 0, 1e-9, r); cr = np.cross(np.asarray(J)[:, None], Rpos.reshape(3, -1), axis=0).reshape(Rpos.shape)
    return -(lam * Gk / 2) * cr / r ** 3


def shift_static(P, kappa, J, lam=LAM):
    n = P.shape[1]; LP = link_pos(n); Gk = kappa / (4 * np.pi); B = np.zeros_like(P)
    for a in range(3):
        bnd = dipole_beta(J, LP[a], Gk, lam)[a]; B[a] = poisson(lam * kappa * P[a], bnd)
    return B


def angmom(P):
    n = P.shape[1]; LP = link_pos(n); J = np.zeros(3)
    for a in range(3):
        v = np.zeros((3,) + P.shape[1:]); v[a] = P[a]; J += np.cross(LP[a].reshape(3, -1), v.reshape(3, -1), axis=0).sum(1)
    return J


# ------------------------------------------------------------------ loops through 3D fields
def square_loop(center, half, plane=(0, 1)):
    """closed lattice square, counterclockwise in the (plane[0], plane[1]) axes; returns list of (site, axis, sign) steps."""
    a, b = plane; p = list(center); p[a] -= half; p[b] -= half; steps = []
    for ax, sg in ((a, 1), (b, 1), (a, -1), (b, -1)):
        for _ in range(2 * half):
            steps.append((tuple(p), ax, sg)); p[ax] += sg
    return steps


def staircase_loop(center, radius, axis):
    """closed lattice path approximating a circle in the plane normal to `axis` (counterclockwise about axis)."""
    nz = np.asarray(axis, float); nz /= np.linalg.norm(nz); e1 = np.cross(nz, [1.0, 0, 0]); e1 /= np.linalg.norm(e1); e2 = np.cross(nz, e1)
    pts = []; c = np.asarray(center, float)
    for t in np.linspace(0, 2 * np.pi, 4000, endpoint=False):
        q = tuple(int(v) for v in np.round(c + radius * (np.cos(t) * e1 + np.sin(t) * e2)))
        if not pts or q != pts[-1]:
            pts.append(q)
    pts.append(pts[0]); steps = []
    for p0, p1 in zip(pts[:-1], pts[1:]):
        p = list(p0)
        for ax in range(3):
            while p[ax] != p1[ax]:
                sg = 1 if p1[ax] > p[ax] else -1; steps.append((tuple(p), ax, sg)); p[ax] += sg
    return steps


def loop_fields(steps, N, B):
    """site N and link shift (signed along the traversal) for a loop of steps."""
    Ns, bl = [], []
    for (x, ax, sg) in steps:
        Ns.append(N[x]); xl = list(x)
        if sg < 0:
            xl[ax] -= 1
        bl.append(sg * B[ax][tuple(xl)])
    return np.array(Ns), np.array(bl)


def loop_dipole(steps, J, Gk, c):
    s = 0.0
    for (x, ax, sg) in steps:
        xl = np.array(x, float) - c; xl[ax] += 0.5 * sg
        s += sg * dipole_beta(J, xl.reshape(3, 1), Gk)[ax, 0]
    return 2 * s


def part_CD(n=129, u=0.01):
    c = n // 2; cv = np.array([c, c, c], float); res = {}
    q0, P0 = ring_sources(n, 0); Q0 = float(q0.sum()); kap = 4 * np.pi * 12 * u / Q0; Gk = kap / (4 * np.pi)
    store = {}
    for ell in (2, 0, -2):
        t0 = time.time(); q, P = (q0, P0) if ell == 0 else ring_sources(n, ell)
        J = angmom(P); N, info = G.static([], [], 0.0, kap, n, (c, c, c), qfix=q); B = shift_static(P, kap, J)
        store[ell] = (q, N, B, J)
        res[str(ell)] = dict(J=J.tolist(), Q=float(q.sum()), info=info, Pmax=float(np.abs(P).max()), wall=time.time() - t0)
        print("53C ell", ell, "J", np.round(J, 6), "Q %.6f" % q.sum(), "%.0fs" % (time.time() - t0), flush=True)
    qp, Np_, Bp, Jp = store[2]; qm, Nm_, Bm, Jm = store[-2]; _, N0_, B0, J0 = store[0]
    res["sym"] = dict(q_diff=float(np.abs(qp - qm).max() / np.abs(qp).max()), N_diff=float(np.abs(Np_ - Nm_).max()))
    out = dict(enclosing={}, inplane={}, axis={})
    for h in (20, 28, 36):
        st = square_loop((c, c, c), h); row = {}
        for ell, (q, N, B, J) in store.items():
            Ns, bl = loop_fields(st, N, B); row[str(ell)] = dict(readout=readout(Ns, bl), wave_energy=loop_dt(Ns, bl, density="energy"),
                                                                 lattice_factor=float(np.cos(0.5) / np.cos(0.25) ** 2))
        row["dipole"] = loop_dipole(st, Jp, Gk, cv); out["enclosing"][str(h)] = row
        print("53D enclosing h", h, {k: (v if k == "dipole" else {kk: float("%.6g" % vv) for kk, vv in v.items()}) for k, v in row.items()}, flush=True)
    for d in (16, 24, 32, 40):
        for key, ctr in (("inplane", (c + d, c, c)), ("axis", (c, c, c + d))):
            st = square_loop(ctr, 4); row = {}
            for ell, (q, N, B, J) in store.items():
                Ns, bl = loop_fields(st, N, B); row[str(ell)] = dict(readout=readout(Ns, bl))
            row["dipole"] = loop_dipole(st, Jp, Gk, cv); out[key][str(d)] = row
        print("53D small d", d, "inplane %.4e (dip %.4e) axis %.4e (dip %.4e)" % (out["inplane"][str(d)]["2"]["readout"], out["inplane"][str(d)]["dipole"],
              out["axis"][str(d)]["2"]["readout"], out["axis"][str(d)]["dipole"]), flush=True)
    res["loops"] = out
    # tilted ring, axis (1,1,1)/sqrt3, ell = +2, one in-plane staircase loop of radius 28
    t0 = time.time(); ax = (1, 1, 1); qt, Pt = ring_sources(n, 2, axis=ax); Jt = angmom(Pt); Nt, it = G.static([], [], 0.0, kap, n, (c, c, c), qfix=qt); Bt = shift_static(Pt, kap, Jt)
    stt = staircase_loop((c, c, c), 28, ax); Ns, bl = loop_fields(stt, Nt, Bt); st_al = staircase_loop((c, c, c), 28, (0, 0, 1)); Na, ba = loop_fields(st_al, Np_, Bp)
    res["tilted"] = dict(J=Jt.tolist(), readout=readout(Ns, bl), dipole=loop_dipole(stt, Jt, Gk, cv), aligned_readout=readout(Na, ba), aligned_dipole=loop_dipole(st_al, Jp, Gk, cv),
                         wall=time.time() - t0)
    res["tilted"]["ratio_tilted"] = res["tilted"]["readout"] / res["tilted"]["dipole"]; res["tilted"]["ratio_aligned"] = res["tilted"]["aligned_readout"] / res["tilted"]["aligned_dipole"]
    print("53D tilted", {k: v for k, v in res["tilted"].items() if k != "J"}, flush=True)
    res["kappa"] = kap; res["G"] = Gk
    return res


# ------------------------------------------------------------------ 53E dynamic co-moving pair with the shift
def evolve_b(psis, ms, f, g, B, gB, kappa, lam, T, dt, record, shift_on=True):
    I = (slice(1, -1),) * 3
    def faces0(a):
        for ax in range(3):
            sl = [slice(None)] * 3; sl[ax] = 0; a[tuple(sl)] = 0; sl[ax] = -1; a[tuple(sl)] = 0
        return a
    def rhs(ps, f, g, B, gB):
        N = Nf(f); dps = [-1j * (G.Hm(p, N, mm, 0.0) + (Hb(p, B) if shift_on else 0)) for p, mm in zip(ps, ms)]
        q = sum(G.charge(p, N, mm, 0.0) for p, mm in zip(ps, ms)); df = faces0(N * g); dg = np.zeros_like(g); dg[I] = N[I] * (F49.op(N, f) - kappa * q[I])
        if shift_on:
            j = sum(link_P(p) for p in ps); dB = np.zeros_like(B); dgB = np.zeros_like(gB)
            for a in range(3):
                Nb = 0.5 * (N + np.roll(N, -1, a)); dB[a] = faces0(Nb * gB[a]); dgB[a][I] = Nb[I] * (F49.op(Nb, B[a]) - lam * kappa * j[a][I])
        else:
            dB = dgB = 0 * B
        return dps, df, dg, dB, dgB
    t = 0.0; hist = [record(psis, f, B, 0.0)]
    def add(st, k, h):
        return [p + h * kp for p, kp in zip(st[0], k[0])], st[1] + h * k[1], st[2] + h * k[2], st[3] + h * k[3], st[4] + h * k[4]
    st = (psis, f, g, B, gB)
    while t < T - 1e-9:
        k1 = rhs(*st); k2 = rhs(*add(st, k1, dt / 2)); k3 = rhs(*add(st, k2, dt / 2)); k4 = rhs(*add(st, k3, dt))
        st = ([p + dt / 6 * (a + 2 * b + 2 * c + d) for p, a, b, c, d in zip(st[0], k1[0], k2[0], k3[0], k4[0])],) + tuple(
            st[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(1, 5))
        t += dt; hist.append(record(st[0], st[1], st[3], t))
    return st, hist


def initial_shift(ps, kap, lam, n, c, V):
    j = sum(link_P(p) for p in ps); Ptot = j.sum(axis=(1, 2, 3)); gam = 1 / np.sqrt(1 - V * V); LP = link_pos(n); B = np.zeros_like(j)
    for a in range(3):
        X = LP[a]; Rg = np.sqrt((gam * X[0]) ** 2 + X[1] ** 2 + X[2] ** 2); Rg = np.where(Rg == 0, 1e-9, Rg)
        B[a] = poisson(lam * kap * j[a], -lam * kap * Ptot[a] * gam / (4 * np.pi * Rg), V)
    gB = np.zeros_like(B)
    for a in range(3):
        gB[a][1:-1] = -V * (B[a][2:] - B[a][:-2]) / 2                            # (1/Nbar) d_t B with N ~ 1 (O(U V) neglected)
    return B, gB


def run_pair(V, chi=UP, n=65, sep=14, u=0.01, T=16.0, dt=0.1, lam=LAM, shift_on=True, kap=None):
    c = n // 2; m = S.M; k0 = S.k_for_v(V, m, 0.0)
    p1, _, _ = S.packet_k(n, m, 0.0, (c, c - sep // 2, c), S.SIG, (k0, 0, 0), chi=chi); p2, _, _ = S.packet_k(n, m, 0.0, (c, c + sep // 2, c), S.SIG, (k0, 0, 0), chi=chi)
    if kap is None:
        p0, _, _ = S.packet_k(n, m, 0.0, (c, c - sep // 2, c), S.SIG, (0, 0, 0)); kap = 4 * np.pi * sep * u / float(G.charge(p0, np.ones((n, n, n)), m, 0.0).sum())
    runs = {}
    for lab, ps in (("pair", [p1, p2]), ("one", [p1])):
        t0 = time.time(); f, g, info = S.initial_field(ps, [m] * len(ps), 0.0, kap, n, c, V)
        if shift_on:
            B, gB = initial_shift(ps, kap, lam, n, c, V)
        else:
            B = gB = np.zeros((3, n, n, n))
        _, hist = evolve_b(ps, [m] * len(ps), f, g, B, gB, kap, lam, T, dt, lambda ps, f, B, t: (t, [S.momentum(p, 1) for p in ps]), shift_on)
        runs[lab] = hist; print("53E V %.1f chi %s %s %.0fs" % (V, chi[:2], lab, time.time() - t0), flush=True)
    t = np.array([h[0] for h in runs["pair"]]); P1 = np.array([h[1][0] for h in runs["pair"]]) - np.array([h[1][0] for h in runs["one"]])
    F1 = float(np.polyfit(t, P1, 2)[1]); print("53E V %.1f chi %s F1 %.6e" % (V, chi[:2], F1), flush=True)
    return dict(V=V, chi=list(chi), kappa=kap, F1=F1, t=t.tolist(), P1=P1.tolist())


def part_E(Vs=(0.0, 0.1, 0.2, 0.3)):
    res = {}
    for V in Vs:
        res[f"{V}_up"] = run_pair(V)
    res["0.2_down"] = run_pair(0.2, chi=DOWN)
    F0 = res["0.0_up"]["F1"]
    for k, v in res.items():
        if v["V"] > 0:
            v["R_dyn"] = v["F1"] / F0 / np.sqrt(1 - v["V"] ** 2); print("53E", k, "R_dyn %.5f" % v["R_dyn"], flush=True)
    return res


def part_C0(n=65, V=0.2):
    """sign/normalization: steady shift of one moving compact source vs -lam kappa P gam/(4 pi R~) at transverse 8..24."""
    c = n // 2; m = S.M; k0 = S.k_for_v(V, m, 0.0); psi, _, _ = S.packet_k(n, m, 0.0, (c, c, c), 2.0, (k0, 0, 0)); kap = 1.0
    B, _ = initial_shift([psi], kap, LAM, n, c, V); Ptot = link_P(psi).sum(axis=(1, 2, 3))[0]; gam = 1 / np.sqrt(1 - V * V)
    rr = np.arange(8, 25); pred = -LAM * kap * Ptot * gam / (4 * np.pi * rr)        # transverse (X = 0, rho = rr): R~ = rho
    meas = B[0][c, c + rr, c]; ratio = meas / pred
    res = dict(ratio_mean=float(ratio.mean()), ratio_min=float(ratio.min()), ratio_max=float(ratio.max()), Ptot=float(Ptot))
    print("53C0", res, flush=True); return res


def save(name, obj, t0):
    OUT.mkdir(parents=True, exist_ok=True)
    rec = dict(study="EX-GC53", part=name, label="EXPLORATORY", platform=platform.platform(), python=platform.python_version(), numpy=np.__version__,
               wall_s=time.time() - t0, result=obj)
    (OUT / f"{name}.json").write_text(json.dumps(rec, indent=1, default=float))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for p in sys.argv[1:]:
        t0 = time.time(); save(p, {"A": part_A, "B": part_B, "C0": part_C0, "CD": part_CD, "E": part_E}[p](), t0); print(p, "done %.0fs" % (time.time() - t0), flush=True)
