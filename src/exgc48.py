"""EX-GC48: Wilson-Dirac lattice fermion under the 44 coupling rules.
48A free inertia R = 2 t_pair E_rest / c^2 (c from a massless packet) + 47 bodies vs the network's max group velocity.
48B packet charges Q_N, Q_a at rest and moving (Wilson readings W1/W2). 48C two Dirac fermions bound by an on-site
contact (no-pair): inertia vs energy, lattice Nordtvedt eta. Plan: spec/EX-GC48_plan.md. EXPLORATORY.
Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "EX-GC48"

# Dirac matrices: Gamma = (I, alpha1, alpha2, alpha3, beta)
s0 = np.eye(2); sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.])
Z = np.zeros((2, 2))
ALPHA = [np.block([[Z, s], [s, Z]]) for s in (sx, sy, sz)]
BETA = np.block([[s0, Z], [Z, -s0]]).astype(complex)
GAM = [np.eye(4, dtype=complex)] + ALPHA + [BETA]
GG = np.array([[np.kron(a, b) for b in GAM] for a in GAM])   # (5,5,16,16)


def parts(k, m, r=1.0, t=1.0, wil="W1"):
    """k (...,3). Returns S-vector components (t sin k_j), mass term M, energy E."""
    sk = t * np.sin(k); W = np.sum(1 - np.cos(k), -1)
    M = m + t * r * W                         # identical for W1/W2 at given t only if t = 1; W2 below
    if wil == "W2":
        M = m + 3 * r - t * r * np.sum(np.cos(k), -1)
    E = np.sqrt(np.sum(sk ** 2, -1) + M ** 2)
    return sk, M, E


def nvec(k, m, r=1.0, t=1.0, wil="W1"):
    sk, M, E = parts(k, m, r, t, wil)
    return np.concatenate([np.ones(E.shape + (1,)), sk / E[..., None], (M / E)[..., None]], -1), E


def Qa1(k, m, r=1.0, wil="W1"):
    """t dE/dt at t = 1 (analytic)."""
    sk, M, E = parts(k, m, r, 1.0, wil)
    S2 = np.sum(sk ** 2, -1)
    if wil == "W1":
        W = np.sum(1 - np.cos(k), -1); return (S2 + r * W * M) / E
    C = np.sum(np.cos(k), -1); return (S2 - r * C * M) / E


def grid(L):
    k = 2 * np.pi * np.fft.fftfreq(L); return np.stack(np.meshgrid(k, k, k, indexing="ij"), -1).reshape(-1, 3)


# ------------------------------------------------------------------ 48A
def part_A(ms=(0.05, 0.1, 0.2, 0.3, 0.5), rs=(1.0, 0.5), h=1e-4):
    rows = {}
    for r in rs:
        for m in ms:
            out = {}
            for name, d in (("axis", np.array([1., 0, 0])), ("diag", np.ones(3) / np.sqrt(3))):
                E = lambda s: parts(s * d, m, r)[2]
                tp = (E(h) + E(-h) - 2 * E(0.0)) / (2 * h * h); out[name] = dict(t_pair=float(tp), R=float(2 * tp * m))
            out["R_pred"] = 1 + r * m; out["rel_axis"] = out["axis"]["R"] / out["R_pred"] - 1; out["rel_diag"] = out["diag"]["R"] / out["R_pred"] - 1
            rows[f"r{r}_m{m}"] = out
            print("48A r", r, "m", m, "R axis %.6f diag %.6f pred %.6f" % (out["axis"]["R"], out["diag"]["R"], out["R_pred"]), flush=True)
    # massless packet speed (1D chain = 3D with transverse k = 0), exact k-space evolution.
    # Deviation D1: the planned k0 = 0.05, sigma = 20 packet straddles k = 0 (massless degeneracy, left movers);
    # use k0 = 0.4, 0.2, 0.1 with k-width k0/4 (sigma = 4/k0) and extrapolate the centroid speed to k0 -> 0.
    L = 4096; k = 2 * np.pi * np.fft.fftfreq(L); K3 = np.stack([k, 0 * k, 0 * k], -1); chi = np.array([1, 0, 0, 0], complex)
    sk, M, E = parts(K3, 0.0, 1.0); Es = np.where(E > 0, E, 1.0)
    H = np.einsum("kj,jab->kab", sk, np.array(ALPHA)) + M[:, None, None] * BETA
    Pi = 0.5 * (np.eye(4) + H / Es[:, None, None]); Pi[E == 0] = 0
    x = np.arange(L); x = np.where(x >= L // 2, x - L, x); hh = 1e-5
    vg = (parts(K3 + [hh, 0, 0], 0.0, 1.0)[2] - parts(K3 - [hh, 0, 0], 0.0, 1.0)[2]) / (2 * hh)
    sig_rows = {}
    for k0 in (0.4, 0.2, 0.1):
        sig = 4 / k0; x0 = -600.0
        dk = (k - k0 + np.pi) % (2 * np.pi) - np.pi; g = np.exp(-dk ** 2 * sig ** 2 / 2) * np.exp(-1j * k * x0)
        a = Pi @ chi; nr = np.linalg.norm(a, axis=1); a = np.where(nr[:, None] > 0, a / np.where(nr > 0, nr, 1)[:, None], 0) * g[:, None]
        ts = np.linspace(0, 400, 21); cen = []
        for tt in ts:
            psi = np.fft.ifft(a * np.exp(-1j * E * tt)[:, None], axis=0); rho = np.sum(np.abs(psi) ** 2, 1); cen.append(float(np.sum(rho * x) / rho.sum()))
        w = np.abs(g) ** 2 * (nr > 0)
        sig_rows[str(k0)] = dict(sigma=sig, c_fit=float(np.polyfit(ts, cen, 1)[0]), v_group_mean=float(np.sum(w * vg) / w.sum()), centroids=cen)
    k0s = np.array([0.4, 0.2, 0.1]); cs = np.array([sig_rows[str(q)]["c_fit"] for q in k0s])
    c_ext = float(np.polyfit(k0s ** 2, cs, 1)[1])
    sig_out = dict(L=L, rows=sig_rows, c_extrapolated=c_ext, c_rel=c_ext - 1)
    print("48A c", {q: round(sig_rows[str(q)]["c_fit"], 6) for q in k0s}, "extrap %.6f" % c_ext, flush=True)
    # 48A-2: 47 bodies vs the 45/47 network's max group velocity (Hellmann-Feynman)
    sys.path.insert(0, str(Path(__file__).parent)); import exgc45 as X
    A47 = json.loads((ROOT / "results" / "EX-GC47" / "A.json").read_text())["result"]; d2 = {}
    for name, m0 in (("B1", 0.5), ("B2", 0.7)):
        e, u, KK = X.eig(64, 0.1, m0=m0, g=0.025); vmax_b = []
        for b in range(7):
            vv = np.zeros(e.shape[:3] + (3,))
            for ax in range(3):
                ub = u[..., :, b]; dH = -1j * np.exp(1j * KK[..., ax])
                vv[..., ax] = 2 * np.real(np.conj(ub[..., 2 * ax]) * dH * ub[..., 2 * ax + 1])
            vmax_b.append(float(np.sqrt(np.sum(vv ** 2, -1)).max()))
        vmax = max(vmax_b); tp = A47[name]["t_pair"]; Eb = A47[name]["E"]
        d2[name] = dict(vmax_by_band=vmax_b, vmax=vmax, t_pair=tp, E=Eb, R47=2 * tp * Eb / vmax ** 2)
        print("48A-2", name, "vmax %.5f R47 %.3e" % (vmax, d2[name]["R47"]), flush=True)
    return dict(free=rows, signal=sig_out, bodies47=d2)


# ------------------------------------------------------------------ 48B
def packet(m, sig, k0, r=1.0, n=81, wil="W1"):
    """Continuous-k quadrature on a local box around k0 (infinite lattice)."""
    hw = min(6.0 / sig, np.pi); s = np.linspace(-hw, hw, n); dk = np.stack(np.meshgrid(s, s, s, indexing="ij"), -1).reshape(-1, 3)
    k = dk + np.asarray(k0); w = np.exp(-np.sum(dk ** 2, -1) * sig ** 2); w /= w.sum()
    E = parts(k, m, r)[2]; Q = Qa1(k, m, r, wil)
    h = 1e-5; vx = (parts(k + [h, 0, 0], m, r)[2] - parts(k - [h, 0, 0], m, r)[2]) / (2 * h)
    vy = (parts(k + [0, h, 0], m, r)[2] - parts(k - [0, h, 0], m, r)[2]) / (2 * h); vz = (parts(k + [0, 0, h], m, r)[2] - parts(k - [0, 0, h], m, r)[2]) / (2 * h)
    v2 = vx ** 2 + vy ** 2 + vz ** 2
    # lapse charge: homogeneous weights => N dE/dN = E per mode
    EN = float(np.sum(w * E)); return dict(E=EN, QN_over_E=float(np.sum(w * E)) / EN, Qa=float(np.sum(w * Q)), Qa_over_E=float(np.sum(w * Q)) / EN, v2=float(np.sum(w * v2)))


def part_B(ms=(0.2, 0.5), r=1.0):
    res = {}
    for m in ms:
        rest = {}
        for sig in (4, 8, 16, 32):
            p1 = packet(m, sig, (0, 0, 0)); p2 = packet(m, sig, (0, 0, 0), wil="W2")
            rest[str(sig)] = dict(W1=p1, W2=p2, pred_W1=3 * (1 + r * m / 2) / (2 * sig ** 2 * m ** 2))
        slope = float(np.log(rest["32"]["W1"]["Qa_over_E"] / rest["16"]["W1"]["Qa_over_E"]) / np.log(2))
        mov = {}; p0 = packet(m, 24, (0, 0, 0)); fac = (1 + r * m / 2) / (1 + r * m) ** 2
        for k0 in (0.05, 0.1, 0.2):
            p = packet(m, 24, (k0, 0, 0)); dQ = p["Qa_over_E"] - p0["Qa_over_E"]; dv = p["v2"] - p0["v2"]
            mov[str(k0)] = dict(pk=p, dQ_over_E=dQ, dv2=dv, pred=dv * fac, rel=dQ / (dv * fac) - 1, ratio_to_v2=dQ / dv)
        res[str(m)] = dict(rest=rest, slope_16_32=slope, rel_32=rest["32"]["W1"]["Qa_over_E"] / rest["32"]["pred_W1"] - 1, moving=mov, rest24=p0,
                           W2_rest_over_pred=rest["32"]["W2"]["Qa_over_E"] / (-3 * r / m))
        print("48B m", m, "slope %.4f rel32 %.4f" % (slope, res[str(m)]["rel_32"]), {k: round(v["rel"], 4) for k, v in mov.items()},
              "W2 Qa/E %.4f" % rest["32"]["W2"]["Qa_over_E"], flush=True)
    return res


# ------------------------------------------------------------------ 48C
class Pair:
    def __init__(self, L, m, r=1.0, t=1.0, K=(0, 0, 0)):
        k = grid(L); self.N = L ** 3; Kv = np.asarray(K, float)
        self.n1, E1 = nvec(k + Kv, m, r, t); self.n2, E2 = nvec(-k, m, r, t); self.D = E1 + E2; self.Dmin = float(self.D.min())
        self.k = k; self.m = m; self.r = r; self.K = Kv

    def Mmat(self, E, U=1.0):
        S = (self.n1 * (1 / (self.D - E))[:, None]).T @ self.n2
        M = np.einsum("mn,mnab->ab", S, GG) * U / (4 * self.N); return 0.5 * (M + M.conj().T)

    def lam(self, E, U=1.0):
        return np.linalg.eigvalsh(self.Mmat(E, U))[-1]

    def solve(self, U):
        f = lambda E: U * self.lam(E) - 1
        lo = self.Dmin - 5
        E0 = brentq(f, lo, self.Dmin - 1e-12, xtol=1e-14, rtol=1e-14); return E0

    def wavefn(self, E0, U):
        w, V = np.linalg.eigh(self.Mmat(E0, U)); phi = V[:, -1]
        g = np.real(np.einsum("a,mnab,b->mn", phi.conj(), GG, phi))
        wk = np.einsum("km,kn,mn->k", self.n1, self.n2, g) / 4 / (self.D - E0) ** 2; return wk / wk.sum(), float(w[-1] - w[-2])


def U_for(L, m, Eb, r=1.0):
    P = Pair(L, m, r); return 1.0 / P.lam(P.Dmin - Eb)


def ddE_axis(k, m, r, ax, h=1e-4):
    e = np.zeros(3); e[ax] = h
    return (parts(k + e, m, r)[2] + parts(k - e, m, r)[2] - 2 * parts(k, m, r)[2]) / h ** 2


def comp_point(L, m, frac, r=1.0, U=None, h=1e-5):
    Eb_t = frac * 2 * m
    if U is None:
        U = U_for(L, m, Eb_t, r)
    P0 = Pair(L, m, r); E0 = P0.solve(U); wk, gap = P0.wavefn(E0, U); Dexp = float(np.sum(wk * P0.D))
    T = Dexp - 2 * m; V = E0 - Dexp; Eb = P0.Dmin - E0
    K1 = 2 * np.pi / L; out = dict(L=L, m=m, frac=frac, U=U, E0=E0, Eb=Eb, T=T, V=V, gap_top2=gap, Dmin=P0.Dmin)
    for name, d in (("axis", np.array([1., 0, 0])), ("diag", np.ones(3))):
        dE = [Pair(L, m, r, K=n * K1 * d).solve(U) - E0 for n in (1, 2)]; K2 = float(np.sum((K1 * d) ** 2))
        out[name] = dict(dE=dE, t_c=(16 * dE[0] - dE[1]) / (12 * K2), t_c1=dE[0] / K2)
    tspec = float(np.sum(wk * (ddE_axis(P0.k, m, r, 0) + ddE_axis(-P0.k, m, r, 0))) / 8)
    Rf = 1 + r * m
    for name in ("axis", "diag"):
        o = out[name]; o["R_c"] = 2 * o["t_c"] * E0; o["dR"] = o["R_c"] / Rf - 1
    out["dR_pred"] = V / Dexp; out["t_spec"] = tspec; out["spec_rel"] = tspec / out["axis"]["t_c"] - 1
    Ep = Pair(L, m, r, t=1 + h).solve(U); Em = Pair(L, m, r, t=1 - h).solve(U); Qa = (Ep - Em) / (2 * h)
    out.update(Qa=Qa, Qa_over_E=Qa / E0, eta=Qa / Eb, QN_over_E=1.0)
    return out


def part_C(pts=((0.5, 96), (0.2, 128)), fracs=(0.02, 0.08), r=1.0):
    res = {}
    for m, L in pts:
        for fr in fracs:
            o = comp_point(L, m, fr, r); Ls = int(round(0.75 * L / 2) * 2); s = comp_point(Ls, m, fr, r, U=o["U"])
            o["size_check"] = dict(L=Ls, Eb=s["Eb"], t_c=s["axis"]["t_c"], dR=s["axis"]["dR"], eta=s["eta"])
            res[f"{m}_{fr}"] = o
            print("48C m", m, "fr", fr, "L", L, "Eb %.5f T %.4f V %.4f | dR axis %.4f diag %.4f pred %.4f | spec %.4f | eta %.3f | L%d: dR %.4f eta %.3f" % (
                o["Eb"], o["T"], o["V"], o["axis"]["dR"], o["diag"]["dR"], o["dR_pred"], o["spec_rel"], o["eta"], Ls, s["axis"]["dR"], s["eta"]), flush=True)
    return res


def save(name, obj, t0):
    OUT.mkdir(parents=True, exist_ok=True)
    rec = dict(study="EX-GC48", part=name, label="EXPLORATORY", platform=platform.platform(), python=platform.python_version(), numpy=np.__version__,
               wall_s=time.time() - t0, result=obj)
    (OUT / f"{name}.json").write_text(json.dumps(rec, indent=1, default=float))


if __name__ == "__main__":
    which = sys.argv[1:] or ["A", "B", "C"]
    for p in which:
        t0 = time.time(); save(p, {"A": part_A, "B": part_B, "C": part_C}[p](), t0); print(p, "done %.1fs" % (time.time() - t0), flush=True)
