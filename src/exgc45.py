"""EX-GC45: heavy on-site matter on the 40A pair network. 45A: 44's exciton total spatial charge at its true dispersion
minimum. 45B: one heavy carrier (on-site mode m0, hybridised by g to the site's link modes). 45C: bound pair of heavy
carriers of two species; lattice Nordtvedt ratio eta = dQ_a/E_b. Lapse-homogeneous weights (44B), hop ~ 1/a.
Plan: spec/EX-GC45_plan.md. Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json, sys, time, itertools
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC45"
M6 = 6


def bloch(L, lam0, t=(1.0, 1.0, 1.0), m0=None, g=0.0, Nsc=1.0):
    """Bloch H: 6 link modes (hop -t_a e^{ik_a}), on-site lam0(J-I); optional heavy mode 6 at m0, coupled g to modes 0..5.
    Every term times Nsc (uniform lapse, homogeneous weights)."""
    M = 6 if m0 is None else 7
    k = 2 * np.pi * np.fft.fftfreq(L); K = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
    H = np.zeros((L, L, L, M, M), complex)
    for a in range(3):
        H[..., 2 * a, 2 * a + 1] = -t[a] * np.exp(1j * K[..., a]); H[..., 2 * a + 1, 2 * a] = -t[a] * np.exp(-1j * K[..., a])
    H[..., :6, :6] += lam0 * (np.ones((6, 6)) - np.eye(6))
    if m0 is not None:
        H[..., 6, 6] = m0; H[..., 6, :6] = g; H[..., :6, 6] = g
    return Nsc * H, K


def eig(L, lam0, **kw):
    H, K = bloch(L, lam0, **kw); e, u = np.linalg.eigh(H); return e, u, K


# ---------------------------------------------------------------- 45A exciton (44 model, bands 5 and 2)
def exc_solve(e, u, U, Ks, bp=5, bv=2):
    L = e.shape[0]; N = L ** 3
    ep = np.roll(e[..., bp], [-s for s in Ks], axis=(0, 1, 2)).reshape(-1); up = np.roll(u[..., :, bp], [-s for s in Ks], axis=(0, 1, 2)).reshape(-1, 6)
    ev = e[..., bv].reshape(-1); uv = u[..., :, bv].reshape(-1, 6)
    D = ep - ev; Phi = (up.conj()[:, :, None] * uv[:, None, :]).reshape(N, 36); Dmin = D.min()
    def lam(E):
        G = Phi.conj().T @ (Phi / (D - E)[:, None]); return np.linalg.eigvalsh(0.5 * (G + G.conj().T))[-1]
    f = lambda E: U * lam(E) / N - 1
    if f(Dmin - 1e-10) < 0:
        return None
    E0 = brentq(f, Dmin - 10, Dmin - 1e-10, xtol=1e-13)
    G = Phi.conj().T @ (Phi / (D - E0)[:, None]); w, V = np.linalg.eigh(0.5 * (G + G.conj().T))
    psi = (Phi @ V[:, -1]) / (D - E0); psi /= np.linalg.norm(psi)
    Dexp = float(np.sum(np.abs(psi) ** 2 * D))
    return dict(E0=float(E0), Dmin=float(Dmin), Eb=float(Dmin - E0), V=float(E0 - Dexp))


def part_A(fracs=(0.05, 0.15, 0.30), Lscan=24, Lev=48, h=1e-5):
    D44 = json.loads((Path(__file__).resolve().parents[1] / "results" / "EX-GC44" / "D.json").read_text())["result"]
    res = {}
    for lam0 in (0.1, 0.2):
        es, us, _ = eig(Lscan, lam0)
        el, ul, _ = eig(Lev, lam0)
        rows = {}
        for fr in fracs:
            U = D44[str(lam0)]["rows"][str(fr)]["U"]
            best = None
            for nx in range(Lscan // 2 + 1):
                for ny in range(nx, Lscan // 2 + 1):
                    for nz in range(ny, Lscan // 2 + 1):
                        r = exc_solve(es, us, U, (nx, ny, nz))
                        if r and (best is None or r["E0"] < best[1]["E0"]):
                            best = ((nx, ny, nz), r)
            Kn = best[0]; Ks = tuple(int(round(n * Lev / Lscan)) for n in Kn)
            r0 = exc_solve(el, ul, U, Ks)
            def E0_t(t):
                e_, u_, _ = eig(Lev, lam0, t=t); rr = exc_solve(e_, u_, U, Ks); return rr["E0"], rr["Dmin"]
            iso = [E0_t((1 + s * h,) * 3) for s in (1, -1)]
            Lb = (iso[0][0] - iso[1][0]) / (2 * h); Le = (iso[0][1] - iso[1][1]) / (2 * h)
            dirs = {}
            for a in range(3):
                tp = [1.0] * 3; tm = [1.0] * 3; tp[a] += h; tm[a] -= h
                dirs["xyz"[a]] = (E0_t(tuple(tp))[0] - E0_t(tuple(tm))[0]) / (2 * h)
            rows[str(fr)] = dict(U=U, K_star_units_2pi_over_24=Kn, K_star_L48=Ks, E0_Kstar=r0["E0"],
                                 E0_K0=D44[str(lam0)]["rows"][str(fr)]["E0"], Eb_Kstar=r0["Eb"], V=r0["V"],
                                 Qa_total_w0=Lb, Qa_binding_only=Lb - Le, Qa_over_dE=Lb / r0["E0"], Qa_dir=dirs,
                                 w_total_cancel=-Lb / r0["V"])
            print("45A", lam0, fr, rows[str(fr)], flush=True)
        res[str(lam0)] = rows
    return res


# ---------------------------------------------------------------- 45B heavy carrier
def heavy_band(L, lam0, m0, g, t=(1.0, 1.0, 1.0)):
    e, u, K = eig(L, lam0, t=t, m0=m0, g=g); return e[..., 3], u[..., :, 3], K


def part_B(L=32, lam0=0.1, m0s=(0.3, 0.5, 0.7), gs=(0.025, 0.05, 0.1, 0.2), h=1e-5):
    res = {}
    for m0 in m0s:
        for g in gs:
            e, u, K = heavy_band(L, lam0, m0, g)
            i = np.unravel_index(np.argmin(e), e.shape); emin = float(e[i]); W = float(e.max() - e.min())
            adm = float(1 - abs(u[i][6]) ** 2)
            ep, _, _ = heavy_band(L, lam0, m0, g, t=(1 + h,) * 3); em, _, _ = heavy_band(L, lam0, m0, g, t=(1 - h,) * 3)
            dl = (ep - em) / (2 * h)
            # curvature (t_eff) along x at the minimum: e(k) ~ emin + t_eff k^2
            ix = list(i); ix[0] = (i[0] + 1) % L; dk = 2 * np.pi / L
            ixm = list(i); ixm[0] = (i[0] - 1) % L
            teff = float((e[tuple(ix)] + e[tuple(ixm)] - 2 * emin) / (2 * dk ** 2))
            km = (L // 8, 0, 0)
            res[f"{m0}_{g}"] = dict(m0=m0, g=g, kmin=[int(x) for x in i], e_min=emin, bandwidth=W, admixture=adm, t_eff=teff,
                                    Qa_min=float(dl[i]), Qa_over_dE_min=float(dl[i] / emin), R_Qa_over_W=float(dl[i] / W),
                                    e_move=float(e[km]), kin_move=float(e[km] - emin), Qa_over_dE_move=float(dl[km] / e[km]))
            print("45B", m0, g, res[f"{m0}_{g}"], flush=True)
    # C2: g = 0
    e0, _, _ = heavy_band(L, lam0, 0.5, 0.0); ep, _, _ = heavy_band(L, lam0, 0.5, 0.0, t=(1 + h,) * 3)
    res["C2_g0"] = dict(flat_spread=float(e0.max() - e0.min()), Qa_max=float(np.abs((ep - e0) / h).max()))
    # C1: homogeneity Q_N = dE for the carrier
    eN = {s: np.linalg.eigvalsh(bloch(L, lam0, m0=0.5, g=0.1, Nsc=1 + s * h)[0])[..., 3].min() for s in (1, -1)}
    e1 = np.linalg.eigvalsh(bloch(L, lam0, m0=0.5, g=0.1)[0])[..., 3].min()
    res["C1_carrier"] = dict(QN=float((eN[1] - eN[-1]) / (2 * h)), dE=float(e1))
    return res


# ---------------------------------------------------------------- 45C bound heavy pair (two species, rank-1 contact)
def pair_solve(e, u, U, Ks=(0, 0, 0)):
    L = e.shape[0]; N = L ** 3
    e1 = np.roll(e, [-s for s in Ks], axis=(0, 1, 2)); u1 = np.roll(u[..., 6], [-s for s in Ks], axis=(0, 1, 2))
    neg = lambda a: np.roll(np.flip(a, axis=(0, 1, 2)), 1, axis=(0, 1, 2))            # f(-k)
    e2 = neg(e); u2 = neg(u[..., 6])
    D = (e1 + e2).reshape(-1); phi = (np.conj(u1) * np.conj(u2)).reshape(-1); Dmin = D.min()
    S = lambda E: float(np.sum(np.abs(phi) ** 2 / (D - E)))
    f = lambda E: U * S(E) / N - 1
    if f(Dmin - 1e-14) < 0:
        return None
    E0 = brentq(f, Dmin - 50, Dmin - 1e-14, xtol=1e-15, rtol=1e-14)
    psi = phi / (D - E0); psi /= np.linalg.norm(psi)
    Dexp = float(np.sum(np.abs(psi) ** 2 * D))
    return dict(E0=float(E0), Dmin=float(Dmin), Eb=float(Dmin - E0), V=float(E0 - Dexp), T=float(Dexp - Dmin))


def U_for(e, u, Eb):
    L = e.shape[0]; N = L ** 3
    neg = lambda a: np.roll(np.flip(a, axis=(0, 1, 2)), 1, axis=(0, 1, 2))
    D = (e + neg(e)).reshape(-1); phi = (np.conj(u[..., 6]) * np.conj(neg(u[..., 6]))).reshape(-1)
    E = D.min() - Eb
    return float(N / np.sum(np.abs(phi) ** 2 / (D - E)))


def part_C(L=32, lam0=0.1, pts=((0.3, 0.2), (0.5, 0.1), (0.7, 0.05)), fracs=(0.05, 0.15, 0.30), h=1e-5, Lscan=16):
    res = {}
    for (m0, g) in pts:
        e, u, K = eig(L, lam0, m0=m0, g=g); eh, uh = e[..., 3], u[..., :, 3]
        W = float(eh.max() - eh.min())
        ep_, up_, _ = eig(L, lam0, m0=m0, g=g, t=(1 + h,) * 3); em_, um_, _ = eig(L, lam0, m0=m0, g=g, t=(1 - h,) * 3)
        es, us, _ = eig(Lscan, lam0, m0=m0, g=g)
        for fr in fracs:
            Eb = fr * W; U = U_for(eh, uh, Eb)
            r0 = pair_solve(eh, uh, U)
            rp = pair_solve(ep_[..., 3], up_[..., :, 3], U); rm = pair_solve(em_[..., 3], um_[..., :, 3], U)
            Qa_pair = (rp["E0"] - rm["E0"]) / (2 * h); Qa_edges = (rp["Dmin"] - rm["Dmin"]) / (2 * h)
            dQa = Qa_pair - Qa_edges
            rU = pair_solve(eh, uh, U * (1 + 1e-6)); rUm = pair_solve(eh, uh, U * (1 - 1e-6))
            VHF = (rU["E0"] - rUm["E0"]) / 2e-6
            # dispersion minimum over K (wedge, L = 16)
            best = None
            for nx in range(Lscan // 2 + 1):
                for ny in range(nx, Lscan // 2 + 1):
                    for nz in range(ny, Lscan // 2 + 1):
                        rr = pair_solve(es[..., 3], us[..., :, 3], U, (nx, ny, nz))
                        if rr and (best is None or rr["E0"] < best[1]):
                            best = ((nx, ny, nz), rr["E0"])
            # homogeneity C1
            eNp = np.linalg.eigh(bloch(L, lam0, m0=m0, g=g, Nsc=1 + h)[0]); eNm = np.linalg.eigh(bloch(L, lam0, m0=m0, g=g, Nsc=1 - h)[0])
            QN = (pair_solve(eNp[0][..., 3], eNp[1][..., :, 3], U * (1 + h))["E0"] - pair_solve(eNm[0][..., 3], eNm[1][..., :, 3], U * (1 - h))["E0"]) / (2 * h)
            row = dict(m0=m0, g=g, W=W, U=U, E0=r0["E0"], Eb=r0["Eb"], V=r0["V"], T=r0["T"], V_HF=VHF, C4=abs(VHF / r0["V"] - 1),
                       Qa_pair=Qa_pair, Qa_edges=Qa_edges, dQa=dQa, eta=dQa / r0["Eb"], dQa_over_T=dQa / r0["T"],
                       Qa_pair_over_dE=Qa_pair / r0["E0"], K_min_L16=best[0], E0_Kmin_L16=best[1], QN=QN, C1=abs(QN / r0["E0"] - 1))
            res[f"{m0}_{g}_{fr}"] = row
            print("45C", m0, g, fr, row, flush=True)
    return res


# ---------------------------------------------------------------- checks C3 (dense L = 10) and C5 (4 upper bands, L = 4/6)
def pair_dense_full(L, lam0, m0, g, U, bands=(3,), t=(1.0, 1.0, 1.0)):
    e, u, K = eig(L, lam0, m0=m0, g=g, t=t); N = L ** 3
    neg = lambda a: np.roll(np.flip(a, axis=(0, 1, 2)), 1, axis=(0, 1, 2))
    E1 = [e[..., b].reshape(-1) for b in bands]; E2 = [neg(e[..., b]).reshape(-1) for b in bands]
    H1 = [u[..., 6, b].reshape(-1) for b in bands]; H2 = [neg(u[..., 6, b]).reshape(-1) for b in bands]
    nb = len(bands); dim = N * nb * nb
    D = np.zeros(dim); A = np.zeros(dim, complex); idx = 0
    for i in range(nb):
        for j in range(nb):
            D[idx * N:(idx + 1) * N] = E1[i] + E2[j]; A[idx * N:(idx + 1) * N] = np.conj(H1[i]) * np.conj(H2[j]); idx += 1
    # V(k',k) = -(U/N) sum_v n_h1 n_h2: pair (k, -k) -> amplitude conj(u1(k')_6) u1(k)_6 conj(u2(-k')_6) u2(-k)_6
    Hm = np.diag(D).astype(complex) - (U / N) * np.outer(A, A.conj())
    return np.linalg.eigvalsh(Hm)[0]


def part_checks(lam0=0.1, m0=0.5, g=0.1, h=1e-5):
    C = json.loads((OUT / "C.json").read_text())["result"][f"{m0}_{g}_0.15"]
    U = C["U"]; out = {}
    Ed = pair_dense_full(10, lam0, m0, g, U)
    e, u, _ = eig(10, lam0, m0=m0, g=g); Es = pair_solve(e[..., 3], u[..., :, 3], U)["E0"]
    out["C3"] = dict(E_dense=float(Ed), E_sep=float(Es), diff=abs(Ed - Es))
    for Lc in (4, 6):
        proj = [pair_dense_full(Lc, lam0, m0, g, U, (3,), t=(1 + s * h,) * 3) for s in (1, -1)]
        full = [pair_dense_full(Lc, lam0, m0, g, U, (3, 4, 5, 6), t=(1 + s * h,) * 3) for s in (1, -1)]
        e_, u_, _ = eig(Lc, lam0, m0=m0, g=g)
        ep_, up_, _ = eig(Lc, lam0, m0=m0, g=g, t=(1 + h,) * 3); em_, um_, _ = eig(Lc, lam0, m0=m0, g=g, t=(1 - h,) * 3)
        neg = lambda a: np.roll(np.flip(a, axis=(0, 1, 2)), 1, axis=(0, 1, 2))
        edge = lambda ee: float((ee[..., 3] + neg(ee[..., 3])).min())
        Qe = (edge(ep_) - edge(em_)) / (2 * h)
        dq_p = (proj[0] - proj[1]) / (2 * h) - Qe; dq_f = (full[0] - full[1]) / (2 * h) - Qe
        out[f"C5_L{Lc}"] = dict(dQa_proj=float(dq_p), dQa_4band=float(dq_f), rel=float(abs(dq_f - dq_p) / abs(dq_p)))
        print("C5", Lc, out[f"C5_L{Lc}"], flush=True)
    print("C3", out["C3"], flush=True)
    return out


if __name__ == "__main__":
    part = sys.argv[1]; OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    res = {"A": part_A, "B": part_B, "C": part_C, "checks": part_checks}[part]()
    res = dict(result=res, wall_seconds=time.time() - t0, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)) + "\n")
    print("done", part, res["wall_seconds"])
