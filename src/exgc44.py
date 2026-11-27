"""EX-GC44: separating the clock from space on the 40A pair network. 44A end-to-end equal-energy benchmark (43's 37B
model); 44B lapse-homogeneous weights (links Nbar, sites N); 44C two fields (lapse N, spatial link factor a; hop ~ Nbar/a);
44D band-projected exciton and its excess strain stress. Plan: spec/EX-GC44_plan.md. Our construction; not LaViolette's
Model G, not Chain Theory.
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc40a as A
import exgc40c as C
import exgc42 as E
import exgc43 as F

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC44"
M = 6
J = E.J


def bloch(L, lam0, c_link=1.0, c_site=1.0):
    """H = c_link * H_link + c_site * lam0 J (Bloch), plus K grid."""
    H, K = A.bloch_h(L, lam0)
    Hs = lam0 * J
    return c_link * (H - Hs) + c_site * Hs, K


def bands_g(L, lam0, c_link=1.0, c_site=1.0, chi=None):
    """Bands with a smooth gauge: each band vector's phase fixed by its overlap with a fixed vector chi; the degenerate
    flat pairs (3,4) and (0,1) are replaced by the normalised projection of chi onto their subspace (slot 3 / slot 0)."""
    H, K = bloch(L, lam0, c_link, c_site)
    e, u = np.linalg.eigh(H)
    if chi is None:
        chi = np.random.default_rng(20261127).normal(size=M) + 1j * np.random.default_rng(7).normal(size=M)
    for b in range(M):
        ov = np.einsum("...i,i->...", u[..., :, b].conj(), chi)
        u[..., :, b] *= np.exp(1j * np.angle(ov))[..., None]          # makes u^dag chi real and positive
    for (b0, b1) in ((0, 1), (3, 4)):
        P = np.einsum("...ib,...jb->...ij", u[..., :, b0:b1 + 1], u[..., :, b0:b1 + 1].conj())
        v = P @ chi; v /= np.linalg.norm(v, axis=-1, keepdims=True); u[..., :, b0] = v
    return e, u, K, H


# ---------------------------------------------------------------- 44A: real-space conjugate-charge profiles (37B model)
def psi_real(u_b, w, K, x0, L):
    ck = np.sqrt(w) * np.exp(-1j * np.einsum("...a,a->...", K, np.array(x0, float)))
    return np.fft.ifftn(ck[..., None] * u_b, axes=(0, 1, 2)) * L ** 1.5


def q_profile(psi, lam0):
    """37B conjugate density at N = 1: q_v = E_site,v + sum of the full energies of the 6 links touching v."""
    site = np.real(np.einsum("...i,ij,...j->...", psi.conj(), lam0 * J, psi))
    q = site.copy(); elink_tot = 0.0
    for a in range(3):
        El = 2 * np.real(psi[..., 2 * a].conj() * (-1) * np.roll(psi[..., 2 * a + 1], -1, axis=a))   # link v -> v + e_a
        q += El + np.roll(El, 1, axis=a); elink_tot += El.sum()
    return q, float(site.sum() + 2 * elink_tot), float(elink_tot)       # Q = E_site + 2 E_link


def band_QE(e, u, H, lam0, b, w):
    Hl = H - lam0 * J
    el = np.real(np.einsum("...i,...ij,...j->...", u[..., :, b].conj(), Hl, u[..., :, b]))
    ee = np.real(np.einsum("...i,...ij,...j->...", u[..., :, b].conj(), H, u[..., :, b]))
    return float(np.sum(w * ee)), float(np.sum(w * (ee + el)))


def gauss_monopole(n, R):
    """G_eff*Q from the exact lattice Gauss law: sum of (Lap n) over the cube |x|,|y|,|z| <= R equals the face flux;
    Lap n = 4 pi G_eff (Q delta-like source - Q/V) beyond the source, so divide by 4 pi (1 - V_in/V). Added after the
    registered Green's-function fit was seen to be biased by the anisotropy of the neutral sources (post hoc)."""
    lap = sum(np.roll(n, 1, a) + np.roll(n, -1, a) for a in range(3)) - 6 * n
    idx = np.r_[0:R + 1, -R:0]; Vin = (2 * R + 1) ** 3
    return float(lap[np.ix_(idx, idx, idx)].sum() / (4 * np.pi * (1 - Vin / n.size)))


def part_44A(L=32, Lg=128, fourpiG=0.01):
    """End-to-end: localized sources built by projecting a fixed orbital chi (Gaussian envelope sigma) onto the gapped
    conduction group (bands 3-5) or valence group (0-2). Single-band packets (the registered recipes) are not localized:
    band 5 touches the flat bands at k = 0 and the band gauges have obstructions (deviation, recorded in the note)."""
    vac = json.loads((Path(__file__).resolve().parents[1] / "results" / "EX-GC43" / "vacuum.json").read_text())["result"]
    chis = {"e0": np.eye(M)[0], "uni": np.ones(M) / np.sqrt(M)}
    res = {}
    for lam0 in (0.1, 0.2):
        H, K = bloch(L, lam0); e, u = np.linalg.eigh(H); Hl = H - lam0 * J
        Pc = np.einsum("...ib,...jb->...ij", u[..., :, 3:], u[..., :, 3:].conj())
        Pv = np.einsum("...ib,...jb->...ij", u[..., :, :3], u[..., :, :3].conj())
        x0 = [L // 2] * 3
        def st(P, chi, sigma):
            w = F.packet_weights(L, [0, 0, 0], sigma); v = P @ chis[chi]
            nrm = np.real(np.einsum("...i,...i->...", v.conj(), v)); Z = np.sum(w * nrm)
            ee = np.real(np.einsum("...i,...ij,...j->...", v.conj(), H, v)); el = np.real(np.einsum("...i,...ij,...j->...", v.conj(), Hl, v))
            psi = psi_real(v / np.sqrt(Z), w, K, x0, L)
            q, Qr, _ = q_profile(psi, lam0)
            return dict(q=q, Q=Qr, Qb=float(np.sum(w * (ee + el)) / Z), dE=float(np.sum(w * ee) / Z))
        def comb(a, b):   # particle a plus hole b
            return dict(q=a["q"] - b["q"], Q=a["Q"] - b["Q"], Qb=a["Qb"] - b["Qb"], dE=a["dE"] - b["dE"])
        srcs = {"neutral_e0_x_uni": comb(st(Pc, "e0", 2.0), st(Pv, "uni", 2.0)),
                "neutral_e0_x_e0": comb(st(Pc, "e0", 2.0), st(Pv, "e0", 2.0)),
                "charged_uni_s2": st(Pc, "uni", 2.0)}
        target = srcs["neutral_e0_x_uni"]["dE"]
        f = lambda sg: st(Pc, "uni", sg)["dE"] - target
        sig = None
        try:
            if f(0.5) * f(20.0) < 0:
                sig = brentq(f, 0.5, 20.0, xtol=1e-12); srcs["charged_uni_equalE"] = st(Pc, "uni", sig)
        except Exception:
            sig = None
        d = vac[str(lam0)]; chi0 = d["uniform"]["chi0"]
        qx = np.array([0.0] + [v["qt2"] for v in d["axis"].values()]); cx = np.array([chi0] + [v["chi"] for v in d["axis"].values()])
        o = np.argsort(qx); fun = lambda q2, qx=qx[o], cx=cx[o]: np.interp(q2, qx, cx)
        c2 = float(np.polyfit(qx[o][1:5], cx[o][1:5] - chi0, 1)[0]); Gpred = fourpiG / (4 * np.pi) / (1 + fourpiG * c2)
        x = np.arange(L) - L // 2; R = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2)
        out = dict(sigma_charged_equalE=sig, G_eff_pred=Gpred, sources={})
        sites = [(a - x0[0], b - x0[1], c - x0[2]) for a in range(L) for b in range(L) for c in range(L)]
        for name, s_ in srcs.items():
            n, g = E.far_field(sites, s_["q"].reshape(-1), fun, fourpiG, Lg=Lg, variant="A", K0v=chi0)
            fa = E.fit_A(n, g, rs=range(6, 21))
            out["sources"][name] = dict(dE=s_["dE"], Q=s_["Q"], Q_band=s_["Qb"], C1=abs(s_["Q"] - s_["Qb"]), Q_over_dE=s_["Q"] / s_["dE"],
                                        frac_Q_within_r6=float(s_["q"][R <= 6].sum() / s_["Q"]),
                                        GeffM=fa["GeffM"], amp_over_Q_over_Gpred=fa["GeffM"] / s_["Q"] / Gpred, resid=fa["max_rel_resid"],
                                        gauss={R: gauss_monopole(n, R) / s_["Q"] / Gpred for R in (10, 16, 24)})
            print("44A", lam0, name, out["sources"][name], flush=True)
        res[str(lam0)] = out
    return res


# ---------------------------------------------------------------- 44B/44C: uniform derivatives of band and hub energies
def band_weights_QE(L, lam0, cl_fun, cs_fun, h=1e-5):
    """d/dx at x = 1 of packet energies for H(x) = cl(x) H_link + cs(x) H_site; returns per-state derivative arrays."""
    ep = np.linalg.eigvalsh(bloch(L, lam0, cl_fun(1 + h), cs_fun(1 + h))[0])
    em = np.linalg.eigvalsh(bloch(L, lam0, cl_fun(1 - h), cs_fun(1 - h))[0])
    return (ep - em) / (2 * h)


def hub_E(L, lam0, dl, hubsites, c_link, c_site, nq=160):
    H = bloch(L, lam0, c_link, c_site)[0]
    e = np.linalg.eigvalsh(H); lo, hi = e[..., M // 2 - 1].max(), e[..., M // 2].min(); mu = 0.5 * (lo + hi); gap = hi - lo
    t, w = np.polynomial.legendre.leggauss(nq); t = 0.5 * (t + 1); w = 0.5 * w
    ys = gap * t / (1 - t); wy = w * gap / (1 - t) ** 2
    hub = np.array(hubsites); D = np.kron(np.eye(len(hub)), c_site * dl * J); dE = 0.0
    for y, wt in zip(ys, wy):
        G = C.block(C.Gr_at(H, mu + 1j * y), hub, hub, L)
        dE += -wt / np.pi * np.log(abs(np.linalg.det(np.eye(len(D)) - G @ D)))
    return float(dE)


def part_44BC(L=32, h=1e-5):
    p43 = json.loads((Path(__file__).resolve().parents[1] / "results" / "EX-GC43" / "packets.json").read_text())["result"]
    h43 = json.loads((Path(__file__).resolve().parents[1] / "results" / "EX-GC43" / "hubs.json").read_text())["result"]
    res = {}
    for lam0 in (0.1, 0.2):
        e0, K = np.linalg.eigvalsh(bloch(L, lam0)[0]), bloch(L, lam0)[1]
        dN_hom = band_weights_QE(L, lam0, lambda x: x, lambda x: x, h)          # 44B: links Nbar -> N, sites N
        dLna = -band_weights_QE(L, lam0, lambda x: 1 / x, lambda x: 1.0, h)     # 44C: -a dH/da, hop ~ 1/a
        dN_37B = band_weights_QE(L, lam0, lambda x: x * x, lambda x: x, h)       # C2 mapping: a = 1/Nbar -> N^2 links
        rows = {}
        kc = p43[str(lam0)]["k_cond_min"]
        defs = [("cond_k0", 3, [0, 0, 0]), ("top_k0", 5, [0, 0, 0])]
        for name, b, k0 in defs:
            for s in (1.0, 2.0, 4.0):
                w = F.packet_weights(L, k0, s)
                dE = float(np.sum(w * e0[..., b])); QN = float(np.sum(w * dN_hom[..., b])); Qa = float(np.sum(w * dLna[..., b]))
                Q37 = float(np.sum(w * dN_37B[..., b])); key = f"{name}_s{s:g}"
                rows[key] = dict(dE=dE, QN=QN, Qa=Qa, QN_over_dE=QN / dE, Qa_over_dE=Qa / dE, bend_vs_dust=1 + Qa / (2 * QN),
                                 C2_37B=abs(Q37 - p43[str(lam0)]["packets"][key]["Q"]))
        # pair (flat cond x valence s2, as 43's pair_s2)
        wc = F.packet_weights(L, kc, 2.0); kv = [0, 0, 0]; wv = F.packet_weights(L, kv, 2.0)
        def pair(arr): return float(np.sum(wc * arr[..., 3]) - np.sum(wv * arr[..., 2]))
        dE = pair(e0); QN = pair(dN_hom); Qa = pair(dLna); Q37 = pair(dN_37B)
        rows["pair_s2"] = dict(dE=dE, QN=QN, Qa=Qa, QN_over_dE=QN / dE, Qa_over_dE=Qa / dE, bend_vs_dust=1 + Qa / (2 * QN),
                               C2_37B=abs(Q37 - p43[str(lam0)]["packets"]["pair_s2"]["Q"]))
        # optics exponent under 44B: group velocity of the top band at k ~ (2pi/L)(1,0,0) vs N
        vel = {}
        for Nv in (0.9, 1.0, 1.1):
            ee = np.linalg.eigvalsh(bloch(L, lam0, Nv, Nv)[0])[..., 5]
            vel[Nv] = float((ee[2, 0, 0] - ee[1, 0, 0]) / (2 * np.pi / L))
        p_opt = float(np.log(vel[1.1] / vel[0.9]) / np.log(1.1 / 0.9))
        hubs = {}
        for key, hs in (("1", [(0, 0, 0)]), ("7", C.NEIGH)):
            for dl in (0.1, 0.2, 0.4):
                e1 = hub_E(16, lam0, dl, hs, 1.0, 1.0)
                eN = {Nv: hub_E(16, lam0, dl, hs, Nv, Nv) for Nv in (0.9, 1.1)}
                ep_ = hub_E(16, lam0, dl, hs, 1 + 1e-3, 1 + 1e-3); em_ = hub_E(16, lam0, dl, hs, 1 - 1e-3, 1 - 1e-3)
                ap = hub_E(16, lam0, dl, hs, 1 / (1 + 1e-3), 1.0); am = hub_E(16, lam0, dl, hs, 1 / (1 - 1e-3), 1.0)
                QN = (ep_ - em_) / 2e-3; Qa = -(ap - am) / 2e-3
                hubs[f"{key}_{dl}"] = dict(dE0=e1, QN=QN, Qa=Qa, QN_over_dE=QN / e1, Qa_over_dE=Qa / e1,
                                           EP_change=float((eN[1.1] / 1.1) / (eN[0.9] / 0.9) - 1),
                                           C2_37B=abs(QN + Qa - h43[str(lam0)][f"{key}_{dl}"]["Q"]))
                print("44BC hub", lam0, key, dl, hubs[f"{key}_{dl}"], flush=True)
        res[str(lam0)] = dict(packets=rows, optics=dict(vel=vel, p=p_opt), hubs=hubs)
        for k, v in rows.items():
            print("44BC", lam0, k, v, flush=True)
    return res


def part_C0(Lg=128, fourpiG=0.01):
    zero = lambda q2: 0.0 * q2
    out = {}
    for name in ("Psi_unitQN", "Phi_unitQN_plus_Qa"):
        n, g = E.far_field([(0, 0, 0)], [1.0], zero, fourpiG, Lg=Lg, variant="B")
        fa = E.fit_A(n, g); out[name] = dict(GeffM=fa["GeffM"], expected=fourpiG / (4 * np.pi), rel=abs(fa["GeffM"] / (fourpiG / (4 * np.pi)) - 1))
    return out


# ---------------------------------------------------------------- 44D: band-projected exciton (separable contact)
def exciton_setup(L, lam0, c_link=1.0):
    H, K = bloch(L, lam0, c_link, 1.0)
    e, u = np.linalg.eigh(H)
    return e, u


def exciton_E0(e, u, U, L, Kshift=(0, 0, 0), bp=5, bv=2, Emin_guess=None):
    """Lowest bound energy of D - (U/N) Phi Phi^dag, D(k) = e_p(k+K) - e_v(k), Phi_ij(k) = conj(u_p(k+K)_i) u_v(k)_j."""
    ep = np.roll(e[..., bp], [-s for s in Kshift], axis=(0, 1, 2)).reshape(-1)
    up = np.roll(u[..., :, bp], [-s for s in Kshift], axis=(0, 1, 2)).reshape(-1, M)
    ev = e[..., bv].reshape(-1); uv = u[..., :, bv].reshape(-1, M)
    D = ep - ev; N = L ** 3
    Phi = (up.conj()[:, :, None] * uv[:, None, :]).reshape(N, M * M)
    Dmin = D.min()
    def lam(E):
        G = Phi.conj().T @ (Phi / (D - E)[:, None]); return np.linalg.eigvalsh(0.5 * (G + G.conj().T))[-1]
    f = lambda E: U * lam(E) / N - 1.0
    lo = Dmin - 10.0; hi = Dmin - 1e-9
    if f(hi) < 0:
        return None
    E0 = brentq(f, lo, hi, xtol=1e-13)
    G = Phi.conj().T @ (Phi / (D - E0)[:, None]); w, V = np.linalg.eigh(0.5 * (G + G.conj().T))
    psi = (Phi @ V[:, -1]) / (D - E0); psi /= np.linalg.norm(psi)
    Dexp = float(np.real(np.sum(np.abs(psi) ** 2 * D)))
    return dict(E0=float(E0), Dmin=float(Dmin), Eb=float(Dmin - E0), D_exp=Dexp, V_exp=float(E0 - Dexp),
                size2=float(np.sum(np.abs(psi) ** 2 * 0)))


def U_for_Eb(e, u, L, Eb, bp=5, bv=2):
    ep = e[..., bp].reshape(-1); ev = e[..., bv].reshape(-1); up = u[..., :, bp].reshape(-1, M); uv = u[..., :, bv].reshape(-1, M)
    D = ep - ev; N = L ** 3; E = D.min() - Eb
    Phi = (up.conj()[:, :, None] * uv[:, None, :]).reshape(N, M * M)
    G = Phi.conj().T @ (Phi / (D - E)[:, None])
    return float(N / np.linalg.eigvalsh(0.5 * (G + G.conj().T))[-1])


def part_44D(L=32, fracs=(0.05, 0.15, 0.30), ws=(0.0, 0.5, 1.0, 2.0), h=1e-5):
    res = {}
    for lam0 in (0.1, 0.2):
        e, u = exciton_setup(L, lam0); ep_, up_ = exciton_setup(L, lam0, 1 + h); em_, um_ = exciton_setup(L, lam0, 1 - h)
        gap = float((e[..., 5] - 0).min() - e[..., 2].max())
        edge = lambda ee: float(ee[..., 5].min() - ee[..., 2].max())
        edge_link = (edge(ep_) - edge(em_)) / (2 * h)            # t d(edge)/dt: link energy of the unbound pair at band edges
        rows = {}
        for fr in fracs:
            Eb = fr * gap; U = U_for_Eb(e, u, L, Eb)
            r0 = exciton_E0(e, u, U, L)
            rp = exciton_E0(ep_, up_, U, L); rm = exciton_E0(em_, um_, U, L)
            Elink = (rp["E0"] - rm["E0"]) / (2 * h)
            rU = exciton_E0(e, u, U * (1 + 1e-5), L); rUm = exciton_E0(e, u, U * (1 - 1e-5), L)
            V_HF = U * (rU["E0"] - rUm["E0"]) / (2e-5 * U)
            dlink = Elink - edge_link
            wstar = -dlink / r0["V_exp"]
            disp = {}
            for Ks in ((1, 0, 0), (2, 0, 0)):
                rk = exciton_E0(e, u, U, L, Kshift=Ks); disp[str(Ks)] = (rk["E0"] - r0["E0"]) if rk else None
            Lc = 24; e24, u24 = exciton_setup(Lc, lam0); r24 = exciton_E0(e24, u24, U, Lc)
            rows[str(fr)] = dict(U=U, Eb=r0["Eb"], E0=r0["E0"], V_exp=r0["V_exp"], V_HF=V_HF, C4=abs(V_HF / r0["V_exp"] - 1),
                                 Elink_bound=Elink, Elink_edges=edge_link, d_link=dlink,
                                 dQa={str(w): dlink + w * r0["V_exp"] for w in ws}, w_star=wstar,
                                 dispersion=disp, Eb_L24=(r24["Eb"] if r24 else None))
            print("44D", lam0, fr, rows[str(fr)], flush=True)
        # C3: U = 0 spectrum (continuum minimum = band-edge sum) and an independent dense solve at L = 10
        Ld = 10; ed, ud = exciton_setup(Ld, lam0); U = rows[str(fracs[1])]["U"]
        ep = ed[..., 5].reshape(-1); ev = ed[..., 2].reshape(-1); upp = ud[..., :, 5].reshape(-1, M); uvv = ud[..., :, 2].reshape(-1, M)
        Amat = upp.conj() @ upp.T                         # A(k',k) = sum_i conj(u_p(k')_i) u_p(k)_i
        Bmat = uvv @ uvv.conj().T                         # B(k',k) = sum_j u_v(k')_j conj(u_v(k)_j)
        Hd = np.diag(ep - ev).astype(complex) - (U / Ld ** 3) * Amat * Bmat
        Edense = float(np.linalg.eigvalsh(Hd)[0]); Esep = exciton_E0(ed, ud, U, Ld)
        res[str(lam0)] = dict(gap=gap, edge_link=edge_link, rows=rows,
                              C3=dict(U0_cont_min=float((ep - ev).min()), edge=edge(ed), E_dense_L10=Edense,
                                      E_separable_L10=(Esep["E0"] if Esep else None),
                                      diff=(abs(Edense - Esep["E0"]) if Esep else None)))
        print("44D C3", lam0, res[str(lam0)]["C3"], flush=True)
    return res


if __name__ == "__main__":
    part = sys.argv[1]; OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    res = {"A": part_44A, "BC": part_44BC, "C0": part_C0, "D": part_44D}[part]()
    res = dict(result=res, wall_seconds=time.time() - t0, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)) + "\n")
    print("done", part, res["wall_seconds"])
