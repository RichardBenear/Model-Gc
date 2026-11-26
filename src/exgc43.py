"""EX-GC43: what the lapse couples to. Under one functional the source is q_v = <dH/dN_v> (37B weights: links N_vN_w,
sites N_v), so link energy counts twice. Vacuum kernel chi(q); hub charges; positive-energy packets of the unchanged
network. Plan: spec/EX-GC43_plan.md. Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc40a as A
import exgc40c as C
import exgc42 as E

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC43"
M = 6
J = E.J


# ---------------------------------------------------------------- uniform vacuum
def vac_derivs(lam0, L=32, h=1e-3):
    f = {s: E.rho_vac(L, lam0, 1 + s * h) for s in (-1, 0, 1)}
    eps, d1, d2 = f[0], (f[1] - f[-1]) / (2 * h), (f[1] - 2 * f[0] + f[-1]) / h ** 2
    return dict(eps=eps, sigma=d1, chi0=d2, eps_link=d1 - eps, eps_site=2 * eps - d1)   # eps = l + s, d1 = 2l + s


# ---------------------------------------------------------------- supercell conjugate density
def supercell_q(lam0, nfield, S, Lfull):
    """q_v = sum_links N_w E_link,vw + E_site,v/N_v... i.e. <dH/dN_v>, per supercell site, lapse N = 1 + nfield."""
    Sx, Sy, Sz = S; ns = Sx * Sy * Sz
    N = 1.0 + nfield.reshape(-1)
    coords = np.array([(x, y, z) for x in range(Sx) for y in range(Sy) for z in range(Sz)])
    sid = {tuple(c): i for i, c in enumerate(coords)}
    grids = [2 * np.pi * np.arange(Lfull[i] // S[i]) / (Lfull[i] // S[i]) for i in range(3)]
    base = np.zeros((ns * M, ns * M), complex)
    links = []
    for i, c in enumerate(coords):
        base[i * M:(i + 1) * M, i * M:(i + 1) * M] += N[i] * lam0 * J
        for a in range(3):
            d = list(c); d[a] += 1; wrap = d[a] // S[a]; d[a] %= S[a]
            links.append((i, sid[tuple(d)], i * M + 2 * a, sid[tuple(d)] * M + 2 * a + 1, a, wrap))
    q = np.zeros(ns); cnt = 0
    for ta in grids[0]:
        for tb in grids[1]:
            for tc in grids[2]:
                th = (ta, tb, tc); h = base.copy(); ph = []
                for (i, j, p, r, a, wrap) in links:
                    f = np.exp(1j * th[a] * wrap); ph.append(f)
                    h[p, r] += -N[i] * N[j] * f; h[r, p] += -N[i] * N[j] * np.conj(f)
                e, v = np.linalg.eigh(h); occ = v[:, :ns * M // 2]
                P = occ @ occ.conj().T                       # P_ij = <c_j^dag c_i>... energy = tr(h P^T) convention below
                Pd = P.reshape(ns, M, ns, M)
                for i in range(ns):
                    q[i] += lam0 * np.real(np.sum(J * Pd[i, :, i, :]))        # site term / N_v
                for (i, j, p, r, a, wrap), f in zip(links, ph):
                    El = 2 * np.real(-f * P[r, p])                             # unweighted link energy h0_pr P_rp + c.c.
                    q[i] += N[j] * El; q[j] += N[i] * El
                cnt += 1
    return q / cnt


def chi_axis(lam0, L=32, ms=range(1, 9), eps=1e-4):
    out = {}; x = np.arange(L)
    for m in ms:
        qq = 2 * np.pi * m / L; res = []
        for s in (+1, -1):
            res.append(supercell_q(lam0, (s * eps * np.cos(qq * x)).reshape(L, 1, 1), (L, 1, 1), (L, L, L)))
        dq = (res[0] - res[1]) / 2
        chi = float((1 if 2 * m == L else 2) / L * np.sum(dq * np.cos(qq * x)) / eps)
        out[m] = dict(q=qq, qt2=2 * (1 - np.cos(qq)), chi=chi); print("chi", lam0, m, chi, flush=True)
    return out


# ---------------------------------------------------------------- packets (Bloch)
def band_data(L, lam0):
    H, K = A.bloch_h(L, lam0)
    Hs = lam0 * J
    Hl = H - Hs
    e, u = np.linalg.eigh(H)
    el = np.real(np.einsum("...ib,...ij,...jb->...b", u.conj(), Hl, u))
    es = np.real(np.einsum("...ib,...ij,...jb->...b", u.conj(), np.broadcast_to(Hs, H.shape), u))
    return e, el, es, u, K


def packet_weights(L, k0, sigma):
    k = 2 * np.pi * np.fft.fftfreq(L)
    K = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
    dk = (K - np.array(k0) + np.pi) % (2 * np.pi) - np.pi
    w = np.exp(-sigma ** 2 * np.sum(dk ** 2, -1)); return w / w.sum()      # |c_k|^2 with c_k ~ exp(-s^2 dk^2 / 2)


def packet_QE(bd, b, k0, sigma, L):
    e, el, es, u, K = bd
    w = packet_weights(L, k0, sigma)
    dE = float(np.sum(w * e[..., b])); dEl = float(np.sum(w * el[..., b])); dEs = float(np.sum(w * es[..., b]))
    return dict(dE=dE, dE_link=dEl, dE_site=dEs, Q=2 * dEl + dEs)


def kmin(bd, b, sign=+1):
    e = bd[0][..., b]; i = np.unravel_index(np.argmin(sign * e), e.shape)
    return [float(bd[4][i][a]) for a in range(3)], float(e[i])


# ---------------------------------------------------------------- real-space checks
def apply_link(psi):
    """H_link psi on (L,L,L,6): (H psi)[x,2a] = -psi[x+e_a,2a+1], (H psi)[x,2a+1] = -psi[x-e_a,2a]."""
    out = np.zeros_like(psi)
    for a in range(3):
        out[..., 2 * a] = -np.roll(psi[..., 2 * a + 1], -1, axis=a)
        out[..., 2 * a + 1] = -np.roll(psi[..., 2 * a], +1, axis=a)
    return out


def check_c3(lam0=0.1, L=16, sigma=2.0, ntimes=200, tmax=50.0):
    bd = band_data(L, lam0); b = M // 2; k0, _ = kmin(bd, b)
    ref = packet_QE(bd, b, k0, sigma, L)
    e, el, es, u, K = bd
    w = packet_weights(L, k0, sigma); ck = np.sqrt(w)
    x = np.arange(L); X = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1)
    x0 = np.array([L // 2] * 3)
    def realspace(phase_t):
        amp = ck * np.exp(-1j * np.einsum("...a,a->...", K, x0)) * phase_t       # centre the packet at x0
        spec = amp[..., None] * u[..., :, b]                                       # (L,L,L,6) in k
        return np.fft.ifftn(spec, axes=(0, 1, 2)) * L ** 1.5                       # unitary: sum |psi|^2 = 1
    def QE(psi):
        nrm = np.sum(np.abs(psi) ** 2)
        hl = np.real(np.vdot(psi, apply_link(psi))) / nrm
        hs = np.real(np.vdot(psi.reshape(-1, M), (psi.reshape(-1, M) @ (lam0 * J).T).reshape(-1))) / nrm
        return dict(dE=hl + hs, Q=2 * hl + hs, norm=float(nrm))
    t0 = QE(realspace(1.0))
    Qs = [QE(realspace(np.exp(-1j * e[..., b] * t)))["Q"] for t in np.linspace(0, tmax, ntimes)]
    return dict(band_formula=ref, realspace_t0=t0, diff_t0=abs(t0["Q"] - ref["Q"]), time_mean_Q=float(np.mean(Qs)),
                time_rel=abs(np.mean(Qs) / ref["Q"] - 1))


def check_c1(L=8, lam0=0.1, dl=0.2, eps=1e-4):
    out = {}
    one = np.ones((L, L, L))
    for name, hub in (("vacuum", None), ("hub", dl)):
        h0 = E.dense_H(L, lam0, one, hub)
        ev, V = np.linalg.eigh(h0); occ = V[:, :h0.shape[0] // 2]; P = occ @ occ.conj().T
        rows = []
        for v in ((0, 0, 0), (1, 0, 0)):
            Np = one.copy(); Np[v] += eps; Nm = one.copy(); Nm[v] -= eps
            hp, hm = E.dense_H(L, lam0, Np, hub), E.dense_H(L, lam0, Nm, hub)
            dH = (hp - hm) / (2 * eps)
            q_hf = float(np.real(np.trace(dH @ P.T)))
            Ep = np.linalg.eigvalsh(hp)[:h0.shape[0] // 2].sum(); Em = np.linalg.eigvalsh(hm)[:h0.shape[0] // 2].sum()
            q_fd = float((Ep - Em) / (2 * eps))
            rows.append(dict(site=v, q_HF=q_hf, q_FD=q_fd, rel=abs(q_hf / q_fd - 1)))
        out[name] = rows; print("C1", name, rows, flush=True)
    return out


def hub_charges(lam0, L=16):
    res = {}
    for key, hs in (("1", [(0, 0, 0)]), ("7", C.NEIGH)):
        for dl in (0.1, 0.2, 0.4):
            r = E.hub_source(L, lam0, dl, hs)
            Q = float(np.sum(r["S"] + r["Slink"]))
            # uniform-N derivative of dE0 (hub on-site weighted N, as in 42A); 7-site hub via Woodbury on the cluster
            dE = {}
            for s in (-1, 1):
                Nv = 1 + s * 1e-3
                dE[s] = hub_dE_uniform(L, lam0, Nv, dl, hs)
            dQ = (dE[1] - dE[-1]) / 2e-3
            res[f"{key}_{dl}"] = dict(dE0=r["dE0"], sumS=r["sumS"], Q=Q, dE0_dN=float(dQ), Q_over_dE=Q / r["dE0"],
                                      rel=abs(Q / dQ - 1))
            print("hub", lam0, key, dl, res[f"{key}_{dl}"], flush=True)
    return res


def hub_dE_uniform(L, lam0, N, dl, hubsites, nq=160):
    H = E.bloch_H(L, lam0, N)
    e = np.linalg.eigvalsh(H); lo, hi = e[..., M // 2 - 1].max(), e[..., M // 2].min(); mu = 0.5 * (lo + hi); gap = hi - lo
    t, w = np.polynomial.legendre.leggauss(nq); t = 0.5 * (t + 1); w = 0.5 * w
    ys = gap * t / (1 - t); wy = w * gap / (1 - t) ** 2
    hub = np.array(hubsites); D = np.kron(np.eye(len(hub)), N * dl * J); dE = 0.0
    for y, wt in zip(ys, wy):
        G = C.block(C.Gr_at(H, mu + 1j * y), hub, hub, L)
        dE += -wt / np.pi * np.log(abs(np.linalg.det(np.eye(len(D)) - G @ D)))
    return float(dE)


def far_point(chi_fun, chi0, fourpiG, Lg, variant):
    n, g = E.far_field([(0, 0, 0)], [1.0], chi_fun, fourpiG, Lg=Lg, variant=variant, K0v=chi0)
    return n, g


if __name__ == "__main__":
    part = sys.argv[1]; OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    if part == "vacuum":
        res = {}
        for lam0 in (0.1, 0.2):
            res[str(lam0)] = dict(uniform=vac_derivs(lam0), axis=chi_axis(lam0))
    elif part == "checks":
        res = dict(C1=check_c1(), C3=check_c3())
    elif part == "hubs":
        res = {str(l): hub_charges(l) for l in (0.1, 0.2)}
    elif part == "packets":
        res = {}
        for lam0 in (0.1, 0.2):
            L = 32; bd = band_data(L, lam0); e, el, es = bd[0], bd[1], bd[2]
            mu = 0.5 * (e[..., M // 2 - 1].max() + e[..., M // 2].min())
            rows = {}
            bc, bt = M // 2, M - 1
            kc, ec = kmin(bd, bc); kt, et = kmin(bd, bt)
            for name, b, k0 in (("cond_k0", bc, [0, 0, 0]), ("cond_min", bc, kc), ("top_k0", bt, [0, 0, 0]), ("top_min", bt, kt)):
                for s in (1.0, 2.0, 4.0):
                    r = packet_QE(bd, b, k0, s, L); r["ratio"] = r["Q"] / r["dE"]; rows[f"{name}_s{s:g}"] = r
            kv, _ = kmin(bd, M // 2 - 1, sign=-1)
            pc = packet_QE(bd, bc, kc, 2.0, L); pv = packet_QE(bd, M // 2 - 1, kv, 2.0, L)
            rows["pair_s2"] = dict(dE=pc["dE"] - pv["dE"], Q=pc["Q"] - pv["Q"], ratio=(pc["Q"] - pv["Q"]) / (pc["dE"] - pv["dE"]))
            bmap = {}
            for b in range(M):
                ratio = 1 + el[..., b] / e[..., b]; ok = np.abs(e[..., b]) > 1e-3
                bmap[b] = dict(e_range=[float(e[..., b].min()), float(e[..., b].max())],
                               ratio_range=[float(ratio[ok].min()), float(ratio[ok].max())],
                               near1=int(np.sum(np.abs(ratio[ok] - 1) < 1e-3)), near2=int(np.sum(np.abs(ratio[ok] - 2) < 1e-3)))
            res[str(lam0)] = dict(mu=float(mu), k_cond_min=kc, e_cond_min=ec, packets=rows, band_map=bmap)
            for k, v in rows.items():
                print(lam0, k, {a: round(b, 5) for a, b in v.items()}, flush=True)
    elif part == "far":
        vac = json.loads((OUT / "vacuum.json").read_text())["result"]; res = {}
        for lam0 in (0.1, 0.2):
            d = vac[str(lam0)]; chi0 = d["uniform"]["chi0"]
            qx = np.array([0.0] + [v["qt2"] for v in d["axis"].values()]); cx = np.array([chi0] + [v["chi"] for v in d["axis"].values()])
            o = np.argsort(qx); qx, cx = qx[o], cx[o]
            fun = lambda q2, qx=qx, cx=cx: np.interp(q2, qx, cx)
            c2 = float(np.polyfit(qx[1:5], cx[1:5] - chi0, 1)[0])
            ext = float(np.polyfit(qx[1:3], cx[1:3], 1)[1])
            r = dict(chi0=chi0, c2chi=c2, chi0_extrap=ext, C2_rel=abs(ext / chi0 - 1),
                     eps_link=d["uniform"]["eps_link"], P1_chi0_le_2link=chi0 <= 2 * d["uniform"]["eps_link"], cases={})
            for fg in (0.01, 0.05):
                G = fg / (4 * np.pi)
                nL, g = far_point(fun, chi0, fg, 128, "A"); fL = E.fit_A(nL, g)
                nV, _ = far_point(fun, chi0, fg, 128, "B"); fV = E.fit_A(nV, g)
                sweep = {}
                for Lg in (96, 112, 128, 144, 160):
                    nn, _ = far_point(fun, chi0, fg, Lg, "B"); sweep[Lg] = float(nn[10, 0, 0])
                kchi = float(np.sqrt(fg * abs(chi0) / (1 + fg * c2)))
                dir_ = {}
                for R in (30, 40, 50, 60):
                    rr = 10.0
                    dir_[R] = float(-(fg / (1 + fg * c2)) / (4 * np.pi) * (np.cos(kchi * rr) - np.sin(kchi * rr) / np.tan(kchi * R)) / rr)
                r["cases"][str(fg)] = dict(L=dict(fL, ratio=fL["GeffM"] / G, pred=1 / (1 + fg * c2)), V_resid_1overr=fV["max_rel_resid"],
                                           V_over_L_resid=fV["max_rel_resid"] / fL["max_rel_resid"], sweep_n10=sweep, kchi=kchi,
                                           dirichlet_n10=dir_, kchi_R_over_pi={R: kchi * R / np.pi for R in (30, 40, 50, 60)})
                print(lam0, fg, r["cases"][str(fg)], flush=True)
            res[str(lam0)] = r
    res = dict(result=res, wall_seconds=time.time() - t0, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("done", part, res["wall_seconds"])
