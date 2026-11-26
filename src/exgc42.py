"""EX-GC42: the 40A pair network coupled to a lapse N by 37B event weights (links N_v N_w, on-site terms N_v), with N
fixed by lattice Poisson sourced by the network's energy density. Variant B: raw source; variant A: local vacuum energy
rho_vac(N_v) subtracted. Plan: spec/EX-GC42_plan.md. Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc40a as A
import exgc40c as C

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC42"
M = 6
J = np.ones((M, M)) - np.eye(M)


# ---------------------------------------------------------------- uniform N (Bloch)
def bloch_H(L, lam0, N):
    k = 2 * np.pi * np.fft.fftfreq(L)
    K = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
    H = np.zeros((L, L, L, M, M), complex)
    for a in range(3):
        ph = np.exp(1j * K[..., a])
        H[..., 2 * a, 2 * a + 1] = -N * N * ph
        H[..., 2 * a + 1, 2 * a] = -N * N * np.conj(ph)
    H += N * lam0 * J
    return H


def rho_vac(L, lam0, N):
    e = np.linalg.eigvalsh(bloch_H(L, lam0, N))
    return float(e[..., :M // 2].sum() / L ** 3)


def hub_uniform(L, lam0, N, dl, nq=160):
    """dE0 and in-gap states of a one-site hub (on-site N*dl*J) in the uniform-N network."""
    H = bloch_H(L, lam0, N)
    e = np.linalg.eigvalsh(H)
    lo, hi = e[..., M // 2 - 1].max(), e[..., M // 2].min(); mu = 0.5 * (lo + hi); gap = hi - lo
    t, w = np.polynomial.legendre.leggauss(nq); t = 0.5 * (t + 1); w = 0.5 * w
    ys = gap * t / (1 - t); wy = w * gap / (1 - t) ** 2
    D = N * dl * J; hub = np.array([(0, 0, 0)]); dE = 0.0
    for y, wt in zip(ys, wy):
        G = C.block(C.Gr_at(H, mu + 1j * y), hub, hub, L)
        dE += -wt / np.pi * np.log(abs(np.linalg.det(np.eye(M) - G @ D)))
    # bound states
    Es = np.linspace(lo + 1e-3 * gap, hi - 1e-3 * gap, 300); Dinv = np.linalg.inv(D); prev = None; cross = []
    for E in Es:
        G = C.block(C.Gr_at(H, E + 0j), hub, hub, L); ev = np.linalg.eigvalsh(0.5 * (Dinv - G + (Dinv - G).conj().T))
        if prev is not None:
            cross += [float(E) for k in range(M) if np.sign(ev[k]) != np.sign(prev[k])]
        prev = ev
    return dict(dE0=float(dE), gap=float(gap), mu=float(mu), bound=cross)


# ---------------------------------------------------------------- supercell response kernel
def supercell_rho(lam0, nfield, S, Lfull):
    """Energy density per supercell site for lapse N = 1 + nfield (shape S), Bloch phases over the full lattice Lfull.
    Supercell S = (Sx,Sy,Sz) sites; theta grids of size Lfull[i]/S[i]."""
    Sx, Sy, Sz = S
    ns = Sx * Sy * Sz
    N = 1.0 + nfield.reshape(-1)
    coords = np.array([(x, y, z) for x in range(Sx) for y in range(Sy) for z in range(Sz)])
    sid = {tuple(c): i for i, c in enumerate(coords)}
    grids = [2 * np.pi * np.arange(Lfull[i] // S[i]) / (Lfull[i] // S[i]) for i in range(3)]
    base = np.zeros((ns * M, ns * M), complex)
    rho = np.zeros(ns); cnt = 0
    links = []   # (i_mode, j_mode, weight, axis, wraps)
    for i, c in enumerate(coords):
        base[i * M:(i + 1) * M, i * M:(i + 1) * M] += N[i] * lam0 * J
        for a in range(3):
            d = list(c); d[a] += 1; wrap = d[a] // S[a]; d[a] %= S[a]
            j = sid[tuple(d)]
            links.append((i * M + 2 * a, j * M + 2 * a + 1, -N[i] * N[j], a, wrap))
    for ta in grids[0]:
        for tb in grids[1]:
            for tc in grids[2]:
                th = (ta, tb, tc)
                h = base.copy()
                for (p, q, wgt, a, wrap) in links:
                    ph = np.exp(1j * th[a] * wrap)
                    h[p, q] += wgt * ph; h[q, p] += wgt * np.conj(ph)
                e, v = np.linalg.eigh(h)
                occ = v[:, :ns * M // 2]; eo = e[:ns * M // 2]
                dens = (np.abs(occ) ** 2) @ eo                    # (hP)_ii = sum_occ e |v_i|^2
                rho += dens.reshape(ns, M).sum(1); cnt += 1
    return rho / cnt


def kernel_axis(lam0, L=32, ms=range(1, 9), eps=1e-4):
    out = {}
    for m in ms:
        q = 2 * np.pi * m / L
        x = np.arange(L)
        res = []
        for s in (+1, -1):
            n = (s * eps * np.cos(q * x)).reshape(L, 1, 1)
            res.append(supercell_rho(lam0, n, (L, 1, 1), (L, L, L)))
        dr = (res[0] - res[1]) / 2
        K = float((1 if 2 * m == L else 2) / L * np.sum(dr * np.cos(q * x)) / eps)   # q = pi: cos^2 = 1, norm 1/L
        out[m] = dict(q=q, qt2=2 * (1 - np.cos(q)), K_B=K)
        print("K", lam0, m, K, flush=True)
    return out


def kernel_diag(lam0, L=16, ms=(1, 2), eps=1e-4):
    out = {}
    for m in ms:
        q = 2 * np.pi * m / L
        X, Y = np.meshgrid(np.arange(L), np.arange(L), indexing="ij")
        res = []
        for s in (+1, -1):
            n = (s * eps * np.cos(q * (X + Y))).reshape(L, L, 1)
            res.append(supercell_rho(lam0, n, (L, L, 1), (L, L, L)).reshape(L, L))
        dr = (res[0] - res[1]) / 2
        K = float(2 / L ** 2 * np.sum(dr * np.cos(q * (X + Y))) / eps)
        out[m] = dict(q=q, qt2=4 * (1 - np.cos(q)), K_B=K)
        print("Kdiag", lam0, m, K, flush=True)
    return out


def K0(lam0, L=32, dN=1e-4):
    return (rho_vac(L, lam0, 1 + dN) - rho_vac(L, lam0, 1 - dN)) / (2 * dN)


# ---------------------------------------------------------------- hub source (flat space, Woodbury)
def hub_source(L, lam0, dl, hubsites, half=6, nq=160):
    """Energy-density change S_v = delta rho_v of a hub (on-site dl*J on hubsites) at all sites within a cube of
    half-width `half`. rho_v = Re sum_{i in v} (h P)_ii."""
    H = A.bloch_h(L, lam0)[0]
    mu, gap, lo, hi = C.mu_gap(L, lam0)
    t, w = np.polynomial.legendre.leggauss(nq); t = 0.5 * (t + 1); w = 0.5 * w
    ys = gap * t / (1 - t); wy = w * gap / (1 - t) ** 2
    hub = np.array(hubsites); D = np.kron(np.eye(len(hub)), dl * J)
    r = np.arange(-half, half + 1)
    sites = np.array([(a, b, c) for a in r for b in r for c in r])
    nbr = [sites + np.eye(3, dtype=int)[a] for a in range(3)] + [sites - np.eye(3, dtype=int)[a] for a in range(3)]
    # modes needed: all modes of each site, and for each site's link mode its partner at the neighbouring site
    dPss = np.zeros((len(sites), M, M), complex); dPlink = np.zeros((len(sites), 6), complex); dE = 0.0
    hub_dP = np.zeros((len(hub) * M, len(hub) * M), complex)
    for y, wt in zip(ys, wy):
        Gr = C.Gr_at(H, mu + 1j * y)
        Ghh = C.block(Gr, hub, hub, L)
        T = D @ np.linalg.solve(np.eye(len(D)) - Ghh @ D, np.eye(len(D)))
        dE += -wt / np.pi * np.log(abs(np.linalg.det(np.eye(len(D)) - Ghh @ D)))
        GSh = C.block(Gr, sites, hub, L); GhS = C.block(Gr, hub, sites, L)
        U = (GSh @ T).reshape(len(sites), M, -1)          # (ns, 6, nh)
        Wm = GhS.reshape(-1, len(sites), M)               # (nh, ns, 6)
        dG = np.einsum("sar,rsb->sab", U, Wm)
        dPss += wt / np.pi * 0.5 * (dG + np.conj(np.swapaxes(dG, 1, 2)))
        for d in range(6):           # link of mode d at site s to partner mode d^1 at site s +/- e_axis
            a = d // 2; sgn = 1 if d % 2 == 0 else -1
            nb = sites + sgn * np.eye(3, dtype=int)[a]
            GNh = C.block(Gr, nb, hub, L).reshape(len(sites), M, -1)[:, d ^ 1, :]
            GhN = C.block(Gr, hub, nb, L).reshape(-1, len(sites), M)[:, :, d ^ 1]
            dG_sn = np.einsum("sr,rs->s", U[:, d, :], GhN)            # dG[(s,d),(nb,d^1)]
            dG_ns = np.einsum("sr,rs->s", GNh @ T, Wm[:, :, d])        # dG[(nb,d^1),(s,d)]
            dPlink[:, d] += wt / np.pi * 0.5 * (dG_sn + np.conj(dG_ns))
        hub_dP += wt / np.pi * 0.5 * ((Ghh @ T @ Ghh) + (Ghh @ T @ Ghh).conj().T)
    # rho_v = Re[ sum_{i in v} sum_j h_ij P_ji ]: on-site part lam0*J (+dl*J on hub), link part -P_{partner,i}
    Pr = A.P_real(*A.bands(L, lam0)[:2]); P0ss = Pr[0, 0, 0]
    S = np.zeros(len(sites))
    Slink = np.real(np.sum(-dPlink, axis=1))           # half of the adjacent link energies at each site (used by EX-GC43)
    hubset = {tuple(h): k for k, h in enumerate(hub)}
    for s in range(len(sites)):
        onsite = lam0 * J
        S[s] = np.real(np.trace(onsite @ dPss[s]))
        S[s] += np.real(np.sum(-dPlink[s]))          # h_{i,partner} = -1, P_{partner,i} = conj(P_{i,partner}); Re equal
        key = tuple(sites[s])
        if key in hubset:
            k = hubset[key]
            Pfull = P0ss + dPss[s]
            S[s] += np.real(np.trace(dl * J @ Pfull))
    return dict(sites=sites, S=S, dE0=float(dE), sumS=float(S.sum()), Slink=Slink)


# ---------------------------------------------------------------- far field
def far_field(S_sites, S, Kfun, fourpiG, Lg=128, variant="A", K0v=0.0):
    src = np.zeros((Lg, Lg, Lg))
    for (a, b, c), v in zip(S_sites, S):
        src[a % Lg, b % Lg, c % Lg] += v
    src -= src.mean()                               # jellium: periodic Poisson needs a zero-mean source
    k = 2 * np.pi * np.fft.fftfreq(Lg)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    qt2 = 2 * (3 - np.cos(KX) - np.cos(KY) - np.cos(KZ))
    Kt = Kfun(qt2) - (K0v if variant == "A" else 0.0)
    den = -(qt2 + fourpiG * Kt)                       # lattice Laplacian eigenvalue is -qt2
    Sq = np.fft.fftn(src)
    nq = np.zeros_like(Sq); m = qt2 > 0
    nq[m] = fourpiG * Sq[m] / den[m]
    n = np.real(np.fft.ifftn(nq))
    # periodic lattice Green's function g: Lap g = delta - 1/V
    gq = np.zeros_like(Sq); gq[m] = 1.0 / (-qt2[m]); g = np.real(np.fft.ifftn(gq))
    return n, g


def fit_A(n, g, rs=range(4, 21)):
    """n = 4 pi (G_eff M) g along the three axes and the (1,1,1) diagonal; returns G_eff M and the residual."""
    vals = []; gs = []
    for r in rs:
        for v in ((r, 0, 0), (0, r, 0), (0, 0, r)):
            vals.append(n[v]); gs.append(g[v])
    vals = np.array(vals); gs = np.array(gs)
    c = np.sum(vals * gs) / np.sum(gs * gs)
    res = float(np.max(np.abs(vals - c * gs) / np.abs(vals)))
    return dict(GeffM=float(c / (4 * np.pi)), max_rel_resid=res)


def fit_B(n, rs=np.arange(4, 41)):
    """r n(r) along +x fitted to a cos(k r + phi) by a scan over k."""
    y = np.array([r * n[r, 0, 0] for r in rs])
    best = None
    for kk in np.linspace(0.01, 1.5, 3000):
        X = np.column_stack([np.cos(kk * rs), np.sin(kk * rs)])
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        r2 = float(np.sum((X @ c - y) ** 2))
        if best is None or r2 < best[1]:
            best = (kk, r2, c)
    return dict(k=float(best[0]), resid=best[1] / float(np.sum(y ** 2)), amp=float(np.hypot(*best[2])))


# ---------------------------------------------------------------- dense self-consistent check (C2)
def dense_H(L, lam0, N, hub_dl=None):
    h, idx = A.real_h(L, 0.0)                       # hopping only
    Nm = np.repeat(N.reshape(-1), M)
    h = h * np.outer(Nm, Nm)                         # link weight N_v N_w (on-site block of real_h is zero at lam = 0)
    for s in range(L ** 3):
        b = s * M; h[b:b + M, b:b + M] += N.reshape(-1)[s] * lam0 * J
    if hub_dl is not None:
        h[0:M, 0:M] += N.reshape(-1)[0] * hub_dl * J
    return h


def dense_rho(h, L):
    e, v = np.linalg.eigh(h); n = h.shape[0] // 2
    dens = (np.abs(v[:, :n]) ** 2) @ e[:n]
    return dens.reshape(L ** 3, M).sum(1).reshape(L, L, L)


def poisson_periodic(src, fourpiG):
    Lg = src.shape[0]
    k = 2 * np.pi * np.fft.fftfreq(Lg)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    qt2 = 2 * (3 - np.cos(KX) - np.cos(KY) - np.cos(KZ))
    Sq = np.fft.fftn(src - src.mean()); nq = np.zeros_like(Sq); m = qt2 > 0
    nq[m] = fourpiG * Sq[m] / (-qt2[m])
    return np.real(np.fft.ifftn(nq))


def check_c2(L=8, lam0=0.1, dl=0.2, fourpiG=0.01, iters=40, eps=1e-4):
    one = np.ones((L, L, L))
    rv = {Nv: rho_vac(L, lam0, Nv) for Nv in (0.98, 0.99, 1.0, 1.01, 1.02)}
    cfit = np.polyfit(list(rv), list(rv.values()), 2)
    rvac = lambda Nf: np.polyval(cfit, Nf)
    rho_flat = dense_rho(dense_H(L, lam0, one), L)
    rho_hub = dense_rho(dense_H(L, lam0, one, dl), L)
    S_hub = rho_hub - rho_flat
    # exact kernel on L^3 from a one-site lapse perturbation (variant B source response)
    Np = one.copy(); Np[0, 0, 0] += eps; Nm_ = one.copy(); Nm_[0, 0, 0] -= eps
    Kr = (dense_rho(dense_H(L, lam0, Np), L) - dense_rho(dense_H(L, lam0, Nm_), L)) / (2 * eps)
    Kq = np.real(np.fft.fftn(Kr))
    K0d = float(Kr.sum())
    k = 2 * np.pi * np.fft.fftfreq(L)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    qt2 = 2 * (3 - np.cos(KX) - np.cos(KY) - np.cos(KZ))
    out = dict(K0_dense=K0d, K0_bloch=float(np.polyval(np.polyder(cfit), 1.0)), sumS_hub=float(S_hub.sum()))
    for variant in ("A", "B"):
        Kt = Kq - (K0d if variant == "A" else 0.0)
        Sq = np.fft.fftn(S_hub - S_hub.mean()); nq = np.zeros_like(Sq); m = qt2 > 0
        nq[m] = fourpiG * Sq[m] / (-(qt2[m] + fourpiG * Kt[m]))
        n_lin = np.real(np.fft.ifftn(nq))
        n = np.zeros((L, L, L))
        for it in range(iters):
            N = 1 + n
            rho = dense_rho(dense_H(L, lam0, N, dl), L)
            src = rho - (rvac(N) if variant == "A" else rvac(np.ones_like(N)))
            n_new = poisson_periodic(src, fourpiG)
            dn = float(np.abs(n_new - n).max())
            n = 0.5 * n + 0.5 * n_new
            if dn < 1e-11:
                break
        out[variant] = dict(iters=it + 1, last_change=dn,
                            rows=[dict(r=r, n_sc=float(n[r, 0, 0]), n_lin=float(n_lin[r, 0, 0]),
                                       rel=float(abs(n[r, 0, 0] / n_lin[r, 0, 0] - 1))) for r in (2, 3, 4)])
        print("C2", variant, out[variant], flush=True)
    return out


if __name__ == "__main__":
    part = sys.argv[1]
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if part == "ep":
        res = {}
        for lam0 in (0.1, 0.2):
            res[str(lam0)] = {str(Nv): dict(rho_vac=rho_vac(32, lam0, Nv), hub={str(dl): hub_uniform(16, lam0, Nv, dl) for dl in (0.1, 0.2, 0.4)})
                              for Nv in (0.9, 1.0, 1.1)}
    elif part == "kernel":
        res = {}
        for lam0 in (0.1, 0.2):
            res[str(lam0)] = dict(K0=K0(lam0), axis=kernel_axis(lam0, ms=range(1, 17)), diag=kernel_diag(lam0), qpi_norm="1/L")
    elif part == "source":
        res = {"cases": {}}
        for lam0 in (0.1, 0.2):
            for key, hs in (("1", [(0, 0, 0)]), ("7", C.NEIGH)):
                for dl in (0.1, 0.2, 0.4):
                    r = hub_source(16, lam0, dl, hs)
                    res["sites"] = r["sites"].tolist()   # S stored to 12 significant digits (repo size)
                    res["cases"][f"{lam0}_{key}_{dl}"] = dict(S=[float("%.12e" % x) for x in r["S"]], dE0=r["dE0"], sumS=r["sumS"])
                    print(lam0, key, dl, r["dE0"], r["sumS"], flush=True)
    elif part == "c2":
        res = check_c2()
    elif part == "far":
        kern = json.loads((OUT / "kernel.json").read_text())["result"]
        src = json.loads((OUT / "source.json").read_text())["result"]
        res = {}
        for lam0 in (0.1, 0.2):
            kd = kern[str(lam0)]; K0v = kd["K0"]
            qx = np.array([0.0] + [v["qt2"] for v in kd["axis"].values()])
            # The committed kernel.json predates the q = pi normalisation fix (cos^2 = 1 there, so 1/L not 2/L) and has no
            # "qpi_norm" flag: halve its m = 16 entry. A kernel.json written by the fixed generator carries the flag.
            fixq = kd.get("qpi_norm") != "1/L"
            Kx = np.array([K0v] + [v["K_B"] * (0.5 if (fixq and int(m) == 16) else 1.0) for m, v in kd["axis"].items()])
            o = np.argsort(qx); qx, Kx = qx[o], Kx[o]
            Kfun = lambda q2, qx=qx, Kx=Kx: np.interp(q2, qx, Kx)      # constant beyond the measured range
            small = qx[1:5]; c2 = float(np.polyfit(small, Kx[1:5] - K0v, 1)[0])
            res[str(lam0)] = dict(K0=K0v, c2=c2, cases={})
            for key in ("1", "7"):
                for dl in (0.1, 0.2, 0.4):
                    sc = dict(src["cases"][f"{lam0}_{key}_{dl}"], sites=src["sites"])
                    for fg in (0.01, 0.05):
                        G = fg / (4 * np.pi); row = {}
                        nA, g = far_field(sc["sites"], sc["S"], Kfun, fg, variant="A", K0v=K0v)
                        fa = fit_A(nA, g)
                        row["A"] = dict(fa, ratio=fa["GeffM"] / (G * sc["dE0"]), pred=1 / (1 + fg * c2),
                                        n_axis=[float(nA[r, 0, 0]) for r in range(0, 41)])
                        nB, _ = far_field(sc["sites"], sc["S"], Kfun, fg, variant="B", K0v=K0v)
                        fb = fit_B(nB)
                        kJ = float(np.sqrt(fg * abs(K0v)))
                        fbA = fit_A(nB, g)
                        row["B"] = dict(fb, kJ_pred=kJ, k_ratio=fb["k"] / kJ, oneover_r_resid=fbA["max_rel_resid"],
                                        n_axis=[float(nB[r, 0, 0]) for r in range(0, 64)])
                        res[str(lam0)]["cases"][f"{key}_{dl}_{fg}"] = row
                        print(lam0, key, dl, fg, "A", fa, row["A"]["ratio"], "B k", fb["k"], kJ, flush=True)
    res = dict(result=res, wall_seconds=time.time() - t0, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("done", part, res["wall_seconds"])
