"""EX-GC39A: entanglement equilibrium on the lattice. Does matter set the geometry through the vacuum's entanglement?
Plan: spec/EX-GC39A_plan.md. Our construction, not LaViolette's Model G.

Free scalar on a periodic 16^3 lattice, metric gamma = (1 + 2 psi) delta:
  H = sum_v [pi^2 / (2 sqrt g) + sqrt g m^2 phi^2 / 2] + sum_links (1 + 2 psi_link)^{1/2} (dphi)^2 / 2.
Geometry response eta_eff (fixed-volume entropy change per unit area change), matter response k_B
(first-order entropy change per unit local energy density, from a smooth non-local squeezing), G_eq.
"""
import json, sys, time
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC39A"
L = 16
C = np.array([7.5, 7.5, 7.5])        # centre between sites (blocks of even side are centred here)


def sym_entropy(X, P):
    nu = np.sqrt(np.clip(np.linalg.eigvals(X @ P).real, 0.25, None))
    a, b = nu + 0.5, nu - 0.5
    return float(np.sum(a * np.log(a) - np.where(b > 1e-15, b * np.log(np.where(b > 1e-15, b, 1)), 0)))


def coords():
    g = np.arange(L)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    return np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)


XYZ = coords()
IDX = lambda x, y, z: ((x % L) * L + (y % L)) * L + (z % L)


def psi_profile(eps, centre=C):
    return -eps * np.sum(np.cos(2 * np.pi * (XYZ - centre) / L), 1)


def hamiltonian(psi, m):
    n = L ** 3
    sg = (1 + 2 * psi) ** 1.5
    V = np.zeros((n, n))
    V[np.arange(n), np.arange(n)] += sg * m * m
    for a in range(3):
        e = np.zeros(3, int); e[a] = 1
        j = IDX(*(XYZ + e).T)
        w = (1 + 2 * 0.5 * (psi + psi[j])) ** 0.5
        i = np.arange(n)
        V[i, i] += w; V[j, j] += w
        V[i, j] -= w; V[j, i] -= w
    return 1 / sg, V                          # kinetic coefficient (diag) and potential matrix


def ground(psi, m):
    kin, V = hamiltonian(psi, m)
    Dh = np.sqrt(kin)
    S = Dh[:, None] * V * Dh[None, :]
    lam, Q = np.linalg.eigh(S)
    A = Dh[:, None] * Q                       # X = A diag(1/(2 sqrt lam)) A^T ; P = B diag(sqrt lam / 2) B^T
    B = Q / Dh[:, None]
    return A, B, lam


def block_sites(l, centre=C):
    lo = np.round(centre - (l - 1) / 2).astype(int) if l % 2 else np.floor(centre - (l - 1) / 2 + 1e-9).astype(int)
    g = [np.arange(lo[a], lo[a] + l) for a in range(3)]
    G = np.stack(np.meshgrid(*g, indexing="ij"), -1).reshape(-1, 3)
    return IDX(*G.T), G


def block_XP(gs, sites):
    A, B, lam = gs
    a, b = A[sites], B[sites]
    return (a / (2 * np.sqrt(lam))) @ a.T, (b * (np.sqrt(lam) / 2)) @ b.T


def flat_XP_fft(m):
    k = 2 * np.pi * np.arange(L) / L
    K2 = sum(np.meshgrid(*(4 * np.sin(k / 2) ** 2,) * 3, indexing="ij"))
    E = np.sqrt(K2 + m * m)
    return np.real(np.fft.ifftn(1 / (2 * E))), np.real(np.fft.ifftn(E / 2))


def cube_dA_fixV(psi_fun, l, centre=C, nq=60):
    """2A(<psi>_A - <psi>_V) for the coordinate cube of side l (faces at half-integers) around the block centre."""
    lo = block_sites(l, centre)[1].min(0) - 0.5
    u = (np.arange(nq) + 0.5) / nq * l
    U1, U2 = np.meshgrid(u, u, indexing="ij")
    face = []
    for a in range(3):
        for side in (0.0, l):
            P = np.zeros((nq * nq, 3))
            others = [b for b in range(3) if b != a]
            P[:, a] = lo[a] + side; P[:, others[0]] = lo[others[0]] + U1.ravel(); P[:, others[1]] = lo[others[1]] + U2.ravel()
            face.append(psi_fun(P))
    psiA = float(np.mean(np.concatenate(face)))
    W = np.stack(np.meshgrid(u, u, u, indexing="ij"), -1).reshape(-1, 3) + lo
    psiV = float(np.mean(psi_fun(W)))
    A = 6 * l * l
    return 2 * A * (psiA - psiV), psiV


def squeeze_F(k0, q, centre=C):
    """F(x,y) = [sum_{|k|<=k0} cos k.(x-y) / V] * (cos q.(x-c) + cos q.(y-c)) / 2 : smooth, non-local, symmetric."""
    k = 2 * np.pi * np.fft.fftfreq(L, 1 / L) / L
    KS = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1).reshape(-1, 3)
    KS = KS[np.linalg.norm(KS, axis=1) <= k0 + 1e-12]
    E = np.exp(1j * (XYZ @ KS.T))                    # plane waves, (sites, nk)
    low = np.real(E @ E.conj().T) / L ** 3
    mod = np.cos(XYZ @ q - centre @ q)
    return low * 0.5 * (mod[:, None] + mod[None, :]), len(KS)


def energy_density(dX, dP, m):
    """First-order site energy density of the flat Hamiltonian: pi^2/2 + m^2 phi^2/2 + half of each attached link."""
    n = L ** 3
    i = np.arange(n)
    rho = 0.5 * dP[i, i] + 0.5 * m * m * dX[i, i]
    for a in range(3):
        e = np.zeros(3, int); e[a] = 1
        for sgn in (1, -1):
            j = IDX(*(XYZ + sgn * e).T)
            rho += 0.25 * (dX[i, i] + dX[j, j] - 2 * dX[i, j])
    return rho


def main(m, ells=(2, 3, 4, 5), eps=2e-3, r=1e-3):
    t0 = time.process_time()
    out = dict(m=m, L=L, ells=list(ells), eps=eps, r=r)
    flat = ground(np.zeros(L ** 3), m)
    # G1: flat eigh vs FFT correlators
    xf, pf = flat_XP_fft(m)
    s, _ = block_sites(4)
    Xb, Pb = block_XP(flat, s)
    D = (XYZ[s][:, None, :] - XYZ[s][None, :, :]) % L
    out["G1_max_dev"] = float(max(np.max(np.abs(Xb - xf[D[..., 0], D[..., 1], D[..., 2]])),
                                  np.max(np.abs(Pb - pf[D[..., 0], D[..., 1], D[..., 2]]))))
    S0 = {l: sym_entropy(*block_XP(flat, block_sites(l)[0])) for l in range(1, max(ells) + 2)}
    out["S_flat"] = S0
    # ---- geometry response
    c0 = 1e-3
    def blockS(psi):
        gs = ground(psi, m)
        return {l: sym_entropy(*block_XP(gs, block_sites(l)[0])) for l in ells}
    Sgp, Sgm = blockS(psi_profile(+eps)), blockS(psi_profile(-eps))
    Sup, Sum = blockS(np.full(L ** 3, c0)), blockS(np.full(L ** 3, -c0))
    lap_psi = 3 * eps * (2 * np.pi / L) ** 2          # at the centre, for psi = -eps sum cos
    geo = {}
    for l in ells:
        s, _ = block_sites(l)
        dS = (Sgp[l] - Sgm[l]) / 2
        dSu = (Sup[l] - Sum[l]) / 2 / c0   # per unit uniform psi
        dA, psiV = cube_dA_fixV(lambda P: -eps * np.sum(np.cos(2 * np.pi * (P - C) / L), 1), l)
        dS_fixV = dS - psiV * dSu
        # alternative fixed-volume correction (scale-based, descriptive): - <psi>_V * l dS/dl from the flat S(l)
        ldSdl = l * (S0[l + 1] - S0[l - 1]) / 2
        dS_fixV_scale = dS - psiV * ldSdl
        lin = Sgp[l] - S0[l], S0[l] - Sgm[l]
        geo[l] = dict(dS=dS, dS_uniform_per_psi=dSu, psiV=psiV, dA_fixV=dA, c_B=dA / lap_psi,
                      eta_eff=dS_fixV / dA, eta_eff_scale_corr=dS_fixV_scale / dA,
                      linearity=float(abs(lin[0] - lin[1]) / max(abs(lin[0]), 1e-300)))
        print("geo", l, geo[l], flush=True)
    out["geometry"] = geo
    # ---- matter response: smooth non-local squeezing phi -> e^{-rF} phi, pi -> e^{rF} pi
    q = np.array([0, 0, 2 * np.pi / L])
    F, nk = squeeze_F(2 * np.pi * 2 / L, q)
    A, B, lam = flat
    Xfull = (A / (2 * np.sqrt(lam))) @ A.T
    Pfull = (B * (np.sqrt(lam) / 2)) @ B.T
    dX1 = -(F @ Xfull + Xfull @ F)
    dP1 = F @ Pfull + Pfull @ F
    rho1 = energy_density(dX1, dP1, m)          # per unit r
    mat = {}
    for l in ells:
        s, _ = block_sites(l)
        Xb, Pb = Xfull[np.ix_(s, s)], Pfull[np.ix_(s, s)]
        dXb, dPb = dX1[np.ix_(s, s)], dP1[np.ix_(s, s)]
        Sp = sym_entropy(Xb + r * dXb, Pb + r * dPb); Sm = sym_entropy(Xb - r * dXb, Pb - r * dPb)
        dS = (Sp - Sm) / (2 * r)
        rho_B = float(np.mean(rho1[s]))
        mat[l] = dict(dS_per_r=dS, rho_block_per_r=rho_B, k_B=dS / rho_B,
                      linearity=float(abs((Sp - S0[l]) + (Sm - S0[l])) / abs(Sp - Sm)))
        print("mat", l, mat[l], flush=True)
    out["matter"] = mat
    out["matter_modes"] = nk
    out["rho_total_per_r"] = float(np.sum(rho1))
    # ---- equilibrium G
    eta38 = {0.3: 0.022005057765912824, 0.5: 0.019608694421226883}[m]
    eq = {}
    for l in ells:
        g, mt = geo[l], mat[l]
        eq[l] = dict(G_eq=mt["k_B"] / (4 * np.pi * g["eta_eff"] * g["c_B"]),
                     G_eq_scale_corr=mt["k_B"] / (4 * np.pi * g["eta_eff_scale_corr"] * g["c_B"]),
                     G_ideal_with_eta38=mt["k_B"] / (4 * np.pi * eta38 * g["c_B"]),
                     eta_eff_over_eta38=g["eta_eff"] / eta38, eta_eff_scale_over_eta38=g["eta_eff_scale_corr"] / eta38)
    out["equilibrium"] = eq
    out["one_over_4eta38"] = 1 / (4 * eta38)
    out["eta38"] = eta38
    out["cpu_seconds"] = time.process_time() - t0
    return out, (flat, dX1, dP1)


def g3(m, gs_flat, dX1, dP1, ells=(2, 3), prec=600):
    """Separate path: delta<K_B> = 1/2 Tr(M dX) + 1/2 Tr(N dP) with K_B from arbitrary-precision flat correlators."""
    import flint
    from flint import arb, arb_mat
    import exgc38b as B38
    flint.ctx.prec = prec
    k = [2 * np.pi * j / L for j in range(L)]
    out = {}
    for l in ells:
        s, Gs = block_sites(l)
        disp = {}
        def corr(d):
            key = tuple(int(v) % L for v in d)
            if key not in disp:
                sx, sp = arb(0), arb(0)
                pi = arb.pi()
                for a in range(L):
                    for b in range(L):
                        for c in range(L):
                            K2 = sum(4 * (arb(2 * j) * pi / L / 2).sin() ** 2 for j in (a, b, c)) + arb(m) ** 2
                            E = K2.sqrt()
                            ph = (2 * pi * (a * key[0] + b * key[1] + c * key[2]) / L).cos()
                            sx += ph / (2 * E); sp += ph * E / 2
                disp[key] = (sx / L ** 3, sp / L ** 3)
            return disp[key]
        n = len(s)
        X = arb_mat(n, n); P = arb_mat(n, n)
        for i in range(n):
            for j in range(n):
                x, p = corr(Gs[i] - Gs[j]); X[i, j] = x; P[i, j] = p
        Nk, Mk, eps, nu = B38.ent_hamiltonian(X, P)
        Mf = np.array([[float(Mk[i, j].mid()) for j in range(n)] for i in range(n)])
        Nf = np.array([[float(Nk[i, j].mid()) for j in range(n)] for i in range(n)])
        dXb, dPb = dX1[np.ix_(s, s)], dP1[np.ix_(s, s)]
        out[l] = dict(dK=float(0.5 * np.sum(Mf * dXb) + 0.5 * np.sum(Nf * dPb)))
    return out


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    m = float(sys.argv[1])
    res, (flat, dX1, dP1) = main(m)
    if m == 0.3:
        res["G3_first_law"] = g3(m, flat, dX1, dP1)
        for l, v in res["G3_first_law"].items():
            v["dS_per_r"] = res["matter"][l]["dS_per_r"]
            v["rel"] = abs(v["dK"] - v["dS_per_r"]) / abs(v["dS_per_r"])
    (OUT / f"m{m:g}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k in ("equilibrium", "G1_max_dev", "G3_first_law")}, indent=1, default=float))
