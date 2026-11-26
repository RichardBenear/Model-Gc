"""EX-GC37A S2 (static scalar sector, amendment A2): induced static kernel for (n, h_T, h_L) from one lattice scalar on Y1.
Exact vacuum energy from normal modes; block-diagonal in (k_x, k_y) for fields depending on z only.
Our construction, not LaViolette's Model G."""
import json, sys, time
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC37A"
K = lambda k: 2 * np.sin(k / 2)


def coeffs(n, hT, hL, m):
    """Per-site coefficients (arrays over z) for N = 1+n, gamma = diag(1+hT, 1+hT, 1+hL)."""
    N = 1 + n
    sg = (1 + hT) * np.sqrt(1 + hL)
    A = N / sg                                   # kinetic (inverse mass)
    Mm = N * sg * m * m
    cx = sg / (1 + hT); cz = sg / (1 + hL)
    Cx = N * (cx + 0.5 * hT); Rx = -0.25 * N * hT   # same for y
    Cz = N * (cz + 0.5 * hL); Rz = -0.25 * N * hL
    return A, Mm, Cx, Rx, Cz, Rz


def E0_blocks(nz, nxy, fields, m):
    n, hT, hL = fields
    A, Mm, Cx, Rx, Cz, Rz = coeffs(n, hT, hL, m)
    # z-chain: g_v = psi_{v+1} - psi_v ; E_z = 1/2 g^T W g, W_vv = (Cz_v + Cz_{v+1})/2, W_{v,v-1} = Rz_v
    D = -np.eye(nz) + np.roll(np.eye(nz), 1, axis=1)
    W = np.diag(0.5 * (Cz + np.roll(Cz, -1)))
    for v in range(nz):
        W[v, (v - 1) % nz] += Rz[v]; W[(v - 1) % nz, v] += Rz[v]
    Mz = D.T @ W @ D
    k = 2 * np.pi * np.arange(nxy) / nxy
    KX, KY = np.meshgrid(k, k, indexing="ij")
    kx2, ky2 = K(KX).ravel() ** 2, K(KY).ravel() ** 2
    cxs, cys = np.cos(KX).ravel(), np.cos(KY).ravel()
    diag = (Mm[None, :] + Cx[None, :] * (kx2 + ky2)[:, None]
            + 2 * Rx[None, :] * (kx2 * cxs + ky2 * cys)[:, None])            # (blocks, nz)
    M = Mz[None, :, :] + diag[:, :, None] * np.eye(nz)[None]
    sa = np.sqrt(A)
    H = sa[None, :, None] * M * sa[None, None, :]
    ev = np.linalg.eigvalsh(H)
    return 0.5 * np.sum(np.sqrt(np.maximum(ev, 0)))


def E0_realspace(L, fields, m):
    """Separate path: full real-space L^3 Hamiltonian, fields depend on z, built site by site."""
    n, hT, hL = fields
    A, Mm, Cx, Rx, Cz, Rz = coeffs(n, hT, hL, m)
    idx = lambda x, y, z: ((x % L) * L + (y % L)) * L + (z % L)
    Vn = L ** 3
    Mt = np.zeros((Vn, Vn)); Av = np.zeros(Vn)
    def addgrad(c, a_, b_):   # adds c * g_a g_b where g = phi[a_[1]] - phi[a_[0]] etc. into 1/2 phi^T Mt phi (Mt += 2c sym)
        (p0, p1), (q0, q1) = a_, b_
        for (i, si) in ((p1, 1), (p0, -1)):
            for (j, sj) in ((q1, 1), (q0, -1)):
                Mt[i, j] += c * si * sj; Mt[j, i] += c * si * sj
    for x in range(L):
        for y in range(L):
            for z in range(L):
                v = idx(x, y, z); Av[v] = A[z]; Mt[v, v] += Mm[z]
                for (C, R, e) in ((Cx[z], Rx[z], (1, 0, 0)), (Cx[z], Rx[z], (0, 1, 0)), (Cz[z], Rz[z], (0, 0, 1))):
                    fwd = (v, idx(x + e[0], y + e[1], z + e[2])); bwd = (idx(x - e[0], y - e[1], z - e[2]), v)
                    addgrad(0.25 * C, fwd, fwd); addgrad(0.25 * C, bwd, bwd)   # 1/2 C * 1/2 (g_f^2 + g_b^2)
                    addgrad(R, fwd, bwd)                                          # R g_f g_b
    sa = np.sqrt(Av)
    ev = np.linalg.eigvalsh(sa[:, None] * Mt * sa[None, :])
    return 0.5 * np.sum(np.sqrt(np.maximum(ev, 0)))


def kernel(nz, nxy, j, m, eps=2e-3):
    z = np.arange(nz); prof = np.cos(2 * np.pi * j * z / nz)
    base = np.zeros(nz)
    def E(a):
        return E0_blocks(nz, nxy, [base + a[0] * prof, base + a[1] * prof, base + a[2] * prof], m)
    e0 = E((0, 0, 0)); V = nz * nxy * nxy
    P = np.zeros((3, 3))
    for i in range(3):
        u = np.zeros(3); u[i] = eps
        P[i, i] = 4 * (E(u) + E(-u) - 2 * e0) / (2 * eps ** 2) / V
    for i in range(3):
        for k in range(i + 1, 3):
            u = np.zeros(3); u[i] = eps; w = np.zeros(3); w[k] = eps
            c = E(u + w) + E(-u - w) - E(u - w) - E(-u + w)
            P[i, k] = P[k, i] = c / (2 * V * eps ** 2)
    if j == 0:
        P *= 0.5   # cos(0) = 1: <cos^2> = 1, not 1/2
    return P


def eh_control():
    """gamma = Pi_Tn / (2 Pi_TT) applied to the 34N Y1 EH kernel (static, q along z)."""
    import exgc34n as N34
    out = []
    for qz in (0.05, 0.2):
        M = N34.kernel(np.array([0, 0, qz]), 0.0).real
        MTn = M[4, 0] + M[5, 0]; MTT = M[4, 4] + M[4, 5] + M[5, 4] + M[5, 5]
        out.append(dict(qz=qz, gamma=float(MTn / (2 * MTT)), M_nn=float(M[0, 0]), M_nL=float(M[0, 6]), M_TL=float(M[4, 6] + M[5, 6]), M_LL=float(M[6, 6])))
    return out


if __name__ == "__main__":
    t0 = time.process_time()
    res = {}
    if sys.argv[1] == "check":
        L, m = 6, 0.5
        z = np.arange(L); prof = 0.03 * np.cos(2 * np.pi * z / L)
        f = [prof, -0.7 * prof, 0.4 * prof]
        a, b = E0_blocks(L, L, f, m), E0_realspace(L, f, m)
        res = dict(L=L, blocks=a, realspace=b, rel=abs(a - b) / abs(a), eh_control=eh_control())
    else:
        m = 0.5
        for nz, nxy in ((32, 32), (48, 48)):
            rows = []
            for j in range(0, 5):
                P = kernel(nz, nxy, j, m)
                rows.append(dict(j=j, q=2 * np.pi * j / nz, Kq2=float(K(2 * np.pi * j / nz) ** 2), Pi=P.tolist()))
                print(nz, j, np.round(P, 7).tolist(), flush=True)
            # q^2 (K^2) coefficients from j = 1..4 relative to j = 0: fit dPi = b K^2 + c K^4
            Kq2 = np.array([r["Kq2"] for r in rows[1:]])
            dP = np.array([np.array(r["Pi"]) - np.array(rows[0]["Pi"]) for r in rows[1:]])
            X = np.column_stack([Kq2, Kq2 ** 2])
            B = np.zeros((3, 3))
            for i in range(3):
                for k in range(3):
                    B[i, k] = np.linalg.lstsq(X, dP[:, i, k], rcond=None)[0][0]
            names = ["n", "T", "L"]
            res[f"N{nz}"] = dict(rows=rows, b_q2={f"{names[i]}{names[k]}": float(B[i, k]) for i in range(3) for k in range(i, 3)},
                                 gamma=float(B[1, 0] / (2 * B[1, 1])), nn_over_Tn=float(B[0, 0] / B[1, 0]),
                                 ward_LL_over_TT=float(B[2, 2] / B[1, 1]), ward_TL_over_TT=float(B[1, 2] / B[1, 1]),
                                 ward_nL_over_Tn=float(B[0, 2] / B[1, 0]))
            print(nz, {k: v for k, v in res[f"N{nz}"].items() if k != "rows"}, flush=True)
    res["cpu_seconds"] = time.process_time() - t0
    (OUT / f"s2_{sys.argv[1]}.json").write_text(json.dumps(res, indent=2) + "\n")
    print({k: v for k, v in res.items() if not k.startswith("N")})
